# Copyright (c) 2022 Sentieon Inc. All rights reserved
from __future__ import print_function
import argparse
import bisect
import collections
import functools
import heapq
import io
import itertools
import operator
import os
import re
import sys
import tempfile
import vcflib

from vcflib import bgzf
from vcflib.compat import *

try:
    import Levenshtein
except:
    pass

Contig = collections.namedtuple('Contig', 'length offset width skip')

HAP_EVAL_HDRS = [
    '##INFO=<ID=CALL,Number=1,Type=String,Description="Evaluation call of the'
    ' comparison VCF relative to the baseline VCF">',
    '##INFO=<ID=CALL_REGION,Number=1,Type=String,Description="Evaluation region'
    ' for the variant">',
]


class Reference(object):
    def __init__(self, path):
        self.path = path
        self.index = collections.OrderedDict()
        with io.open(self.path+'.fai', 'rb') as fp:
            for line in fp:
                flds = line.rstrip().decode().split('\t')
                if len(flds) != 5:
                    raise RuntimeError('Corrupted index')
                self.index[flds[0]] = Contig._make(map(int,flds[1:]))
        self.fp = io.open(path, 'rb')

    def __getstate__(self):
        odict = self.__dict__.copy()
        odict['fp'] = self.fp.tell()
        return odict

    def __setstate__(self, ndict):
        path = ndict['path']
        fp = io.open(path, 'rb')
        fp.seek(ndict.pop('fp'))
        ndict['fp'] = fp
        self.__dict__.update(ndict)

    def get(self, c, s, e):
        ci = self.index[c]
        e = min(e, ci.length)
        if s > e:
            raise ValueError
        seq = b''
        while s < e:
            o = ci.offset + s // ci.width * ci.skip + s % ci.width
            n = ci.width - s % ci.width
            n = min(n, e-s)
            self.fp.seek(o)
            seq += self.fp.read(n)
            s += n
        return seq.decode()

    def __iter__(self):
        return iter(self.index.items())

class IntervalList(object):
    def __init__(self, path):
        self.regions = self.load_bed(path)

    def load_bed(self, path):
        regions = {}
        if path.endswith('.gz'):
            fp = bgzf.open(path, 'rb')
        else:
            fp = io.open(path, 'rb')
        for line in fp:
            if line.startswith(b'#') or line.startswith(b'track'):
                continue
            cols = line.rstrip().decode().split('\t')
            if len(cols) < 3:
                continue
            chrom,start,end = cols[:3]
            regions.setdefault(chrom,[]).append([int(start),int(end)])
        fp.close()
        for chrom in regions:
            v = []
            s0 = e0 = None
            for (s,e) in sorted(regions[chrom]):
                if e0 is None:
                    s0 = s
                    e0 = e
                elif s > e0:
                    v.extend((s0,e0))
                    s0 = s
                    e0 = e
                else:
                    e0 = e
            if e0 is not None:
                v.extend((s0,e0))
            regions[chrom] = v
        return regions

    def get(self, c, s, e):
        r = self.regions.get(c)
        v = []
        if r is not None:
            i = bisect.bisect_right(r, s)
            if i % 2:
                v.append((r[i-1], r[i]))
                i += 1
            while i < len(r) and r[i] < e:
                v.append((r[i], r[i+1]))
                i += 2
        return v

    def __contains__(self, cs):
        if isinstance(cs, basestring):
            return cs in self.regions
        c, s = cs
        r = self.regions.get(c)
        return r and bisect.bisect_right(r, s) % 2 != 0

class VCFEvaluator(vcflib.Shardable, vcflib.ShardResult):
    params = ( # name, defval, descr, kwargs
        ('maxdist',     1000,   'Maximum distance to cluster variants', {}),
        ('minsize',     50,     'Minimum size of variants to consider', {}),
        ('maxdiff',     0.2,    'Haplotype difference theshold', {}),
        ('metric',      'Levenshtein',
            'Distance metric', {"choices": ['Levenshtein', 'Length']}),
    )

    def __init__(self, ref, vcf1, vcf2, bed, args):
        self.ref = Reference(ref)
        self.vcf1 = vcflib.VCF(vcf1)
        self.vcf2 = vcflib.VCF(vcf2)
        self.bed = bed and IntervalList(bed) or None
        self.args = args
        self.log = None

    def __shard__(self, cse):
        self.shard = cse
        self.log = tempfile.NamedTemporaryFile(mode='w', suffix='.txt',
            delete=False, dir=os.getenv('SENTIEON_TMPDIR'))
        return self

    def __getdata__(self):
        self.log.close()
        return self.log.name

    def __accum__(self, tlog):
        log = self.log or sys.stdout
        with open(tlog) as fp:
            for line in fp:
                log.write(line)
        os.unlink(tlog)

    @property
    def contigs(self):
        contigs = [(c,0,t.length) for c,t in self.ref]
        if self.bed:
            contigs = [(c,s,e) for c,s,e in contigs if c in self.bed]
        return contigs

    @classmethod
    def add_arguments(cls, parser):
        for k,v,h,kwargs in cls.params:
            h += ' (default: %(default)s)'
            parser.add_argument('--'+k, default=v, type=type(v), help=h, **kwargs)

    @staticmethod
    def cluster(vcfs, c, s, e, maxdist, pred):
        q = []
        for k, vcf in enumerate(vcfs):
            if vcf is None:
                continue
            i = iter(vcf.range(c, s-maxdist))
            v = next(i, None)
            while v and not pred(v):
                v = next(i, None)
            if v:
                heapq.heappush(q, (v.pos, v.end, k, v, i))
        g = []
        while q:
            pos, end, k, v, i = heapq.heappop(q)
            if g and pos >= g[-1][1].end + maxdist:
                if g[0][1].pos >= s and g[0][1].pos < e:
                    yield g
                g = []
                if pos >= e:
                    break
            g.append((k, v))
            v = next(i, None)
            while v and not pred(v):
                v = next(i, None)
            if v:
                heapq.heappush(q, (v.pos, v.end, k, v, i))
        if g and g[0][1].pos >= s and g[0][1].pos < e:
            yield g

    @staticmethod
    def maxsize(g):
        sz = [0, 0]
        for v in g:
            svlen = v.info.get('SVLEN') or [len(a)-len(v.ref) for a in v.alt]
            sz[0] += min(svlen)
            sz[1] += max(svlen)
        return max(-sz[0], sz[1])

    @staticmethod
    def fixup(v):
        # correct various reporting problems
        if any(a == '<DUP>' for a in v.alt):
            # many tools make the mistake of setting END=POS+SVLEN!
            v.end = v.pos + len(v.ref)
        svlen = v.info.get('SVLEN')
        if not isinstance(svlen, list):
            if isinstance(svlen, int):
                svlen = [svlen] * len(v.alt)
            elif any(a[0] == '<' for a in v.alt):
                return False        # symbolic yet w/o SVLEN
            else:
                svlen = [len(a)-len(v.ref) for a in v.alt]
            v.info['SVLEN'] = svlen
        return True

    def filter(self, v):
        if len(v.filter) > 0 and v.filter[0] != 'PASS':
            return False
        if '.' in v.samples[0]['GT']:
            return False
        svtype = v.info.get('SVTYPE')
        if svtype == 'BND':
            return False
        if svtype is not None:
            return self.fixup(v)
        return True

    def reassemble(self, c, s, e, r, g):
        combos = []
        phased = 0
        for v in g:
            gt = re.split(r'([/|])', v.samples[0]['GT'])
            ia = tuple(int(a) for i,a in enumerate(gt) if i % 2 == 0)
            if all(a == ia[0] for a in ia):
                combos.append([ia])
            elif gt[1] == '/':
                combos.append(list(itertools.permutations(ia)))
            elif gt[1] == '|':
                combos.append([ia])
                phased += 1
            else:
                combos.append([])       # invalid
        # if none is phased, fix the first het
        if phased == 0:
            for b in combos:
                if len(b) > 1:
                    del b[1:]
                    break
        #print(combos)
        n = 1
        for b in combos:
            n *= len(b)
        if n > 1024:
            print('%s:%d-%d too many phasing combos %d' %
                (c, s, e, n), file=sys.stderr)
            sys.stderr.flush()
            return
        for b in itertools.product(*combos):
            if not b:
                continue
            h = [r] * len(b[0])
            for j,p in enumerate(zip(*b)):
                t = s
                for i in range(len(p)):
                    if p[i] == 0:
                        continue
                    v = g[i]
                    a = v.alt[p[i]-1]
                    if v.pos+1 == t and v.ref[0] == a[0]:
                        # deletion followed by insertion
                        pass
                    elif v.pos < t:
                        t = -1          # incompatible
                        break
                    if a == '<DEL>':
                        a = v.ref[0]    # imprecise
                        t = v.end
                    elif a == '<DUP>':
                        n = v.info['SVLEN'][p[i]-1]
                        a = self.ref.get(c, v.pos, v.pos+1+n)
                        t = v.pos+1
                    else:
                        t = v.pos+1 + max(len(v.ref)-len(a),0)
                    h[j] = h[j][:v.pos-e] + a + h[j][v.end-e:]
                if t < 0:
                    h = None
                    break
            if h:
                yield h

    @staticmethod
    def diff_by_len(a, b):
        la, lb = len(a), len(b)
        return abs(la-lb)/float(max(la,lb))

    @staticmethod
    def diff_by_lev(a, b):
        return 1. - Levenshtein.ratio(a, b)

    def evaluate(
        self, base_out, comp_out, contig=None, start=0, end=0x7fffffff
    ):
        if contig is None:
            contig, start, end = getattr(self, 'shard')
        log = self.log or sys.stdout
        summary = collections.Counter()
        vcfs = (self.vcf1, self.vcf2)
        pred = self.filter
        diff = self.diff_by_len
        if self.args.metric == 'Levenshtein':
            if 'Levenshtein' in sys.modules:
                diff = self.diff_by_lev
            else:
                print(
                    'Levenshtein module not available. Falling back to length'
                    ' comparison.',
                    file=sys.stderr
                )
        maxdist = self.args.maxdist
        maxdiff = self.args.maxdiff
        minsize = self.args.minsize
        for g in self.cluster(vcfs, contig, start, end, maxdist, pred):
            c, s, e = contig, g[0][1].pos, g[-1][1].end
            call_region = "{}:{}-{}".format(c, s, e)
            if self.bed and not self.bed.get(c, s, e):
                continue
            g0 = [x[1] for x in g if x[0] == 0]
            g1 = [x[1] for x in g if x[0] == 1]
            if self.maxsize(g0) < minsize and self.maxsize(g1) < minsize:
                continue
            s = max(s-100, 0)
            r = self.ref.get(c, s, e+100)
            e = s + len(r)
            desc = [len(g0), len(g1)]
            if g0 and g1:
                best = (999, None, None, None)
                for h0 in self.reassemble(c, s, e, r, g0):
                    for h1 in self.reassemble(c, s, e, r, g1):
                        for hh in itertools.permutations(h1):
                            d = [diff(a,b) for a,b in zip(h0,hh)]
                            dd = min(d)
                            if best[0] > dd:
                                best = (dd, h0, hh, d)
                dd, h0, h1, d = best
                if d is None:
                    call = '??'
                elif min(d) > maxdiff:
                    call = 'XX'
                elif max(d) > maxdiff:
                    call = 'MX'
                else:
                    call = 'MM'
                if d is not None:
                    desc.extend((len(r),
                        '('+','.join(str(len(h)-len(r)) for h in h0)+')',
                        '('+','.join(str(len(h)-len(r)) for h in h1)+')',
                        '%.3f' % max(d)))
            elif g0:
                call = 'FN'
            elif g1:
                call = 'FP'

            for out_vcf, g, err in zip(
                (base_out, comp_out), (g0, g1), ('FN', 'FP')
            ):
                if out_vcf:
                    for v in g:
                        v.line = None  # Update the string representation
                        v.info['CALL'] = 'TP' if call == 'MM' else err
                        v.info['CALL_REGION'] = call_region
                        out_vcf.emit(v)

            print(call, '%s:%d-%d' % (c,s+1,e), *map(str, desc), file=log)
            summary[call] += 1
        return summary

class MixedHelpFormatter(argparse.HelpFormatter):
    def _metavar_formatter(self, action, default):
        if action.metavar is None and action.type is not None:
            action.metavar = action.type.__name__.upper()
        return argparse.HelpFormatter._metavar_formatter(
            self, action, default)

def main():
    parser = argparse.ArgumentParser(formatter_class=MixedHelpFormatter)
    parser.add_argument('-r', '--reference', required=True, metavar='FASTA',
        help='Reference file', dest='ref')
    parser.add_argument('-b', '--base', required=True, metavar='VCF',
        help='Baseline vcf file')
    parser.add_argument('-c', '--comp', required=True, metavar='VCF',
        help='Comparison vcf file')
    parser.add_argument('-i', '--interval', required=False, metavar='BED',
        help='Evaluation region file', dest='bed')
    parser.add_argument('-t', '--thread_count', required=False, type=int,
        help='Number of threads', dest='nthr')
    parser.add_argument('--base_out', metavar='VCF',
        help='Annotated baseline vcf file')
    parser.add_argument('--comp_out', metavar='VCF',
        help='Annotated comparison vcf file')
    parser.add_argument('--step_size', default=100*1000*1000, type=int,
        help=argparse.SUPPRESS, dest='step')
    VCFEvaluator.add_arguments(parser)

    args = parser.parse_args()
    eval = VCFEvaluator(args.ref, args.base, args.comp, args.bed, args)

    out_vcfs = []
    for vcf, out_vcf in zip(
        (eval.vcf1, eval.vcf2), (args.base_out, args.comp_out)
    ):
        if out_vcf:
            out_vcf = vcflib.VCF(out_vcf, 'wb')
            out_vcf.copy_header(
                vcf,
                update=HAP_EVAL_HDRS,
            )
            out_vcf.emit_header()
        out_vcfs.append(out_vcf)

    sharder = vcflib.Sharder(args.nthr)
    shards = sharder.cut(eval.contigs, args.step)
    result = sharder.run(shards, VCFEvaluator.evaluate, [], eval, *out_vcfs)
    summary = sum(result, collections.Counter())
    print(summary)

    for out_vcf in out_vcfs:
        if out_vcf:
            out_vcf.close()

    tp = float(summary['MM'])
    ff = summary['MX'] + summary['XX'] + summary['??']
    fp = summary['FP'] + ff
    fn = summary['FN'] + ff

    if tp < 1:
        print('precision 0.0 recall 0.0 f1 0.0')
        return 0
    print('precision %.4f recall %.4f f1 %.4f' %
        (tp/(tp+fp), tp/(tp+fn), tp/(tp+(fp+fn)/2)))
    return 0

if __name__ == '__main__':
    sys.exit(main())

# vim: ts=4 sw=4 expandtab
