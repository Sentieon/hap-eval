# Copyright (c) 2022 Sentieon Inc. All rights reserved
from __future__ import print_function
import argparse
import bisect
import collections
from collections.abc import Callable
import functools
import heapq
import io
import itertools
import operator
import os
import re
import sys
from typing import Generator, Iterable, List, Tuple
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
    '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Difference in length'
    ' between REF and ALT alleles">',
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

def extend_event(
    event_size: int,
    base_extend: int = 20,
    extend_frac: float = 0.2,
    max_extend: int = 10000,
) -> int:
    '''Extension around events for grouping'''
    return min(event_size * extend_frac + base_extend, max_extend)

class VariantGroup(object):
    '''A group of variants'''
    def __init__(self, extend_func: Callable[[int], int]):
        self.data: List[Tuple[int, vcflib.Variant]] = []
        self.end = -1
        self.end_extend = -1
        self.extend_func = extend_func

    def __iter__(self):
        return iter(self.data)

    def __bool__(self):
        return bool(self.data)

    def append(self, d: Tuple[int, vcflib.Variant]):
        self.data.append(d)
        self.end = max(self.end, d[1].end)
        variant_size = max([len(x) for x in [d[1].ref] + d[1].alt])
        extend_dist = self.extend_func(variant_size)
        self.end_extend = max(self.end_extend, d[1].end + extend_dist)

    def tr_extend(self, _tr_s: int, tr_e: int):
        '''Extend the group across a tandem repeat'''
        self.end = max(self.end, tr_e)

class VariantGrouper(object):
    '''Group variants from multiple VCFs by distance'''
    def __init__(
        self,
        vcfs: Iterable[vcflib.VCF],
        contig: str,
        start: int,
        end: int,
        pred: Callable[[vcflib.Variant], bool],
        extend_func: Callable[[int], int] = functools.partial(extend_event, base_extend=20, extend_frac=0.2, max_extend=10000),
        tr_bed: IntervalList | None = None,
    ):
        self.vcfs = vcfs
        self.contig = contig
        self.start = start
        self.end = end
        self.pred = pred
        self.tr_bed = tr_bed
        self.extend_func = extend_func

        self.maxdist = extend_func(999999999999)

    def cluster(self) -> Generator[VariantGroup, None, None]:
        '''Cluster variants in the VCF'''
        q = []
        for k, vcf in enumerate(self.vcfs):
            if vcf is None:
                continue
            i = iter(vcf.range(self.contig, self.start-self.maxdist))
            v = next(i, None)
            while v and not self.pred(v):
                v = next(i, None)
            if v:
                heapq.heappush(q, (v.pos, v.end, k, v, i))

        p1: List[Tuple[int, int, int, vcflib.Variant]] = []
        p2: List[Tuple[int, int, int, vcflib.Variant]] = []
        g = VariantGroup(self.extend_func)
        while q or p1:
            if p1:
                pos, end, k, v = p1.pop(0)
            else:
                pos, end, k, v, i = heapq.heappop(q)
                v2 = next(i, None)
                while v2 and not self.pred(v2):
                    v2 = next(i, None)
                if v2:
                    heapq.heappush(q, (v2.pos, v2.end, k, v2, i))

            if not g:
                if pos >= self.end:
                    break
                else:
                    g.append((k, v))
                    # extend the group if overlapping a tandem repeat
                    if self.tr_bed:
                        tr_regions = self.tr_bed.get(self.contig, g.data[0][1].pos, g.end_extend)
                        for tr_region in tr_regions:
                            g.tr_extend(*tr_region)
                    continue

            variant_size = max([len(x) for x in [v.ref] + v.alt])
            extend_dist = self.extend_func(variant_size)
            if pos >= g.end + self.maxdist:
                # This variant is too far from the group
                if g and g.data[0][1].pos >= self.start and g.data[0][1].pos < self.end:
                    yield g
                g = VariantGroup(self.extend_func)
                p1 = p2 + p1
                p1.append((pos, end, k, v))
                p2 = []
            elif pos <= g.end_extend or g.end >= (pos - extend_dist):
                # This variant is close enough to add to the group after the variants in the queue
                for _, _, k2, v2 in p2:
                    g.append((k2, v2))
                p2 = []
                g.append((k, v))

                # extend the group if overlapping a tandem repeat
                if self.tr_bed:
                    tr_regions = self.tr_bed.get(self.contig, g.data[0][1].pos, g.end_extend)
                    for tr_region in tr_regions:
                        g.tr_extend(*tr_region)
            else:
                # This variant is an intermediate distance from the group, add to a queue
                p2.append((pos, end, k, v))
        if g and g.data[0][1].pos >= self.start and g.data[0][1].pos < self.end:
            yield g

class VCFEvaluator(vcflib.Shardable, vcflib.ShardResult):
    params = ( # name, defval, descr, kwargs
        ('minsize',     50,     'Minimum size of variants to consider', {}),
        ('maxdiff',     0.2,    'Haplotype difference theshold', {}),
        ('metric',      'Levenshtein',
            'Distance metric', {"choices": ['Levenshtein', 'Length']}),
    )

    def __init__(self, ref, vcf1, vcf2, bed, tr_bed, args):
        self.ref = Reference(ref)
        self.vcf1 = vcflib.VCF(vcf1)
        self.vcf2 = vcflib.VCF(vcf2)
        self.bed = bed and IntervalList(bed) or None
        self.tr_bed = tr_bed and IntervalList(tr_bed) or None
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
            # We represent DUPs as left-aligned insertions for a consistent
            # representation with simple insertion alleles
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
        if svtype in ['BND','INV']:
            return False
        if svtype is not None:
            return self.fixup(v)
        return True

    @staticmethod
    def reassemble(c, s, e, r, g, ref):
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
            h = [r] * len(b[0])  # multiple haplotypes for diploid or polyploid samples
            for j,p in enumerate(zip(*b)):  # j = haplotype index
                t = s
                for i in range(len(p)):  # each variant allele
                    if p[i] == 0:
                        continue
                    v = g[i]
                    a = v.alt[p[i]-1]
                    if a == "*":  # Skip upstream deletion alleles
                        continue
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
                        a = ref.get(c, v.pos, v.pos+1+n)
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


    @staticmethod
    def evaluate_cluster(
        g, contig, ref, diff, bed=None, maxdiff=0.2, minsize=50, log=sys.stdout
    ):
        c, s, e = contig, g.data[0][1].pos, g.end
        if bed and not bed.get(c, s, e):
            return None
        g0 = [x[1] for x in g if x[0] == 0]
        g1 = [x[1] for x in g if x[0] == 1]
        if (
            VCFEvaluator.maxsize(g0) < minsize
            and VCFEvaluator.maxsize(g1) < minsize
        ):
            return None
        s = max(s-100, 0)
        r = ref.get(c, s, e+100)
        e = s + len(r)
        desc = [len(g0), len(g1)]
        if g0 and g1:
            best = (999, None, None, None)
            for h0 in VCFEvaluator.reassemble(c, s, e, r, g0, ref):
                for h1 in VCFEvaluator.reassemble(c, s, e, r, g1, ref):
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
            best = (999999999, None)
            for h0 in VCFEvaluator.reassemble(c, s, e, r, g0, ref):
                dd = max(abs(len(h)-len(r)) for h in h0)
                if best[0] > dd:
                    best = (dd, h0)
            dd, h0 = best
            if dd < minsize:
                return None
            if h0:
                desc.extend((len(r),
                    '('+','.join(str(len(h)-len(r)) for h in h0)+')',
                    '(0,0)',
                    1.))
        elif g1:
            call = 'FP'
            best = (999999999, None)
            for h1 in VCFEvaluator.reassemble(c, s, e, r, g1, ref):
                dd = max(abs(len(h)-len(r)) for h in h1)
                if best[0] > dd:
                    best = (dd, h1)
            dd, h1 = best
            if dd < minsize:
                return None
            if h1:
                desc.extend((len(r),
                    '(0,0)',
                    '('+','.join(str(len(h)-len(r)) for h in h1)+')',
                    1.))
        print(call, '%s:%d-%d' % (c,s+1,e), *map(str, desc), file=log)
        return call


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
        grouper = VariantGrouper(vcfs, contig, start, end, pred, tr_bed=self.tr_bed)
        for g in grouper.cluster():
            call = self.evaluate_cluster(
                g,
                contig,
                self.ref,
                diff,
                bed=self.bed,
                maxdiff=self.args.maxdiff,
                minsize=self.args.minsize,
                log=log,
            )
            if not call:
                continue

            g0 = [x[1] for x in g if x[0] == 0]
            g1 = [x[1] for x in g if x[0] == 1]
            call_region = "{}:{}-{}".format(contig, g.data[0][1].pos, g.end)
            for out_vcf, g, err in zip(
                (base_out, comp_out), (g0, g1), ('FN', 'FP')
            ):
                if out_vcf:
                    for v in g:
                        v.line = None  # Update the string representation
                        v.info['CALL'] = 'TP' if call == 'MM' else err
                        v.info['CALL_REGION'] = call_region
                        if any(a == '<DUP>' for a in v.alt):
                            v.end = v.info.get('END', v.pos + len(v.ref))
                        out_vcf.emit(v)

            summary[call] += 1
        return summary

class MixedHelpFormatter(argparse.HelpFormatter):
    def _metavar_formatter(self, action, default):
        if action.metavar is None and action.type is not None:
            action.metavar = action.type.__name__.upper()
        return argparse.HelpFormatter._metavar_formatter(
            self, action, default)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(formatter_class=MixedHelpFormatter)
    parser.add_argument('-r', '--reference', required=True, metavar='FASTA',
        help='Reference file', dest='ref')
    parser.add_argument('-b', '--base', required=True, metavar='VCF',
        help='Baseline vcf file')
    parser.add_argument('-c', '--comp', required=True, metavar='VCF',
        help='Comparison vcf file')
    parser.add_argument('-i', '--interval', required=False, metavar='BED',
        help='Evaluation region file', dest='bed')
    parser.add_argument('--tr_bed', metavar='BED',
        help='Tandem repeat BED file')
    parser.add_argument('-t', '--thread_count', required=False, type=int,
        help='Number of threads', dest='nthr')
    parser.add_argument('--base_out', metavar='VCF',
        help='Annotated baseline vcf file')
    parser.add_argument('--comp_out', metavar='VCF',
        help='Annotated comparison vcf file')
    parser.add_argument('--step_size', default=100*1000*1000, type=int,
        help=argparse.SUPPRESS, dest='step')
    VCFEvaluator.add_arguments(parser)
    return parser.parse_args(argv)


def main():
    args = parse_args()
    eval = VCFEvaluator(
        args.ref, args.base, args.comp, args.bed, args.tr_bed, args
    )

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
