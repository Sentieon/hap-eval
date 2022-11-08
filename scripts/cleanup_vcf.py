#!/usr/bin/env python

# Copyright (c) 2022 Sentieon Inc. All rights reserved

"""
Remove extra INFO and FORMAT fields from input VCFs
"""

import argparse
import itertools
import sys

import vcflib

KEEP_INFOS = ["SVLEN", "SVTYPE", "END"]
KEEP_FORMATS = ["GT"]
HEADER_LINE_FMT = "##{hdr_type}=<ID={fld}>"

# replace vcfilb's `parse_field` with more flexible VCF handling
class VCF(vcflib.VCF):
    @staticmethod
    def parse_field(desc, kv):
        k,v = kv
        d = desc.get(k)
        if d is None:
            return (k, v)
        cvt = VCF.decoders.get(d['Type'], str)
        if v == '.':
            v = None
        elif d['Number'] != '0' and d['Number'] != '1':
            try:
                v = list(map(cvt, v.split(',')))
            except:
                v = None
        elif cvt:
            try:
                v = cvt(v)
            except:
                v = None
        return (k,v)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", required=True, metavar="VCF", help="The input VCF"
    )
    parser.add_argument(
        "--output", required=True, metavar="VCF", help="The output VCF"
    )
    return parser.parse_args(argv)


def main(args):
    vcf = VCF(args.input)
    output_suffixes = [".vcf", ".vcf.gz"]
    if not any([args.output.endswith(suffix) for suffix in output_suffixes]):
        err = (
            "Unrecognized output suffix. Supported suffixes are: "
            + ", ".join(output_suffixes)
        )
        sys.exit(err)
    out_vcf = vcflib.VCF(args.output, "wb")

    rm_infos = [n for n in vcf.infos.keys() if n not in KEEP_INFOS]
    rm_fmt = [n for n in vcf.formats.keys() if n not in KEEP_FORMATS]

    remove = []
    for hdr_type, fld_list in zip(("INFO", "FORMAT"), (rm_infos, rm_fmt)):
        for fld in fld_list:
            remove.append(HEADER_LINE_FMT.format(hdr_type=hdr_type, fld=fld))

    out_vcf.copy_header(
        vcf,
        remove=remove,
    )
    out_vcf.emit_header()

    for variant in vcf:
        variant.line = None
        drop_info = [k for k in variant.info.keys() if k not in out_vcf.infos]
        for k in drop_info:
            del variant.info[k]

        drop_fmt = set(
            k
            for k in itertools.chain(*(s.keys() for s in variant.samples))
            if k not in out_vcf.formats
        )
        for sample in variant.samples:
            for k in drop_fmt:
                del sample[k]

        out_vcf.emit(variant)
    return 0


if __name__ == "__main__":
    args = parse_args()
    sys.exit(main(args))

# vim: ts=4 sw=4 expandtab
