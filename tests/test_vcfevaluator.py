#!/usr/bin/env python

import random
from unittest.mock import Mock

import hap_eval
import vcflib

REGION1_REF = "AACCATTTTTTAAACCAAATCCAGGCTCTGAGACACCCCTGCTATGACTTTACTGGACCTATATACTACGGAGAAGAGGGCTTGTTTGGCAGCACAGCTTTTGTGTGTGTGTGTGTGGGGGGGGGGGGGGGGGGTGGGGGTGGGGGGTGTCAGACCCCGAGGCTACCCTTCTGGGGTTTTGTTTAGATTCTTCTGTGCAAGTTCCCACCAAGTTTTCAGTGAGTTAGACTTTCAGGGTAATCCCGGTTTGCCTTGCAAAGCTGAGGAAATGAGGAATGTCATGGGTATCAATATCTTTCCTGAGAACTGAAATTTACCAGTAGCACATTCCTATCGGGCATTTATCAAAACAAACATATAGGTCAAATGACCGTTTACTCAGCAGCCACCTGTCAGTGCCTTGGGACGTCAGCAGTCCCGGGGATGCCCGTCAAAGCAGCTGAGTGATCTCGGCCCTGAACATTGCCAACTTCCCGCCTGATGCCCTGGAGCTGAGGCTCCAGCTTGAGGCTTCCCAGTGCTTTGGATAAAAAGATGGGTCTGCACCTGAGAGCCATATACATCTCGCTCTGGGGAGAAAACCCAGGCAGGGCTTTCACACAGAGGAAGAGTGAGGCCAGGATTTATACTCAGAGAGGATGGGCAAATAGCATTCACGCGTGTGTCCTCCAGGGTGGAATCTGCCAGGCTGACTTACAGTGCAGGAATTTTGCCAGTCAGTGGAGGTAGATGTGTTTATAGATAACGTTAAAAAAAAAAAAAAAGGCAAAAACCAAGACCTCTGCAACCTTGGTATTTGAAACACAGAACACGGGTTTGTGTACTAAGTAAAATGTGGTTAAATAAGATTAATCTATTGATGGTGGGGAAGAGAGCCTGTAAGGCAGATTGGGTGTCTCAGAGGCATAATTATTTATCAAACTAATTAAAATGTAACCACCGGCTCATCGACTCTGGAAACCGGCTTTCCTGCCCAAGAAGCAGCCCTCTCCTGAGCACCTGTGATTTTATTTTTCATTCATTATCCCACTGAAACACAAAACAAAGACCCAAATTGGCCACCTGGAGTGGCCGCTCTCACCCTCTGCCTGGCCATGGCCTGCCTGGACTGGGTGTTGGGGAGGGATGCATGGGGACAGTGTTGTTCACAGTGTGGGTGAGTGAGTGAGCCCACCTGGTCCCCCATGAAGCCCTGAGATGCCAGGAGCTCTACTGGCTGTGGCCTGGGGCTGTGTGATGGAGAGCACGTGTGTTCATCAGCGCCCAGAGGCCTGGGCTATTCAGTTGGGCCATCTCTGCAGCCCAGAGCTGCCTGGTAGGGCTTTGCCATGGATTCCTGGCTTGGAGGAGTCATGATTAGCATGAACTTGCAAACCCAGACCAGGCGAATTGACTTGGGGG"


def test_evaluate_cluster_phased():
    contig = "1"
    ref = Mock()
    ref.get = lambda a, b, c: REGION1_REF
    diff = hap_eval.VCFEvaluator.diff_by_len

    # Variant alleles are < minsize
    call_variants = [
        vcflib.Variant(
            "1", 7395183, '.', "T", ["TTTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG"],
            100.0, [], {}, [{'GT': '1|0'}], 7395184, None
        ),
        vcflib.Variant(
            "1", 7395198, '.', "TGGGGGGGGGGGGGG", ["T"], 100.0, [], {},
            [{'GT': '1|0'}], 7395213, None
        ),
        vcflib.Variant(
            "1", 7395199, '.', "G", ["GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT", "*"],
            100.0, [], {}, [{'GT': '2|1'}], 7395200, None
        ),
    ]
    g = [(1, x) for x in call_variants]
    call = hap_eval.VCFEvaluator.evaluate_cluster(g, contig, ref, diff)
    assert call is None

    # Phased variants, together, are larger than minsize
    call_variants = [
        vcflib.Variant(
            "1", 7395183, '.', "T", ["TTTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG"],
            100.0, [], {}, [{'GT': '1|0'}], 7395184, None
        ),
        vcflib.Variant(
            "1", 7395198, '.', "TGGGGGGGGGGGGGG", ["T"], 100.0, [], {},
            [{'GT': '0|1'}], 7395213, None
        ),
        vcflib.Variant(
            "1", 7395199, '.', "G", ["GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT", "*"],
            100.0, [], {}, [{'GT': '1|2'}], 7395200, None
        ),
    ]
    g = [(1, x) for x in call_variants]
    call = hap_eval.VCFEvaluator.evaluate_cluster(g, contig, ref, diff)
    assert call is "FP"

    # Variant is an FP - size > minsize
    call_variants = [
        vcflib.Variant(
            "1", 7395183, '.', "T", ["TTTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGATGTCGATACGATGCTAGCTACTGATCG"],
            100.0, [], {}, [{'GT': '1|0'}], 7395184, None
        ),
    ]
    g = [(1, x) for x in call_variants]
    call = hap_eval.VCFEvaluator.evaluate_cluster(g, contig, ref, diff)
    assert call is "FP"

def test_evaluate_cluster_unphased():
    contig = "1"
    ref = Mock()
    ref.get = lambda a, b, c: REGION1_REF
    diff = hap_eval.VCFEvaluator.diff_by_len

    # Haplotype combination exists where variants are < minsize - no call
    call_variants = [
        vcflib.Variant(
            "1", 7395183, '.', "T", ["TTTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG"],
            100.0, [], {}, [{'GT': '1/0'}], 7395184, None
        ),
        vcflib.Variant(
            "1", 7395198, '.', "TGGGGGGGGGGGGGG", ["T"], 100.0, [], {},
            [{'GT': '0/1'}], 7395213, None
        ),
        vcflib.Variant(
            "1", 7395199, '.', "G", ["GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT", "*"],
            100.0, [], {}, [{'GT': '1/2'}], 7395200, None
        ),
    ]
    g = [(1, x) for x in call_variants]
    call = hap_eval.VCFEvaluator.evaluate_cluster(g, contig, ref, diff)
    assert call is None

    # Unphased variant is > minsize - FP
    call_variants = [
        vcflib.Variant(
            "1", 7395183, '.', "T", ["TTTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGATGTCGATACGATGCTAGCTACTGATCG"],
            100.0, [], {}, [{'GT': '1/0'}], 7395184, None
        ),
    ]
    g = [(1, x) for x in call_variants]
    call = hap_eval.VCFEvaluator.evaluate_cluster(g, contig, ref, diff)
    assert call is "FP"


class MockVCF:
    def __init__(self, variants, *args, **kwargs):
        self.variants = variants

    def range(self, *args, **kwargs):
        return self

    def __iter__(self):
        return iter(self.variants)


def test_altfrac_cluster():
    pred = hap_eval.VCFEvaluator.filter

    variants = [
        vcflib.Variant(
            "1", 10, "1", "A", ["ATGTC"], 100.0, [], {}, [{'GT': '1/0'}], 11, None,
        ),
        vcflib.Variant(
            "1", 20, "2", "A", ["ATGTC"], 100.0, [], {}, [{'GT': '1/0'}], 21, None,
        ),
    ]
    vcf = MockVCF(variants)

    # Variants alt_bases < altfrac
    g = list(
        hap_eval.VCFEvaluator.altfrac_cluster(
            [vcf], '1', 1, 1000, pred, altfrac=0.53
        )
    )
    assert g
    assert len(g) == 2
    assert g[0][0][1].id == "1"
    assert g[1][0][1].id == "2"

    # Variants alt_bases > altfrac - clustered
    g = list(
        hap_eval.VCFEvaluator.altfrac_cluster(
            [vcf], '1', 1, 1000, pred, altfrac=0.51
        )
    )
    assert g
    assert len(g) == 1
    assert g[0][0][1].id == "1"
    assert g[0][1][1].id == "2"

    variants = [
        vcflib.Variant(
            "1", 10, "1", "A", ["ATGTCATGCGTGTTGTAC"], 100.0, [], {},
            [{'GT': '1/0'}], 11, None,
        ),
        vcflib.Variant(
            "1", 20, "2", "A", ["T"], 100.0, [], {}, [{'GT': '1/0'}], 21, None,
        ),
        vcflib.Variant(
            "1", 30, "3", "A", ["T"], 100.0, [], {}, [{'GT': '1/0'}], 31, None,
        ),
        vcflib.Variant(
            "1", 70, "4", "A", [''.join(random.choices("ACGT", k=40))], 100.0,
            [], {}, [{'GT': '1/0'}], 71, None,
        ),
    ]
    vcf = MockVCF(variants)

    # Two large SVs with intervening events
    g = list(hap_eval.VCFEvaluator.altfrac_cluster([vcf], '1', 1, 1000, pred))
    assert g
    assert len(g) == 1
    for i in range(4):
        assert g[0][i][1].id == str(i + 1)

    # Decrease maxdist
    g = list(
        hap_eval.VCFEvaluator.altfrac_cluster(
            [vcf], '1', 1, 1000, pred, maxdist=15
        )
    )
    assert g
    assert len(g) == 2
    assert g[0][0][1].id == "1"
    assert g[0][1][1].id == "2"
    assert g[0][2][1].id == "3"
    assert g[1][0][1].id == "4"

    variants = [
        vcflib.Variant(
            "1", 20, "2", "A", ["T"], 100.0, [], {}, [{'GT': '1/0'}], 21, None,
        ),
        vcflib.Variant(
            "1", 30, "3", "A", ["T"], 100.0, [], {}, [{'GT': '1/0'}], 31, None,
        ),
        vcflib.Variant(
            "1", 40, "4", "A", [''.join(random.choices("ACGT", k=50))], 100.0,
            [], {}, [{'GT': '1/0'}], 71, None,
        ),
    ]
    vcf = MockVCF(variants)

    # Earlier variants are clustered by the large SV
    g = list(hap_eval.VCFEvaluator.altfrac_cluster([vcf], '1', 1, 1000, pred))
    assert g
    assert len(g) == 1
    assert g[0][0][1].id == "2"
    assert g[0][1][1].id == "3"
    assert g[0][2][1].id == "4"

    # Test case from #2
    variants = [
        vcflib.Variant(
            "1", 3045231, "1", "CTCCCACCATGCTACTGGAGACCA", ["C"], 100.0, [],
            {}, [{'GT': '1|1'}], 3045255, None,
        ),
        vcflib.Variant(
            "1", 3045364, "2", "A", ["G"], 100.0, [], {}, [{'GT': '1|1'}], 3045365, None,
        ),
        vcflib.Variant(
            "1", 3045614, "3", "T", ["C"], 100.0, [], {}, [{'GT': '1|1'}], 3045615, None,
        ),
        vcflib.Variant(
            "1", 3045744, "4", "A", ["AGT"], 100.0, [], {}, [{'GT': '1|1'}], 3045745, None,
        ),
        vcflib.Variant(
            "1", 3045864, "5", "AGTGTCCTGCCTGGTGTGTGTGTGGGTGT", ["A"], 100.0,
            [], {}, [{'GT': '1|1'}], 3045893, None,
        ),
    ]
    vcf = MockVCF(variants)
    g = list(
        hap_eval.VCFEvaluator.altfrac_cluster(
            [vcf], '1', 3040000, 3050000, pred
        )
    )
    assert g
    assert len(g) == 5  # Variants do not cluster
