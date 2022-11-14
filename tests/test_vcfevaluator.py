#!/usr/bin/env python

from unittest.mock import Mock

import hap_eval
import vcflib

REGION1_REF = "AACCATTTTTTAAACCAAATCCAGGCTCTGAGACACCCCTGCTATGACTTTACTGGACCTATATACTACGGAGAAGAGGGCTTGTTTGGCAGCACAGCTTTTGTGTGTGTGTGTGTGGGGGGGGGGGGGGGGGGTGGGGGTGGGGGGTGTCAGACCCCGAGGCTACCCTTCTGGGGTTTTGTTTAGATTCTTCTGTGCAAGTTCCCACCAAGTTTTCAGTGAGTTAGACTTTCAGGGTAATCCCGGTTTGCCTTGCAAAGCTGAGGAAATGAGGAATGTCATGGGTATCAATATCTTTCCTGAGAACTGAAATTTACCAGTAGCACATTCCTATCGGGCATTTATCAAAACAAACATATAGGTCAAATGACCGTTTACTCAGCAGCCACCTGTCAGTGCCTTGGGACGTCAGCAGTCCCGGGGATGCCCGTCAAAGCAGCTGAGTGATCTCGGCCCTGAACATTGCCAACTTCCCGCCTGATGCCCTGGAGCTGAGGCTCCAGCTTGAGGCTTCCCAGTGCTTTGGATAAAAAGATGGGTCTGCACCTGAGAGCCATATACATCTCGCTCTGGGGAGAAAACCCAGGCAGGGCTTTCACACAGAGGAAGAGTGAGGCCAGGATTTATACTCAGAGAGGATGGGCAAATAGCATTCACGCGTGTGTCCTCCAGGGTGGAATCTGCCAGGCTGACTTACAGTGCAGGAATTTTGCCAGTCAGTGGAGGTAGATGTGTTTATAGATAACGTTAAAAAAAAAAAAAAAGGCAAAAACCAAGACCTCTGCAACCTTGGTATTTGAAACACAGAACACGGGTTTGTGTACTAAGTAAAATGTGGTTAAATAAGATTAATCTATTGATGGTGGGGAAGAGAGCCTGTAAGGCAGATTGGGTGTCTCAGAGGCATAATTATTTATCAAACTAATTAAAATGTAACCACCGGCTCATCGACTCTGGAAACCGGCTTTCCTGCCCAAGAAGCAGCCCTCTCCTGAGCACCTGTGATTTTATTTTTCATTCATTATCCCACTGAAACACAAAACAAAGACCCAAATTGGCCACCTGGAGTGGCCGCTCTCACCCTCTGCCTGGCCATGGCCTGCCTGGACTGGGTGTTGGGGAGGGATGCATGGGGACAGTGTTGTTCACAGTGTGGGTGAGTGAGTGAGCCCACCTGGTCCCCCATGAAGCCCTGAGATGCCAGGAGCTCTACTGGCTGTGGCCTGGGGCTGTGTGATGGAGAGCACGTGTGTTCATCAGCGCCCAGAGGCCTGGGCTATTCAGTTGGGCCATCTCTGCAGCCCAGAGCTGCCTGGTAGGGCTTTGCCATGGATTCCTGGCTTGGAGGAGTCATGATTAGCATGAACTTGCAAACCCAGACCAGGCGAATTGACTTGGGGG"


def test_evaluate_cluster_phased1():
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
    contig = "1"
    ref = Mock()
    ref.get = lambda a, b, c: REGION1_REF
    diff = hap_eval.VCFEvaluator.diff_by_len
    call = hap_eval.VCFEvaluator.evaluate_cluster(g, contig, ref, diff)
    assert call is None

def test_evaluate_cluster_phased2():
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
    contig = "1"
    ref = Mock()
    ref.get = lambda a, b, c: REGION1_REF
    diff = hap_eval.VCFEvaluator.diff_by_len
    call = hap_eval.VCFEvaluator.evaluate_cluster(g, contig, ref, diff)
    assert call is "FP"

def test_evaluate_cluster_unphased1():
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
    contig = "1"
    ref = Mock()
    ref.get = lambda a, b, c: REGION1_REF
    diff = hap_eval.VCFEvaluator.diff_by_len
    call = hap_eval.VCFEvaluator.evaluate_cluster(g, contig, ref, diff)
    assert call is None

def test_evaluate_cluster_phased3():
    call_variants = [
        vcflib.Variant(
            "1", 7395183, '.', "T", ["TTTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGATGTCGATACGATGCTAGCTACTGATCG"],
            100.0, [], {}, [{'GT': '1|0'}], 7395184, None
        ),
    ]
    g = [(1, x) for x in call_variants]
    contig = "1"
    ref = Mock()
    ref.get = lambda a, b, c: REGION1_REF
    diff = hap_eval.VCFEvaluator.diff_by_len
    call = hap_eval.VCFEvaluator.evaluate_cluster(g, contig, ref, diff)
    assert call is "FP"

def test_evaluate_cluster_unphased2():
    call_variants = [
        vcflib.Variant(
            "1", 7395183, '.', "T", ["TTTTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGATGTCGATACGATGCTAGCTACTGATCG"],
            100.0, [], {}, [{'GT': '1/0'}], 7395184, None
        ),
    ]
    g = [(1, x) for x in call_variants]
    contig = "1"
    ref = Mock()
    ref.get = lambda a, b, c: REGION1_REF
    diff = hap_eval.VCFEvaluator.diff_by_len
    call = hap_eval.VCFEvaluator.evaluate_cluster(g, contig, ref, diff)
    assert call is "FP"
