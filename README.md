# Hap-Eval
An open-source VCF comparison engine for structural variant benchmarking

*Note: hap-eval is pre-release software and is under active development. If hap-eval does
not work with your VCF, please file an issue* 

## Installation
You can use the following commands to clone the repository and vcflib submodule
to your machine and install the tool.
```
git clone --recurse-submodules https://github.com/Sentieon/hap-eval.git
pip install ./hap-eval
```

## Usage
```
usage: hap_eval [-h] -r FASTA -b VCF -c VCF [-i BED] [-t INT] [--base_out VCF]
                [--comp_out VCF] [--maxdist INT] [--minsize INT]
                [--maxdiff FLOAT] [--metric STR]

optional arguments:
  -h, --help            show this help message and exit
  -r FASTA, --reference FASTA
                        Reference file
  -b VCF, --base VCF    Baseline vcf file
  -c VCF, --comp VCF    Comparison vcf file
  -i BED, --interval BED
                        Evaluation region file
  -t INT, --thread_count INT
                        Number of threads
  --base_out VCF        Annotated baseline vcf file
  --comp_out VCF        Annotated comparison vcf file
  --maxdist INT         Maximum distance to cluster variants (default: 1000)
  --minsize INT         Minimum size of variants to consider (default: 50)
  --maxdiff FLOAT       Haplotype difference threshold (default: 0.2)
  --metric STR          Distance metric (default: Levenshtein)
```

### Clustering, `--maxdist` and `--minsize`

hap-eval evaluates groups (clusters) of variants that are nearby on the reference genome. Variants are grouped together if the distance between two variants is less than `--maxdist`. The entire group will be evaluated if the largest possible haplotype that can be constructed from the group has an insertion or deletion variant larger than `--minsize`.

Some examples are shown below. For the purpose of demonstration, these examples use a baseline VCF (`--base`), but an empty comparison VCF `--comp`.

#### Basic grouping

The baseline VCF has the following variants:
```txt
chr1    1541864 .       T       C       30      .       LCR=0   GT:AD   0|1:1,1
chr1    1542773 .       T       C       30      .       LCR=0   GT:AD   0|1:1,1
chr1    1542793 .       C       G       30      .       LCR=0   GT:AD   0|1:1,1
chr1    1542800 .       T       C       30      .       LCR=0   GT:AD   0|1:1,1
chr1    1543500 .       T       G       30      .       LCR=0   GT:AD   0|1:1,1
chr1    1543725 .       C       CCTCTGTCACGGCCTCGGCCACTCCCCACTGTCACGGCCTCGGCCACTCCCCTCTGTCACGGCCTCGGCCACTCCCCTCTGTCACGGCCTCGGCCACTCCCCACTGTCACGGCCTCGGCCACTCCCCACTGTCACGGCCTCGGCCACTCCCCTCTGTCACGGCCTCGGCCACTCCCCACTGTCACGGCCTCGGCCACTCCCCT       30      .       TRF;TRFdiff=8.1;TRFrepeat=TCACGGCCTCGGCCACTCCCCTCTG;TRFovl=1;TRFstart=1543581;TRFend=1543840;TRFperiod=25;TRFcopies=18.7;TRFscore=672;TRFentropy=1.75;TRFsim=0.97;SVTYPE=INS;SVLEN=202;LCR=0.864133;REMAP=tandem        GT:AD   0|1:1,1
chr1    1543953 .       A       G       30      .       LCR=0   GT:AD   0|1:1,1
```

These variants are grouped together as they are all less than `--maxdist` (1000bp) from each other, evaluated as a single region by hap-eval, and contributed one `FN` to the overall evaluation (as there is no comparison VCF).

```txt
FN chr1:1541764-1544053 7 0 2290 (0,202) (0,0) 1.0
```

#### Variants further than `--maxdist`

The baseline VCF has the following variants:
```txt
chr1    1981820 .       A       ACCGCGGACAGACACGGGGGCACGCAGGACACCCAGCCGCGGACAGACACGGGGGCACGCAGGACACCCAG         30      .       TRF;TRFdiff=2;TRFrepeat=ACACCCAGCCGCGGACAGACACGGGGGCACGCAGG;TRFovl=1;TRFstart=1981391;TRFend=1982440;TRFperiod=35;TRFcopies=31.9;TRFscore=2263;TRFentropy=1.71;TRFsim=1;SVTYPE=INS;SVLEN=70;LCR=0.785201;REMAP=tandem   GT:AD   1|1:0,2
chr1    1986432 .       T       A       30      .       LCR=0   GT:AD   1|1:0,2
chr1    1986564 .       G       A       30      .       LCR=0   GT:AD   0|1:1,1
chr1    1991395 .       T       C       30      .       LCR=0   GT:AD   1|1:0,2
chr1    1991761 .       T       C       30      .       LCR=0   GT:AD   1|0:1,1
chr1    1991778 .       A       G       30      .       LCR=0   GT:AD   0|1:1,1
chr1    1991822 .       G       A       30      .       LCR=0   GT:AD   1|1:0,2
chr1    1991835 .       C       T       30      .       LCR=0   GT:AD   1|1:0,2
chr1    1991841 .       T       C       30      .       LCR=0   GT:AD   1|1:0,2
chr1    1991867 .       A       G       30      .       LCR=0   GT:AD   0|1:1,1
chr1    1992226 .       A       G       30      .       LCR=0   GT:AD   1|1:0,2
chr1    1992242 .       A       G       30      .       LCR=0   GT:AD   0|1:1,1
chr1    1993704 .       A       AGGGCACAGTGGCTCATGCCTGTAATCCCAGCAACATGGGAGCCTGAGGTGGGAGGCTCTCTTGAGGCCAGGAGTTTGAGACCAGCCTGGGCAACATAGTGAGACCCCCCACCCCCCGCCATTTCTAGGAAAAAAAAAAAAAGTGGCC    30      .       SVTYPE=INS;SVLEN=147;RM_score=113;RM_repeat=FLAM_C;RM_clsfam=SINE/Alu;LCR=0.98366;REMAP=partial GT:AD   1|1:0,2
chr1    1993887 .       C       A       30      .       LCR=0   GT:AD   1|1:0,2
chr1    1993898 .       C       T       30      .       LCR=0   GT:AD   1|1:0,2
chr1    1994025 .       C       T       30      .       LCR=0   GT:AD   1|1:0,2
```

These variants are split into two groups, as the distance between the variant at `chr1:1981820` and the variant at `chr1:1986432` greater than `--maxdist` (default 1000). Many of the SNVs are not included in either group.

```txt
FN chr1:1981720-1981920 1 0 201 (70,70) (0,0) 1.0
FN chr1:1993604-1994125 4 0 522 (147,147) (0,0) 1.0
```

#### Grouping of deletions

The baseline VCF has the following variants:
```txt
chr1    3643863 .       CGGCCCCGCTGAGCCCCGGACATGGGGTCGCGCCTTCGGAAGGGGACCCAGCGGCCCCGCTGAGCCCCGGACATGGGGTCGCGCCTTCGGAAGGGGACCCAGCGGCCCCGCTGAGCCCCGGACATGGGGTCGCGCCTTCGGAAGGGGACCCAGCGGCCCCGCTGAGCCCCGGACATGGGGTCGCGCCTTCGGAAGGGGACCCAGCGGCCCCGCTGAGCCCCGGACATGGGGTCGCGCCTTCGGAAGGGGACCCAGCGGCCCCGCTGAGCCCCGGACATGGGGTCGCGCCTTCGGAAGGGGACCCAGCGGCCCCGCTGAGCCCCGGACATGGGGTCGCGCCTTCGGAAGGGGACCCAGCGGCCCCGCTGAGCCCCGGACATGGGGTCGCGCCTTCGGAAGGGGACCCAGCGGCCCCGCTGAGCCCCGGACATGGGGTCGCGCCTTCGGAAGGGGACCCAGCGGCCCCGCTGAGCCCCGGACATGGGGTCGCGCCTTCGGAAGGGGATCCAGT C       30      .       TRF;TRFdiff=-10;TRFrepeat=TGGGGTCGCGCCTTCGGAAGGGGACCCAGCGGCCCCGCTGAGCCCCGGACA;TRFovl=1;TRFstart=3643528;TRFend=3644432;TRFperiod=51;TRFcopies=7.8;TRFscore=2578;TRFentropy=1.8;SVTYPE=DEL;SVLEN=510;LCR=0.894375;REMAP=tandem   GT:AD   1|1:0,2
chr1    3645341 .       AGCGGCGCAGCGGGATGGGCGGGTTGCCCCGTGGTGTGCGTGGCGCAGCGGGACGGGCGGGTTGCCCCGTGGTGTGC   A       30      .       TRF;TRFdiff=-2;TRFrepeat=CGTGGTGTGCGCGGCGCAGCGGGACGGGCGGGTTGCCC;TRFovl=1;TRFstart=3645293;TRFend=3645554;TRFperiod=38;TRFcopies=4.9;TRFscore=666;TRFentropy=1.66;SVTYPE=DEL;SVLEN=76;LCR=0.83505;REMAP=tandem   GT:AD   1|1:0,2
chr1    3645458 .       C       G       30      .       TRF;LCR=0       GT:AD   1|1:0,2
chr1    3645460 .       T       C       30      .       TRF;LCR=0       GT:AD   1|1:0,2
chr1    3646061 .       C       G       30      .       LCR=0   GT:AD   0|1:1,1
```

The first deletion extends 510bp over the reference genome. This extension allows the two deletions to be grouped, even if `chr1:3645341` is more than `--maxdist` from `chr1:3643863`.

```txt
FN chr1:3643763-3646161 5 0 2399 (-586,-586) (0,0) 1.0
```

#### Multiple variants exceeding `--minsize`

The baseline VCF has the following variants:
```txt
chr1    5700440 .       A       G       30      .       LCR=0   GT:AD   1|1:0,2
chr1    5700987 .       A       T       30      .       LCR=0   GT:AD   1|1:0,2
chr1    5701753 .       G       GAGAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAA     30      .       TRF;TRFdiff=0;TRFrepeat=GAAA;TRFovl=1;TRFstart=5701732;TRFend=5701949;TRFperiod=4;TRFcopies=54.2;TRFscore=297;TRFentropy=1.01;SVTYPE=INS;SVLEN=34;RM_score=34;RM_repeat=(AGAA)N;RM_clsfam=Simple_repeat;LCR=0.43156;REMAP=interspersed  GT:AD   1|1:0,2
chr1    5701872 .       G       GAGGGAGGAAGGAAGGAAGGA   30      .       TRF;TRFdiff=5;TRFrepeat=GAAG;TRFovl=1;TRFstart=5701732;TRFend=5701921;TRFperiod=4;TRFcopies=50;TRFscore=226;TRFentropy=1.06;TRFsim=0.925;SVTYPE=INS;SVLEN=20;LCR=0.492614;REMAP=interspersed    GT:AD   0|1:1,1
chr1    5701908 .       A       AAGGAAGGAAGGAAGGAAGGAAGGAAGGG   30      .       TRF;TRFdiff=7;TRFrepeat=GAAG;TRFovl=1;TRFstart=5701732;TRFend=5701921;TRFperiod=4;TRFcopies=52;TRFscore=226;TRFentropy=1.06;TRFsim=0.946;SVTYPE=INS;SVLEN=28;RM_score=30;RM_repeat=(AAGG)N;RM_clsfam=Simple_repeat;LCR=0.499571;REMAP=interspersed      GT:AD   1|0:0,1
chr1    5701908 .       A       G       30      .       TRF;LCR=0       GT:AD   0|1:0,1
chr1    5702582 .       G       A       30      .       LCR=0   GT:AD   1|0:1,1
```

Individually, each of these variants is smaller than `--minsize`. However, they can be combined together into a haplotype larger than `--minsize`. Accordingly, they are evaluated as a group by hap-eval.

```txt
FN chr1:5700340-5702682 7 0 2343 (62,54) (0,0) 1.0
```

#### Variants less than `--minsize`

The baseline VCF has the following variants:
```txt
chr1    10575271        .       C       T       30      .       LCR=0   GT:AD   1|1:0,2
chr1    10576386        .       G       A       30      .       LCR=0   GT:AD   1|1:0,2
chr1    10576402        .       G       A       30      .       LCR=0   GT:AD   1|1:0,2
chr1    10576759        .       C       CT      30      .       LCR=0.5 GT:AD   1|1:0,2
chr1    10576898        .       A       T       30      .       LCR=0   GT:AD   1|1:0,2
chr1    10577189        .       C       CAGG    30      .       LCR=0.75        GT:AD   1|1:0,2
chr1    10577410        .       C       CAAAA   30      .       LCR=0.360964    GT:AD   1|1:0,2
chr1    10577546        .       ATATATATATATATATATTTTTTT        A       30      .       TRF;TRFdiff=-10.5;TRFrepeat=TA;TRFovl=0.913;TRFstart=10577524;TRFend=10577562;TRFperiod=2;TRFcopies=9;TRFscore=107;TRFentropy=1.15;SVTYPE=DEL;SVLEN=23;LCR=0.477217;REMAP=interspersed  GT:AD   1|0:1,1
chr1    10577560        .       ATATTTTT        *       30      .       TRF;LCR=0.405639        GT:AD   1|0:0,1
chr1    10577560        .       ATATTTTT        A       30      .       TRF;LCR=0.405639        GT:AD   0|1:0,1
```

Individually, each of these variants is smaller than `--minsize`. In this case, these variants cannot be combined into an insertion or deletion larger than `--minsize` as some variants are deletions while others are insertions. The total SV size is less than `--minsize`, and this site does not pass hap-eval's `--minsize` filter.

### Understanding the `--maxdiff` argument

hap-eval will use the `--maxdiff` (0.2 by default) argument to determine a distance threshold between the assembled haplotypes in the baseline and comparison VCFs. Initially, hap-eval will generate all possible haplotypes in the baseline and comparison VCFs. The "best" haplotype combination (containing a haplotype pair with the smallest distance between the baseline and comparison VCFs) is selected for evaluation.

Evaluation compares the two haplotypes from the baseline sample, with the two haplotypes from the comparison sample. If the two pairs both have a distance smaller than `--maxdiff`, the site will be evaluated as `MM`. If one pair has a distance greater than `--maxdiff`, the site is evaluated as `MX`. If both pairs have a distance greater than `--maxdist`, the site is evaluated as `XX`.

By default, hap-eval will use the Levenshtein distance if the `Levenshtein` python module is available. If the `Levenshtein` module is not available, hap-eval will fallback to using a simple length metric. **Using the Levenshtein distance is highly recommended as the length metric will not evaluate base differences between the two haplotypes.**

#### Example `MX` evaluation

The baseline VCF has the following variants:
```txt
chr2    130061662       .       CAGCCATCAGGAAGCACCAATCAGGCTGTTGGAACCTGGCCATCCGTGCTTCCTTCAGTGGCTGAAGTTGGTGGCTGGTCCACCTGCTCCTGGCACACCTTTGCAGAGGTGGCTGGTTGCTCTTTGAGCCAGCTTGGCCCTGCCTGGCATGCATAGGCCTCAGCTACTGACACACTGCTCCAAGTGAGCTTGTCCTGCATTGGCACAAATTCTGAGTCTGGCCAGGGTCACAGAAGGCCAAGTCCCCTGGATGGTTATCCGGGCTGCTTTCTGCACTTGAACATAAAGTCTTCCTCAAGACAGCCTGTGGTCTGCCTCTTGGCAACCAAGAAGCCCACAGTGCCATACGACCCGAGGCATGGACTGGAGCCCCAAAGGCAGCACACACCCGGCTCCTGAGCCTACTGCTCGTTTCCTCTCTGTGGCTCCATTTGTAGCACAGTTGTTGCACTGAGGCTTGTGCATGCCAGGGAAGGCCAAGCTGGCTCAAAGAGCAACCAGCCACCTCTGCAAGGGTGTGCCAGGAGCAGGTGCACCAGTCACCAACTAGTGGCCGGACATGGTACATCAGTTCTTCTACCCTGAAGGTGGGCCACAGGGCCATCTGCTTTTCCTAAGGCCTCTGCTCCATG        C       30      .       TRF;TRFdiff=-4.2;TRFrepeat=GCAT;TRFovl=0.027;TRFstart=130061807;TRFend=130061814;TRFperiod=4;TRFcopies=-2.2;TRFscore=24;TRFentropy=2;SVTYPE=DEL;SVLEN=631;LCR=0.99135;REMAP=tandem      GT:AD   1|1:0,2
chr2    130062520       .       T       C       30      .       TRF;LCR=0       GT:AD   1|1:0,2
chr2    130062768       .       AG      A       30      .       TRF;LCR=0.5     GT:AD   1|1:0,2
```

The comparison VCF has the following variants:
```txt
chr2    130061662       .       CAGCCATCAGGAAGCACCAATCAGGCTGTTGGAACCTGGCCATCCGTGCTTCCTTCAGTGGCTGAAGTTGGTGGCTGGTCCACCTGCTCCTGGCACACCTTTGCAGAGGTGGCTGGTTGCTCTTTGAGCCAGCTTGGCCCTGCCTGGCATGCATAGGCCTCAGCTACTGACACACTGCTCCAAGTGAGCTTGTCCTGCATTGGCACAAATTCTGAGTCTGGCCAGGGTCACAGAAGGCCAAGTCCCCTGGATGGTTATCCGGGCTGCTTTCTGCACTTGAACATAAAGTCTTCCTCAAGACAGCCTGTGGTCTGCCTCTTGGCAACCAAGAAGCCCACAGTGCCATACGACCCGAGGCATGGACTGGAGCCCCAAAGGCAGCACACACCCGGCTCCTGAGCCTACTGCTCGTTTCCTCTCTGTGGCTCCATTTGTAGCACAGTTGTTGCACTGAGGCTTGTGCATGCCAGGGAAGGCCAAGCTGGCTCAAAGAGCAACCAGCCACCTCTGCAAGGGTGTGCCAGGAGCAGGTGCACCAGTCACCAACTAGTGGCCGGACATGGTACATCAGTTCTTCTACCCTGAAGGTGGGCCACAGGGCCATCTGCTTTTCCTAAGGCCTCTGCTCCATG        C       60      .       SVTYPE=DEL;SVLEN=-631;AF=0.2;END=130062293
      GT:AD:PID:PGT   0/1:20,5:chr2_130060323_130063933:0|1
chr2    130062308       .       A       AGCCAAGGCAGGACAAGCTCACTCAGAGCAGCGTGTTAGTATCTGGGGCCTGTGCATGCCAGGAAGGCCAAGCTGGCTCAAAGAGCAACCAGCCACCTCTGCAAGGTCTGGCCAGGGCTACAGAAGGCCCCCCTGGATGGTAATCCTGGCTGCTTTCTGTACTTGAACATAAAGTCTTCCTCAAGACGGCCTGTGGTCTGCCTCCAGGCAAGTAAGAAGCCCGCGGTGCTATACAACTCGAGGCATGGACTGGAGCCCCAAAGGCAGCGCACACACTGCTCCTGAGCCTGCTGCTCGTTTCCTCTCTATGGCTCCATTTGTAGCACTGTTCTTGCACTGAGGCTTGTGCATGCCGGGCAAGGCCAAGCTGACTCAAAGAGCAACCAGTCACTTCTGCAAGGGTGCGCCAGGAGCCGGTGCACCAGCCACCAACCTCACTTGCTACCGGACATGGCACATCAGTACTTCTACCCCAAAGGTAGGGCCACAGGGCCATCTGCTTTTCCTAAGGCCTCTGCTCCATCAGCCATCAGGAGGCAG    60      .       SVTYPE=INS;SVLEN=539;AF=0.8;END=130062308       GT:AD:PID:PGT   0/1:5,20:chr2_130060323_130063933:1|0
chr2    130062642       .       A       ACCCGAGGCATGGACTGGAGCCCCAAAGGCAGCGCACATCCTGCTCCTAACTCTGCCGCTCATTTCCTCTCTGTGGCTCCATTTGTAGCACAGTTGTTGCACTGAGACTTGTGCATGCCAGGCAAGGCCAAGCTGGCTCAAACAGCAACCAGCCACCTCTGCAAGTGTGTGCCAGGAGCAGCCAGACCAGCCACCAACCTCACTCGCTGCGGGACATGGTACATTGGTTCTTCTACCCTAAAGGCAGGGCCAAGAGGCAGACCACAGGCCGTCTTGAGGAGGACTTTATGTTCAAGTGCAGAAAGCAGCCAGGATTGCCACCAAGGGGACTCAGCCTTCTGTGGTCTGCACTGCCATACAAGCTCTGAGACGTGGACTAGTGACATCTGCTTTATAGAAAAATTAACCTAAGATCTATTAAAGAGTTAAACATGCCATCTGATTTTACTCAGGCCTCTGCTCCATCAGCCCTCAGGTGGCAGCCACTCAGGCTGTGGGAACCTGGCCATCCCTGCTTCCTTCAGTGGGTGAGGTTGGTGGCTGGTCCAACTGGTCCAGGCGCACCCTTGTAGAGGTGGCTGGTTGCTCTTTGAGCCAGCTTGGCCTTCCCGGACATGCACAGGCCCCAGGTACTAACACGCTGCTCTGAGTGAGCTTCTCCTGCCTTGACACAAATTCTAAGTCTCGCCACGGCCACAAAAGGCCGAGTCCCCTGCGTGGCAATCATGGCTGCTTTCTGCACTTGAACATAAAGTCTTCCTCAAGACAGCCTGTGGTCTGCCTCTTGGCAACCAAGAAGCCCACAGTGCCATAC  60      .       SVTYPE=INS;SVLEN=813;AF=0.8;END=130062642       GT:AD:PID:PGT   0/1:5,20:chr2_130060323_130063933:1|0
```

hap-eval assembles the variants into haplotypes. Using the distance metric, it decides that these are the best haplotype combinations:
```txt
# base-hap1 and base-hap2 contain the followig variants
chr2:130061662
chr2:130062520
chr2:130062768
# Note that base-hap1 and base-hap2 are identical as all variants are homozygous for the alternate allele

# comp-hap1
chr2:130061662

# comp-hap2
chr2:130062308
chr2:130062642
```

Here are the full assembled haplotype sequences for these haplotypes across this region. Note that 100bp of padding before the first variant and after the last variant is included in the sequences.

```txt
# base-hap1 and base-hap2
AGCCACCAACCTCACTTGCTGCTGCACATGGCACATCAGTACTTCTACCCTAAAGGTAGGGCCACAGGGCCATCTGCTTTTCCTAAGGCCTCTGCTCCATCAGCCATCAGGAAGCACCAATCAGGCTGTTGGAACCTGGCCATCCGTGCTTCCTTCAGTGGCTGAAGTTGGTGGCTGGTCCACCTGCTCCTGGCACACCTTTGCAGAGGTGGCTGGTTGCTCTTTGAGCCAGCTTGGCCCTGCCTGGCATGCATAGGCCTCAGCTACTGACACACTGCTCCAAGTGAGCTTGTCCTGCATTGGCACAAATTCTGAGTCTGGCCAGGGCCACAGAAGGCCAAGTCCCCTGGATGGTTATCCGGGCTGCTTTCTGCACTTGAACATAAAGTCCTCCTCAAGATGGCCTGTTGTCTGCCTCTTGGCAACCAAGAAGCCTGCAGTGCCATACGAGTTCTGAGGCATGGACTAGAGGCCCAAAGGCAGAGCACACCCTGCTCCTGACCCTGCGGCTCGTTTCCTCTCTGTGGCTCCATTTGTAGCACAGTTGTTGCACTGAGGCTTGTGCATGCCAGGGAAGGCCAAGCTGGCTCAAAGAGCAACCAGCCACCTCTGCAAGGGTGTGCCAGGAGCAGGTGCACCACTCACCCACCAGCGGCCGGACATGGTACATCAGTTC

# comp-hap1
AGCCACCAACCTCACTTGCTGCTGCACATGGCACATCAGTACTTCTACCCTAAAGGTAGGGCCACAGGGCCATCTGCTTTTCCTAAGGCCTCTGCTCCATCAGCCATCAGGAAGCACCAATCAGGCTGTTGGAACCTGGCCATCCGTGCTTCCTTCAGTGGCTGAAGTTGGTGGCTGGTCCACCTGCTCCTGGCACACCTTTGCAGAGGTGGCTGGTTGCTCTTTGAGCCAGCTTGGCCCTGCCTGGCATGCATAGGCCTCAGCTACTGACACACTGCTCCAAGTGAGCTTGTCCTGCATTGGCACAAATTCTGAGTCTGGCCAGGGTCACAGAAGGCCAAGTCCCCTGGATGGTTATCCGGGCTGCTTTCTGCACTTGAACATAAAGTCTTCCTCAAGACAGCCTGTGGTCTGCCTCTTGGCAACCAAGAAGCCCACAGTGCCATACGACCCGAGGCATGGACTGGAGCCCCAAAGGCAGCACACACCCGGCTCCTGAGCCTACTGCTCGTTTCCTCTCTGTGGCTCCATTTGTAGCACAGTTGTTGCACTGAGGCTTGTGCATGCCAGGGAAGGCCAAGCTGGCTCAAAGAGCAACCAGCCACCTCTGCAAGGGTGTGCCAGGAGCAGGTGCACCAGTCACCAACTAGTGGCCGGACATGGTACATCAGTTCTTCTACCCTGAAGGTGGGCCACAGGGCCATCTGCTTTTCCTAAGGCCTCTGCTCCATGAGCCATCAGGAAGCAGCCAAGGCAGGACAAGCTCACTCAGAGCAGCGTGTTAGTATCTGGGGCCTGTGCATGCCAGGAAGGCCAAGCTGGCTCAAAGAGCAACCAGCCACCTCTGCAAGGTCTGGCCAGGGCTACAGAAGGCCCCCCTGGATGGTAATCCTGGCTGCTTTCTGTACTTGAACATAAAGTCTTCCTCAAGACGGCCTGTGGTCTGCCTCCAGGCAAGTAAGAAGCCCGCGGTGCTATACAACTCGAGGCATGGACTGGAGCCCCAAAGGCAGCGCACACACTGCTCCTGAGCCTGCTGCTCGTTTCCTCTCTATGGCTCCATTTGTAGCACTGTTCTTGCACTGAGGCTTGTGCATGCCGGGCAAGGCCAAGCTGACTCAAAGAGCAACCAGTCACTTCTGCAAGGGTGCGCCAGGAGCCGGTGCACCAGCCACCAACCTCACTTGCTACCGGACATGGCACATCAGTACTTCTACCCCAAAGGTAGGGCCACAGGGCCATCTGCTTTTCCTAAGGCCTCTGCTCCATCAGCCATCAGGAGGCAGCCAATCAGGCTGTTGGAACCTGGCCATCCGTGCTTCCTTCAGTGGCTGAAGTTGGTGGCTGGTCCACCTGCTCCTGGCACACCTTTGCAGAGGTGGCTGGTTGCTCTTTGAGCCAGCTTGGCCCTGCCTGGCATGCATAGGCCTCAGCTACTGACACACTGCTCCAAGTGAGCTTGTCCTGCATTGGCACAAATTCTGAGTCTGGCCAGGGTCACAGAAGGCCAAGTCCCCTGGATGGTTATCCGGGCTGCTTTCTGCACTTGAACATAAAGTCCTCCTCAAGATGGCCTGTTGTCTGCCTCTTGGCAACCAAGAAGCCTGCAGTGCCATACGACCCGAGGCATGGACTGGAGCCCCAAAGGCAGCGCACATCCTGCTCCTAACTCTGCCGCTCATTTCCTCTCTGTGGCTCCATTTGTAGCACAGTTGTTGCACTGAGACTTGTGCATGCCAGGCAAGGCCAAGCTGGCTCAAACAGCAACCAGCCACCTCTGCAAGTGTGTGCCAGGAGCAGCCAGACCAGCCACCAACCTCACTCGCTGCGGGACATGGTACATTGGTTCTTCTACCCTAAAGGCAGGGCCAAGAGGCAGACCACAGGCCGTCTTGAGGAGGACTTTATGTTCAAGTGCAGAAAGCAGCCAGGATTGCCACCAAGGGGACTCAGCCTTCTGTGGTCTGCACTGCCATACAAGCTCTGAGACGTGGACTAGTGACATCTGCTTTATAGAAAAATTAACCTAAGATCTATTAAAGAGTTAAACATGCCATCTGATTTTACTCAGGCCTCTGCTCCATCAGCCCTCAGGTGGCAGCCACTCAGGCTGTGGGAACCTGGCCATCCCTGCTTCCTTCAGTGGGTGAGGTTGGTGGCTGGTCCAACTGGTCCAGGCGCACCCTTGTAGAGGTGGCTGGTTGCTCTTTGAGCCAGCTTGGCCTTCCCGGACATGCACAGGCCCCAGGTACTAACACGCTGCTCTGAGTGAGCTTCTCCTGCCTTGACACAAATTCTAAGTCTCGCCACGGCCACAAAAGGCCGAGTCCCCTGCGTGGCAATCATGGCTGCTTTCTGCACTTGAACATAAAGTCTTCCTCAAGACAGCCTGTGGTCTGCCTCTTGGCAACCAAGAAGCCCACAGTGCCATACGTTCTGAGGCATGGACTAGAGGCCCAAAGGCAGAGCACACCCTGCTCCTGACCCTGCGGCTCGTTTCCTCTCTGTGGCTCCATTTGTAGCACAGTTGTTGCACTGAGGCTTGTGCATGCCAGGGAAGGGCCAAGCTGGCTCAAAGAGCAACCAGCCACCTCTGCAAGGGTGTGCCAGGAGCAGGTGCACCACTCACCCACCAGCGGCCGGACATGGTACATCAGTTC

# comp-hap2
AGCCACCAACCTCACTTGCTGCTGCACATGGCACATCAGTACTTCTACCCTAAAGGTAGGGCCACAGGGCCATCTGCTTTTCCTAAGGCCTCTGCTCCATCAGCCATCAGGAAGCACCAATCAGGCTGTTGGAACCTGGCCATCCGTGCTTCCTTCAGTGGCTGAAGTTGGTGGCTGGTCCACCTGCTCCTGGCACACCTTTGCAGAGGTGGCTGGTTGCTCTTTGAGCCAGCTTGGCCCTGCCTGGCATGCATAGGCCTCAGCTACTGACACACTGCTCCAAGTGAGCTTGTCCTGCATTGGCACAAATTCTGAGTCTGGCCAGGGTCACAGAAGGCCAAGTCCCCTGGATGGTTATCCGGGCTGCTTTCTGCACTTGAACATAAAGTCCTCCTCAAGATGGCCTGTTGTCTGCCTCTTGGCAACCAAGAAGCCTGCAGTGCCATACGAGTTCTGAGGCATGGACTAGAGGCCCAAAGGCAGAGCACACCCTGCTCCTGACCCTGCGGCTCGTTTCCTCTCTGTGGCTCCATTTGTAGCACAGTTGTTGCACTGAGGCTTGTGCATGCCAGGGAAGGGCCAAGCTGGCTCAAAGAGCAACCAGCCACCTCTGCAAGGGTGTGCCAGGAGCAGGTGCACCACTCACCCACCAGCGGCCGGACATGGTACATCAGTTC
```

Using the `Levenshtein` module, hap-eval evaluates the distance between the haplotypes. The pair, "base-hap1" and "comp-hap1" have a distance of 0.5947. The pair, "base-hap2" and "comp-hap2" have a distance of 0.0022. 0.5947 is greater than `--maxdiff`, while 0.0022 is less than `--maxdiff`, so the site is evaluated as `MX`:
```txt
MX chr2:130061562-130062869 3 3 1308 (-632,-632) (1352,-631) 0.595
```

## Interpreting the output

### Region evaluation

For each variant region, hap-eval will write the region's evaluation. Here is the example `MX` call:
```txt
MX chr2:130061562-130062869 3 3 1308 (-632,-632) (1352,-631) 0.595
```

Here is the interpretation of each field:
- `MX` - The region is evaluated as an `MX` call
- `chr2:130061562-130062869` - The region covers `chr2:130061562-130062869`
- `3` - The base VCF has three variants
- `3` - The comp VCF has three variants
- `1308` - The region size
- `(-632,-632)` - The best assembled haplotypes in the base VCF have (in total) a 632bp deletion on one haplotype and a 632bp deletion on the other haplotype.
- `(1352,-631)` - The base assembled haplotypes in the comp VCF have (in total) a 1352bp insertion on one haplotype and a 631bp deletion on the other haplotype.
- `0.595` - The maximum of the distance metrics between the haplotype pairs.

### Evaluation Summary
Upon successful evaluation, hap-eval will print an output summary similar to the following:
```
Counter({'MM': 9413, 'FP': 280, 'FN': 183, 'MX': 58, 'XX': 18, '??': 9})
precision 0.9627 recall 0.9723 f1 0.9675
```

The `Counter` contains the number of regions evaluated along with evaluation
result.
* `MM` - A true-positive. The two haplotypes in the query match the two
  haplotypes in the truth VCF. There is at least one variant in both VCFs.
* `FN` - A simple false-negative. There is an event in the truth but no
  corresponding event in the query in the nearby region.
* `FP` - A simple false-positive. There is an event in the query but no
  corresponding event in the truth in the nearby region.
* `??` - The region cannot be evaluated. This is likely due to the large number
  of potential haplotypes at the region.
* `MX` and `XX` - both the truth and the query have some type of event, but the
  difference between at least one of the best assembled haplotype pairs is
  greater than `--maxdiff`.
  * `XX` - both of the best assembled haplotypes in the query have a difference
    greater than `--maxdiff` from the truth.
  * `MX` - one of the best assembled haplotypes in the query has a distance
    greater than `--maxdiff` from the truth. The other has a difference less
    than `--maxdiff` (considered a match).

hap-eval will use the evaluation counts to output a global evaluation on the
matching between the query and truth VCFs using the following formula:
```python
tp = float(summary['MM'])
ff = summary['MX'] + summary['XX'] + summary['??']
fp = summary['FP'] + ff
fn = summary['FN'] + ff

if tp < 1:
  print('precision 0.0 recall 0.0 f1 0.0')
  return 0
print('precision %.4f recall %.4f f1 %.4f' %
  (tp/(tp+fp), tp/(tp+fn), tp/(tp+(fp+fn)/2)))
```

`??`, `MX`, and `XX` events are reported as a combined `FP`/`FN` due to the
presence of variants in both the truth and query VCFs at the evaluated region.