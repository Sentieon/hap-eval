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
  --tr_bed BED          Tandem repeat BED file
  -t INT, --thread_count INT
                        Number of threads
  --base_out VCF        Annotated baseline vcf file
  --comp_out VCF        Annotated comparison vcf file
  --minsize INT         Minimum size of variants to consider (default: 50)
  --maxdiff FLOAT       Haplotype difference threshold (default: 0.2)
  --metric STR          Distance metric (default: Levenshtein)
```

### Clustering, `--tr_bed`, and `--minsize`

hap-eval evaluates groups (clusters) of variants that are nearby on the reference genome. Variants are grouped together if the distance between two variants is small relative to the variant size. If the variant group extends into a known tandem repeat (from the `--tr_bed` file), the group will be extended to the end of the tandem repeat. The entire group will be evaluated if the largest possible haplotype that can be constructed from the group has an insertion or deletion variant larger than `--minsize`.

Using a `--tr_bed` bed file highly recommended as different structural variant callers may place structural variants at different locations along a tandem repeat. A tandem repeat bed file for hg38 is provided in the [/data](/data) directory.

Some examples are shown below. For the purpose of demonstration, these examples use a baseline VCF (`--base`), but an empty comparison VCF `--comp`.

#### Basic grouping

The baseline VCF has the following variants:
```txt
chr1    1948926 .       T       C       30      .       TRF;LCR=0       GT:AD   1|0:1,1
chr1    1948934 .       T       TCCCTCCCTTCTTTCCTTCCCTTTCCCTCCCTCCCTTCCTTCCTCTTTCCTTCCTTCCTTTCCCTCCCTTACTCCTTCCTTCCTTCCCTTCCCCTTCCTTCTTCCTTCTCTC        30      .       TRF;TRFdiff=0;TRFrepeat=CCTTCCTTCCTTC;TRFovl=1;TRFstart=1948876;TRFend=1949289;TRFperiod=13;TRFcopies=32;TRFscore=633;TRFentropy=1.04;SVTYPE=INS;SVLEN=111;RM_score=35;RM_repeat=(TCCC)N;RM_clsfam=Simple_repeat;LCR=0.529985;REMAP=interspersed        GT:AD   1|0:1,1
chr1    1948942 .       T       C       30      .       TRF;LCR=0       GT:AD   1|0:1,1
chr1    1948947 .       T       C       30      .       TRF;LCR=0       GT:AD   1|0:1,1
```

These variants are grouped together as they are relatively close to the insertion allele. They are evaluated as a single region by hap-eval, and contributed one `FN` to the overall evaluation (as there is no comparison VCF).

```txt
FN chr1:1948826-1949047 4 0 222 (111,0) (0,0) 1.0
```

#### Distant variants are not grouped

The baseline VCF has the following variants:
```txt
chr1    1981820 .       A       ACCGCGGACAGACACGGGGGCACGCAGGACACCCAGCCGCGGACAGACACGGGGGCACGCAGGACACCCAG 30      .       TRF;TRFdiff=2;TRFrepeat=ACACCCAGCCGCGGACAGACACGGGGGCACGCAGG;TRFovl=1;TRFstart=1981391;TRFend=1982440;TRFperiod=35;TRFcopies=31.9;TRFscore=2263;TRFentropy=1.71;TRFsim=1;SVTYPE=INS;SVLEN=70;LCR=0.785201;REMAP=tandem   GT:AD   1|1:0,2
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
```

These variants are split into two groups, as the distance between the variant at `chr1:1981820` and the variant at `chr1:1993704` is relatively large and the insertions are relatively small. The intervening SNVs are not included in either group.

```txt
FN chr1:1981720-1982541 1 0 822 (70,70) (0,0) 1.0
FN chr1:1993604-1993804 1 0 201 (147,147) (0,0) 1.0
```

#### Multiple variants exceeding `--minsize`

The baseline VCF has the following variants:
```txt
chr1    16166161        .       G       GCCATCCATCCAGCCATCCAT   30      .       TRF;TRFdiff=0;TRFrepeat=TCCA;TRFovl=1;TRFstart=16166148;TRFend=16166283;TRFperiod=4;TRFcopies=33.8;TRFscore=325;TRFentropy=1.65;SVTYPE=INS;SVLEN=20;LCR=0.890708;REMAP=interspersed     GT:AD   1|0:1,1
chr1    16166181        .       T       G       30      .       TRF;LCR=0       GT:AD   1|0:0,1
chr1    16166181        .       TCCATCCATCCATCCATCCAG   T       30      .       TRF;TRFdiff=-2.5;TRFrepeat=TCCCTCCT;TRFovl=1;TRFstart=16166164;TRFend=16166283;TRFperiod=8;TRFcopies=12.6;TRFscore=73;TRFentropy=1.6;SVTYPE=DEL;SVLEN=20;LCR=0.852383;REMAP=interspersed        GT:AD   0|1:0,1
chr1    16166229        .       T       TCCATCCATCCATCCATCCATCCATCCAGCCATCCATCCATCCATCCAG       30      .
       TRF;TRFdiff=0;TRFrepeat=TCCA;TRFovl=1;TRFstart=16166148;TRFend=16166311;TRFperiod=4;TRFcopies=40.8;TRFscore=399;TRFentropy=1.65;SVTYPE=INS;SVLEN=48;RM_score=48;RM_repeat=(TCCA)N;RM_clsfam=Simple_repeat;LCR=0.83682;REMAP=tandem      GT:AD   1|0:1,1
```

Individually, each of these variants is smaller than `--minsize`. However, they can be combined together into a haplotype larger than `--minsize`. Accordingly, they are evaluated as a group by hap-eval.

```txt
FN chr1:16166061-16166329 4 0 269 (68,-20) (0,0) 1.0
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