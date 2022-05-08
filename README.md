# Hap-Eval
An open-source VCF comparison engine for structual variant benchmarking

*Note: Hap-Eval is pre-release software and is under active development. If Hap-Eval does
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
  --maxdiff FLOAT       Haplotype difference theshold (default: 0.2)
  --metric STR          Distance metric (default: Levenshtein)
```
