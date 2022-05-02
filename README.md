# Hap-Eval
An open-source VCF comparison engine for structual variant benchmarking

*Note: Hap-Eval is pre-release software and is under active development. If Hap-Eval does
not work with your VCF, please file an issue* 

Hap-Eval requires access to the vcflib library contained in the Sentieon
software package. The library is located under
$SENTIEON_INSTALL_DIR/lib/python/sentieon in the Sentieon software package.
Calling `hap_eval.py` with `sentieon pyexec` properly loads the vcflib library
into the PYTHONPATH.

## Usage
```
usage: sentieon pyexec hap_eval.py [-h] -r FASTA -b VCF -c VCF [-i BED]
                                   [-t INT] [--maxdist INT] [--minsize INT]
                                   [--maxdiff FLOAT]

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
  --maxdist INT         Maximum distance to cluster variants (default: 1000)
  --minsize INT         Minimum size of variants to consider (default: 50)
  --maxdiff FLOAT       Haplotype difference theshold (default: 0.2)
```
