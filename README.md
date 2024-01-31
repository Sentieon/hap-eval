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

## Interpreting the output
Upon successful evaluation, hap-eval will print output similar to the following:
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