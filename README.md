## 200602 - TRAIL/IZI in Melanoma

Principal investigators: Pr.
Thomas Sauter (Luxembourg) / Pr.
Dagmar Kulms (Dresden)

### Data

Experimenter: Sebastian Schindler (Dresden)

Sequencing on total RNA: Rashi Halder (Luxembourg)

24 FASTQ generated (deposited in [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/) E-MTAB-10669)

```{bash}
10-1-Lasse_S10_R1_001.fastq.gz
10-3-Lasse_S16_R1_001.fastq.gz
11-1-Lasse_S11_R1_001.fastq.gz
11-3-Lasse_S17_R1_001.fastq.gz
1-1-Lasse_S1_R1_001.fastq.gz
12-1-Lasse_S12_R1_001.fastq.gz
12-3-Lasse_S18_R1_001.fastq.gz
1-3-Lasse_S7_R1_001.fastq.gz
2-1-Lasse_S2_R1_001.fastq.gz
2-3-Lasse_S8_R1_001.fastq.gz
3-1-Lasse_S3_R1_001.fastq.gz
3-3-Lasse_S9_R1_001.fastq.gz
4-1-Lasse_S4_R1_001.fastq.gz
4-3-Lasse_S10_R1_001.fastq.gz
5-1-Lasse_S5_R1_001.fastq.gz
5-3-Lasse_S11_R1_001.fastq.gz
6-1-Lasse_S6_R1_001.fastq.gz
6-3-Lasse_S12_R1_001.fastq.gz
7-1-Lasse_S7_R1_001.fastq.gz
7-3-Lasse_S13_R1_001.fastq.gz
8-1-Lasse_S8_R1_001.fastq.gz
8-3-Lasse_S14_R1_001.fastq.gz
9-1-Lasse_S9_R1_001.fastq.gz
9-3-Lasse_S15_R1_001.fastq.gz
```

### Snakemake

The quality controls, read trimming, feature counting were performed using [Snakemake](https://snakemake.github.io/) and the tailor-made [RNA-seq template](https://git-r3lab.uni.lu/aurelien.ginolhac/snakemake-rna-seq/-/tree/master).

The 3 config files are provided in the `config` folder:

-   `config.yaml` main parameters
-   `samples.tsv` sample annotations
-   `units.tsv` link samples and FASTQ files

#### Outputs

Main output is the `DESeqDataSet` object from [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) saved a binary file \`deseq2/all.rds\` after successful completion of the `snakemake` pipeline.

### Differential gene analysis and plots

See the file `TRAIL_IZI.Rmd` , routines are loaded from `R/utils.R`
