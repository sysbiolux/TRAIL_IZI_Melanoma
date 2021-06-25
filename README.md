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

### Tools used

|   software and algorithms        |    			source           |    			product number/identifier       |
|:--------------------------------:|:-------------------------:|:----------------------------------------:|
|        			AdapterRemoval       |    			https://github.com/MikkelSchubert/adapterremoval  |  v2.3.1     |
|     			apeglm                 |     			R bioconductor                    |           			v1.10.0    |
|       			ComplexHeatmap       |     			https://github.com/jokergoo/ComplexHeatmaps    |    v2.7.8.1000|
|    			DESeq2                    |       			R bioconductor                    |         	v1.128.1      |
|   			FALCON             		   |   			De Landtsheer et al. 2017, https://github.com/sysbiolux/FALCON>  | |
|     			GeneWalk               |  Ietswaart et al. 2021, https://churchman.med.harvard.edu/genewalk    | 	v1.5.1   |
|  		ggplot2                  |        	R CRAN |			v3.3.2 	|  
|  			ggforce                |         	R CRAN | 		v0.3.3      |
|  			multiqc                |  https://multiqc.info/    | 		v1.9      |
|  			python                |  https://python.org/    | 		v3.8.2      |
|		R Project for Statistical Computing |   		  https://www.r-project.org   |    v4.0.0 | 
|   			RStudio            |    	https://www.rstudio.com/     | v1.0.143   |  
|   Rsubread                 | R bioconductor | v2.2.2  |
| 			SAMtools             | https://www.htslib.org   |   v1.10      |
|		Snakemake                |   https://snakemake.github.io  |  v5.20.1    |
|    STAR                    |  https://github.com/alexdobin/STAR | v2.7.4a |
|	Tidyverse                  | R CRAN, http://tidyverse.org/ |  v1.1.1   |

