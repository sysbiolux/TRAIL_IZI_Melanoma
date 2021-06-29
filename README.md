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

-   [`config.yaml`](config/config.yaml) main parameters
-   [`samples.tsv`](config/samples.tsv) sample annotations
-   [`units.tsv`](config/units.tsv) link samples and FASTQ files

#### Outputs

Main output is the `DESeqDataSet` object from [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) saved a binary file \`deseq2/all.rds\` after successful completion of the `snakemake` pipeline.

### Differential gene analysis and plots

See the file [`TRAIL_IZI.Rmd`](TRAIL_IZI.Rmd), routines are loaded from [`R/utils.R`](R/utils.R)


### GeneWalk

Performed on 1,248 genes using the BASH script [`run_genewalk.sh`](run_genewalk.sh) on the https://hpc.uni.lu `iris` machine.

The input file [`same_pdirection.txt`](same_pdirection.txt) is also part of this repository.

### Tools used

The experiments presented in this paper were carried out using the HPC facilities of the University of Luxembourg (Varrette et al. 2014, https://hpc.uni.lu).

Varrette S., Bouvry P., Cartiaux H. and Georgatos F. Management of an Academic HPC Cluster: The UL Experience. _Proc. of the 2014 Intl. Conf. on High Performance Computing & Simulation_. **HPCS 2014**. 959-967.

|   software and algorithms        |    			source           |    			product number/identifier       |
|:--------------------------------:|:-------------------------:|:----------------------------------------:|
|        			AdapterRemoval       |    			https://github.com/MikkelSchubert/adapterremoval  |  v2.3.1     |
|     			apeglm                 |     			R bioconductor                    |           			v1.10.0    |
|       			ComplexHeatmap       |     			https://github.com/jokergoo/ComplexHeatmaps    |    v2.7.8.1000|
|    			DESeq2                    |       			R bioconductor                    |         	v1.128.1      |
|   			FALCON             		   |   			De Landtsheer et al. 2017, https://github.com/sysbiolux/FALCON  | |
|     			GeneWalk               |  Ietswaart et al. 2021, https://churchman.med.harvard.edu/genewalk    | 	v1.5.1   |
|  		ggplot2                  |        	R CRAN |			v3.3.2 	|  
|  		ggVennDiagram            |        	R CRAN |			v1.1 	|  
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


<details><summary>Tools citation</summary>
<p>
- Constantin Ahlmann-Eltze, Peter Hickey and Hervé Pagès (2021). MatrixGenerics: S4 Generic Smmary Statistic Functions that Operate on Matrix-Like Objects. R package version 1.4.0. htps://bioconductor.org/packages/MatrixGenerics
- Gu, Z. (2016) Complex heatmaps reveal patterns and correlations in multidimensional genomic dta. Bioinformatics.
- H. Pagès, M. Lawrence and P. Aboyoun (2021). S4Vectors: Foundation of vector-like and list-like cntainers in Bioconductor. R package version 0.30.0. https://bioconductor.org/packages/S4Vectors
- H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
- Hadley Wickham (2019). stringr: Simple, Consistent Wrappers for Common String Operations. R pckage version 1.4.0. https://CRAN.R-project.org/package=stringr
- Hadley Wickham (2021). forcats: Tools for Working with Categorical Variables (Factors). R pckage version 0.5.1. https://CRAN.R-project.org/package=forcats
- Hadley Wickham (2021). tidyr: Tidy Messy Data. R package version 1.1.3. htps://CRAN.R-project.org/package=tidyr
- Hadley Wickham and Jim Hester (2020). readr: Read Rectangular Text Data. R package version 14.0. https://CRAN.R-project.org/package=readr
- Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021). dplyr: A Grammar of Data Mnipulation. R package version 1.0.7. https://CRAN.R-project.org/package=dplyr
- Henrik Bengtsson (2021). matrixStats: Functions that Apply to Rows and Columns of Matrices (and t Vectors). R package version 0.59.0. https://CRAN.R-project.org/package=matrixStats
- Jim Hester and Hadley Wickham (2021). vroom: Read and Write Rectangular Text Data Quickly. R pckage version 1.5.1. https://CRAN.R-project.org/package=vroom
- Kamil Slowikowski (2021). ggrepel: Automatically Position Non-Overlapping Text Labels with 'gplot2'. R package version 0.9.1. https://CRAN.R-project.org/package=ggrepel
- Kirill Müller and Hadley Wickham (2021). tibble: Simple Data Frames. R package version 3.1.2. htps://CRAN.R-project.org/package=tibble
- Lawrence M, Huber W, Pag`es H, Aboyoun P, Carlson M, et al. (2013) Software for Computing and Anotating Genomic Ranges. PLoS Comput Biol 9(8): e1003118. doi:10.1371/journal.pcbi.1003118
- Lawrence M, Huber W, Pag`es H, Aboyoun P, Carlson M, et al. (2013) Software for Computing and Annotating Genomic Ranges. PLoS Comput Biol 9(8): e1003118. doi:10.1371/journal.pcbi.1003118
- Lionel Henry and Hadley Wickham (2020). purrr: Functional Programming Tools. R package version 03.4. https://CRAN.R-project.org/package=purrr
- Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq dta with DESeq2 Genome Biology 15(12):550 (2014)
- Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2021). SummarizedExperiment: SmmarizedExperiment container. R package version 1.22.0. htps://bioconductor.org/packages/SummarizedExperiment
- Orchestrating high-throughput genomic analysis with Bioconductor. W. Huber, V.J. Carey, R. Gntleman, ..., M. Morgan Nature Methods, 2015:12, 115.
- Orchestrating high-throughput genomic analysis with Bioconductor. W. Huber, V.J. Carey, R. Gntleman, ..., M. Morgan Nature Methods, 2015:12, 115.
- R Core Team (2021). R: A language and environment for statistical computing. R Foundation for Satistical Computing, Vienna, Austria. URL https://www.R-project.org/.
- Sonali Arora, Martin Morgan, Marc Carlson and H. Pagès (2021). GenomeInfoDb: Utilities for mnipulating chromosome names, including modifying them to follow a particular naming style. R pckage version 1.28.0. https://bioconductor.org/packages/GenomeInfoDb
- Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686

</p>
</details>