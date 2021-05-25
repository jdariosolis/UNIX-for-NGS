
Contents
1. [Short URL](#1-short-url)
2. [Detecting differential expressed genes by cuffdiff](#2-detecting-differential-expressed-genes-by-cuffdiff)
3. [More visualization (cummeRbund R package)](#3-more-visualization-cummerbund-r-package)
  * [Advanced exercise: Detecting differential expressed genes by edgeR](#advanced-exercise-detecting-differential-expressed-genes-by-edger)
  * [Advanced exercise: Filtering low expressed genes and edgeR exact test](#advanced-exercise-filtering-low-expressed-genes-and-edger-exact-test)
4. [GSEA (Gene Set Enrichment Analysis, topGO)](#4-gsea-gene-set-enrichment-analysis-topgo)
  * [Advanced exercise: topGO analysis for CC and BP](#advanced-exercise-topgo-analysis-for-cc-and-bp)

----

# 1. Short URL

* https://git.io/JU4wp

# 2. Detecting differential expressed genes by cuffdiff

Preparation

1. STAR mapping all
	* w2271_1_hal.fastq
	* w2271_2_hal.fastq
	* w2271_3_hal.fastq
	* w2271_1_lyr.fastq
	* w2271_2_lyr.fastq
	* w2271_3_lyr.fastq
2. Convert SAM file into a sorted BAM file 

Cheating
```sh
$ cp -r /scratch/bio373_2020/data/star_out.tgz ./ 
```

Detecting differentially expressed genes by cuffdiff
```sh
 $ cd /scratch/bio373_2020/your_directory
 $ source /usr/local/ngseq/etc/lmod_profile
 $ module load Tools/Cufflinks/2.2.1
 $ tar zxvf star_out.tgz
 $ mkdir sorted_bams
 $ cp star_out/*/*.sorted.bam sorted_bams/
 $ cd sorted_bams
 $ cuffdiff -p 2 -o cuffdiff_out /scratch/bio373_2020/data/athal_genome.gtf --library-type fr-unstranded -L hal,lyr w2271_1_hal.sorted.bam,w2271_2_hal.sorted.bam,w2271_3_hal.sorted.bam  w2271_1_lyr.sorted.bam,w2271_2_lyr.sorted.bam,w2271_3_lyr.sorted.bam
```

* **It will take 15 minutes**
* Please have a coffee break

Reference
* cuffdiff2: Cole Trapnell, et. al., Nature Biotechnology 31, 46–53, Differential analysis of gene regulation at transcript resolution with RNA-seq
* http://www.nature.com/nbt/journal/v31/n1/abs/nbt.2450.html

Extract differential expression (Pick up the genes (transcripts) that have statistically significant difference)
```sh
 $ cd cuffdiff_out
 $ ls
 $ head gene_exp.diff
 $ grep 'yes' gene_exp.diff
```

* cuffdiff manual: http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/

**Question1**
* How many genes are detected as a differentially expressed gene

**Question2**
* Is HMA4 included in the differentially expressed genes?

**Question3**
* How high are the HMA4 genes expressed?

Plot differential expression
* Copy cuffdiff_out directory to ~/bio373_2020/
```sh
 $ mkdir -p ~/bio373_2020
 $ cp -r cuffdiff_out ~/bio373_2020/
```

Log-in the RStudio
* http://fgcz-genomics.uzh.ch/

Note
* **The science cloud server and fgcz-genomics.uzh.ch share your home directory but they do not share /scratch directory.**


In RStudio 
```R
 > setwd("~/bio373_2020/cuffdiff_out")
 > stat <- read.table("gene_exp.diff", sep = "\t", header = T)
 > kam_hal <- stat[, "value_1"]
 > kam_lyr <- stat[, "value_2"]
 > cols <- ifelse(stat[, "significant"] == "yes", "red", "grey")
 > plot(kam_hal, kam_lyr, col = cols, pch = 19, cex = 0.5, log="xy", xlim=c(0.1,10000), ylim=c(0.1, 10000))
```

![scatter_plot_kam.png](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/scatter_plot_kam.png)


**Question4**
* If you set FDR (q_value) < 0.01, how many differentially expressed genes are detected?

# 3. More visualization (cummeRbund R package)

First please copy */scratch/bio373_2020/data/metal_genes_ids.dat* to your *~/bio373_2020/cuffdiff_out/* on the server
```sh
  $ cp /scratch/bio373_2020/data/metal_genes_ids.txt ~/bio373_2020/cuffdiff_out/
```

* You can see the gene annotation as follows: 
```sh
 $ cd ~/bio373_2020/cuffdiff_out/
 $ /scratch/bio373_2020/data/express_list metal_genes_ids.txt
 FPKM	GeneID	Description
 0.0	AT5G62160	 AtZIP12, ZIP12   zinc transporter 12 precursor   chr5:24960107-24961263 FORWARD LENGTH=1157
 0.0	AT5G44790	 RAN1, HMA7   copper-exporting ATPase / responsive-to-antagonist 1 / copper-transporting ATPase (RAN1)   chr5:18075610-18079988 REVERSE LENGTH=4379
 0.0	AT5G21930	 PAA2, HMA8, ATHMA8   P-type ATPase of Arabidopsis 2   chr5:7243129-7248895 FORWARD LENGTH=5767
 0.0	AT4G37270	 HMA1, ATHMA1   heavy metal atpase 1   chr4:17541771-17546484 REVERSE LENGTH=4714
 0.0	AT4G33520	 PAA1, HMA6   P-type ATP-ase 1   chr4:16118896-16126151 FORWARD LENGTH=7256
 0.0	AT4G30120	 HMA3, ATHMA3   heavy metal atpase 3   chr4:14730401-14733510 REVERSE LENGTH=3110
 0.0	AT4G30110	 HMA2, ATHMA2   heavy metal atpase 2   chr4:14720241-14724584 REVERSE LENGTH=4344
 0.0	AT3G20870	 ZTP29   ZIP metal ion transporter family   chr3:7309279-7312752 REVERSE LENGTH=3474
 0.0	AT3G12750	 ZIP1   zinc transporter 1 precursor   chr3:4051651-4053201 REVERSE LENGTH=1551
 0.0	AT2G46800	 ZAT, ATMTP1, MTP1, ZAT1, ATCDF1   zinc transporter of Arabidopsis thaliana   chr2:19237490-19239405 FORWARD LENGTH=1916
 0.0	AT2G39450	 MTP11, ATMTP11   Cation efflux family protein   chr2:16471632-16473818 REVERSE LENGTH=2187
 0.0	AT2G32270	 ZIP3   zinc transporter 3 precursor   chr2:13704221-13706860 FORWARD LENGTH=2640
 0.0	AT2G30080	 ZIP6, ATZIP6   ZIP metal ion transporter family   chr2:12838730-12840112 REVERSE LENGTH=1383
 0.0	AT2G29410	 MTPB1, ATMTPB1   metal tolerance protein B1   chr2:12616564-12617987 FORWARD LENGTH=1424
 0.0	AT2G19110	 HMA4, ATHMA4   heavy metal atpase 4   chr2:8278881-8286445 FORWARD LENGTH=7565
 0.0	AT2G17790	 VPS35A, ZIP3   VPS35 homolog A   chr2:7733635-7739673 FORWARD LENGTH=6039
 0.0	AT2G04032	 ZIP7   zinc transporter 7 precursor   chr2:1289832-1291383 FORWARD LENGTH=1552
 0.0	AT1G68100	 IAR1   ZIP metal ion transporter family   chr1:25521204-25524175 FORWARD LENGTH=2972
 0.0	AT1G63440	 HMA5   heavy metal atpase 5   chr1:23527655-23531109 FORWARD LENGTH=3455
 0.0	AT1G56590	 ZIP4   Clathrin adaptor complexes medium subunit family protein   chr1:21202020-21204742 REVERSE LENGTH=2723
 0.0	AT1G10970	 ZIP4, ATZIP4   zinc transporter 4 precursor   chr1:3665087-3667139 REVERSE LENGTH=2053
 0.0	AT1G10970	 ZIP4, ATZIP4   zinc transporter 4 precursor   chr1:3665087-3667139 REVERSE LENGTH=2053
```

CummeRbund
* An R package that is designed to aid and simplify the task of analyzing Cufflinks RNA-Seq output.
* http://compbio.mit.edu/cummeRbund/

Install
* https://bioconductor.org/packages/release/bioc/html/cummeRbund.html

Log-in the RStudio
* http://fgcz-genomics.uzh.ch/

```R
 > setwd("~/bio373_2020/cuffdiff_out")
 > library(cummeRbund)

 # read data
 > cuff <- readCufflinks()
 > cuff

 # scatter plot
 > scat <- csScatter(genes(cuff), "hal", "lyr", smooth=T)
 > scat

 # density plot
 > dens <- csDensity(genes(cuff))
 > dens

 # box plot
 > my.genes <- genes(cuff)
 > boxp <- csBoxplot(my.genes)
 > boxp
 > boxp.rep <- csBoxplot(my.genes, replicates = T)
 > boxp.rep

 # clustering
 > dend <- csDendro(my.genes, replicates = T)

 # make a gene set
 > metal_genes <- read.table("metal_genes_ids.txt")
 > metal_genes[,1]
 > my.geneset <- getGenes(cuff, metal_genes[,1])

 # heatmap
 > h <- csHeatmap(my.geneset, cluster = "both")
 > plot(h)
 > h.rep <- csHeatmap(my.geneset, cluster = "both", replicates = T)
 > plot(h.rep)

 # bar plot
 > bar <- expressionBarplot(my.geneset)
 > print(bar)

 # kmean clustering
 > k.means <- csCluster(my.geneset, k = 4)
 > k.means.plot <- csClusterPlot(k.means)
 > print(k.means.plot)
```

## Advanced exercise: Detecting differential expressed genes by edgeR

Let's try if you have time!

edgeR
* One of popular differential expression analysis softwares
* Empirical Analysis of Digital Gene Expression Data in R
* https://bioconductor.org/packages/release/bioc/html/edgeR.html
* User guide https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf


Preparation
* Copy count data file (kam_fpkm.tsv) to your home directory
```sh
 $ cp /scratch/bio373_2020/data/kam_fpkm.tsv ~/bio373_2020/
```

In RStudio
```R
 > setwd("~/bio373_2020") # if you are using fgcz-rstudio.uzh.ch
 > library(edgeR)

 # load count data
 > counts <- read.table("kam_fpkm.tsv", header=T, row.names=1)
 > head(counts)

 # set groups
 > group <- factor(c("kam_hal", "kam_hal", "kam_hal", "kam_lyr", "kam_lyr", "kam_lyr"))

 # according to UserGuide
 > design <- model.matrix(~ group)
 > d <- DGEList(counts = counts, group = group)
 > d

 # calculation normalization factor
 > d <- calcNormFactors(d)
 > d

 # calculate common disparsion
 > d <- estimateGLMCommonDisp(d, design)
 > d

 # calculate trended disparsion
 > d <- estimateGLMTrendedDisp(d, design)
 > d

 # calculate tagwise disparsion
 > d <- estimateGLMTagwiseDisp(d, design)
 > d

 # plot Multi-dimensional scaling
 > plotMDS(d)

 # estimate distribution
 > fit <- glmFit(d, design)

 # estimate likelihood ratio
 > lrt <- glmLRT(fit, coef = 2)
 > topTags(lrt)

 # export the result to a file
 > table <- as.data.frame(topTags(lrt, n = nrow(counts)))
 > write.table(table, file = "edgeR_result.txt", col.names = T, row.names = T, sep = "\t")

 # detect DEGs with FDR < 0.01
 > is.DEG <- as.logical(table$FDR < 0.01)
 > length(which(is.DEG==TRUE))
 > DEG.names <- rownames(table)[is.DEG]
 > plotSmear(lrt, de.tags = DEG.names)
```

![edgeR_kam.png](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/edgeR_kam.png)

Appendix: TMM (Trimmed Mean of M values normalization)
1. calc M and A values
2. remove the genes of upper/lower 30% M values
3. remove the genes of upper/lower 5% A values
4. calc normalization factor with rest genes
5. shift M values

![tmm.png](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/tmm.png)

Reference
* Robinson MD, Oshlack A. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biol. 2010, 11(3):R25. 
* https://www.ncbi.nlm.nih.gov/pubmed/20196867

## Advanced exercise: Filtering low expressed genes and edgeR exact test

* You can skip this exercise, but let's try if you have time
```R
 > data <- read.table("kam_fpkm.tsv", header=T, row.names=1)
 > filter <- apply(data, 1, prod) > 0.0
 > filter_data <- data[filter,]
 > count <- as.matrix(filter_data)
 > head(count)
 > group <- factor(c("kam_hal", "kam_hal", "kam_hal", "kam_lyr", "kam_lyr", "kam_lyr"))

 > d <- DGEList(counts = count, group = group)
 > d <- calcNormFactors(d)
 > d <- estimateCommonDisp(d)
 > d <- estimateTagwiseDisp(d)
 > result <- exactTest(d)
 > topTags(result, 20)
 > table <- as.data.frame(topTags(result, n = nrow(count)))
 > write.table(table, file = "result1.txt", col.names = T, row.names = T, sep = "\t")
 > is.DEG <- as.logical(table$FDR < 0.01)
 > length(which(is.DEG==TRUE))
 > DEG.names <- rownames(table)[is.DEG]
 > plotSmear(result, de.tags = DEG.names)
```

**Question5**
* What happens if you filter the low expressed genes?

**Question6**
* How different is between cuffdiff result and EdgeR result?
* If they are different, why do you think they are different?

# 4. GSEA (Gene Set Enrichment Analysis, topGO)

Reference
* topGO http://www.bioconductor.org/packages/release/bioc/html/topGO.html

Preparation
* Login the server
* Go to ~/bio373_2020/cuffdiff_out directory
* (assuming that you have copied the 'cuffdiff_out'' directory to your home directory)

Extract gene list (top 1000 significantly differentially expressed genes)
```sh
 $ cut -f 1,13 gene_exp.diff|sort -k2 -n|grep "AT" | head -n 1000 > sig_genes_top1000.dat
 $ cp sig_genes_top1000.dat ~/bio373_2020/cuffdiff_out/
```

Prepare GO database
```sh
 $ cp /scratch/bio373_2020/data/AT_GO_all.txt ~/bio373_2020/cuffdiff_out/
```

in R environment
```R
 > setwd("~/bio373_2020/cuffdiff_out") # if you are using fgcz-genomics.uzh.ch
 
 # load library
 > library(topGO)
 
 # Whole set of A. thaliana GO data
 > geneID2GO <- readMappings("AT_GO_all.txt")    
 
 # Check data
 > head(geneID2GO)
 
 # Specify the gene list you want to analyze
 > myInterestedGenes <- names(readMappings("sig_genes_top1000.dat")) 
 > geneNames <- names(geneID2GO)
 > geneList <- factor(as.integer(geneNames %in% myInterestedGenes))
 > names(geneList) <- geneNames
 
 # Check data
 > geneList
 
 # Test which molecular functions (MFs) are enriched using normal Fishers exact test
 > GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
 > GOdata
 
 > resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
 > resultFisher
 
 > resultTable  <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 20)
 > resultTable
 
 > showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'all')
```

![gotree.png](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/gotree.png)

**Question7**
* Which GO terms (IDs) are related to heavy metal ion transfer?

## Advanced exercise: topGO analysis for CC and BP

* Let's try if you have time

* Run topGO analysis with CC (Cellular Component) and BP (Biological Process)
* Run topGO analysis with top 1,000 high expressed genes in halleri-subgenome
