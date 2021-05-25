
Contents

1. [Short URL](#1-short-url)
2. [Expression Level (Count mapped reads using Cufflinks)](#2-expression-level-count-mapped-reads-using-cufflinks)
3. [Cufflinks other samples](#3-cufflinks-other-samples)
4. [Make a FPKM table](#4-make-a-fpkm-table)
 * [Advanced exercise1 (Kallisto quantification)](#advanced-exercise1-kallisto-quantification)
 * [Advanced exercise2 (Scatter plot and correlation matrix in R)](#advanced-exercise2-scatter-plot-and-correlation-matrix-in-r)
5. [Hierarchical Clustering](#5-hierarchical-clustering)
6. [Log transform (scaling)](#6-log-transform-scaling)
7. [MDS (Multi\-Dimensional Scaling) plot (by ezRun)](#7-mds-multi-dimensional-scaling-plot-by-ezrun)
----

# 1. Short URL

* https://git.io/JU4wx

# 2. Expression Level (Count mapped reads using Cufflinks)

Reference
* cufflinks http://cufflinks.cbcb.umd.edu/

Count reads as FPKM
* Login the server
* Go to your working directory
```sh
 $ cd /scratch/bio373_2020/your_directory
 $ source /usr/local/ngseq/etc/lmod_profile
 $ module load Tools/Cufflinks/2.2.1
 $ cufflinks -p 2 -o cufflinks_out_w2271_1_hal -G /scratch/bio373_2020/data/athal_genome.gtf --library-type fr-unstranded star_out_w2271_1_hal/w2271_1_hal.sorted.bam
```

* It will take 3 minutes

Options
* -p: number of processes
* -o: output directory
* -G: gene annotation file
* --library-type: direction of the library

Note
* default: --library-type fr-unstranded (for usual paired-end reads)
* Illumina sense: --library-type fr-secondstrand (for strand specific single-end reads)
* Illumina antisense: --library-type fr-firststrand (for strand specific single-end reads)

Reference
* Cufflinks: http://cole-trapnell-lab.github.io/cufflinks/
* GTF Gene Annotation Format: http://mblab.wustl.edu/GTF22.html

Check expression level (FPKM)
```sh
 $ cd cufflinks_out_w2271_1_hal
 $ ls
 $ head genes.fpkm_tracking
 $ head isoforms.fpkm_tracking
```

* .fpkm_tracking is a tsv (Tab Separated Values) file.

View the expression level on Excel
* log out from the server
```sh
 $ scp your_account_name@172.23.xx.xx:/scratch/bio373_2020/your_directory/cufflinks_out_w2271_1_hal/genes.fpkm_tracking ./genes.fpkm_tracking.xls
```
 
* Open genes.fpkm_tracking on Excel

![excel_genes.png](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/excel_genes.png)

References
* Cufflinks: http://cole-trapnell-lab.github.io/cufflinks/
* Cufflinks output: http://cole-trapnell-lab.github.io/cufflinks/file_formats/index.html#output-formats-used-in-the-cufflinks-suite
* Cufflinks paper: Trapnell C, et.al, Nature Biotech., 2010, http://www.ncbi.nlm.nih.gov/pubmed/20436464

Show only Gene ID and FPKM
* Login server again
* Go back to cufflinks output directory
```sh
 $ cut -f 1,10 genes.fpkm_tracking | head
```

Note
* **Please think where the cufflinks output directory is, and please think always where you are working on**

Sort data by FPKM
```sh
 $ sort -k 10 -nr genes.fpkm_tracking | head
```

Extract top 20 expression data
```sh
 $ cut -f 1,10 genes.fpkm_tracking|sort -k 2 -nr | head -n 20
```

Commands
* *cut*: show selected columns, -f[column number]
* *sort*: sort output, -k [key column] -n numerical order -r reverse
* *head*: show the first some lines, -n [number of lines]

Add annotation information from BLAST dataset of TAIR10 (Arabidopsis thaliana) to the expression data
* Login server
* Go to your work directory
```sh
 $ cd /scratch/bio373_2020/your_directory
```

* Run the script
```sh
 $ source /usr/local/ngseq/etc/lmod_profile
 $ module load Dev/Ruby/2.6.6
 $ /scratch/bio373_2020/data/express_list cufflinks_out_w2271_1_hal/isoforms.fpkm_tracking
 (You can use more (less) command if you see the result one screen by screen)
 $ /scratch/bio373_2020/data/express_list cufflinks_out_w2271_1_hal/isoforms.fpkm_tracking | more
```

/scratch/bio373_2020/data/express_list (Ruby script)
* You do not have to understand the code, but you should understand what this script does
```Ruby
 #!/usr/bin/env ruby
 # encoding: utf-8

 unless isoforms_fpkm = ARGV[0]
   puts "Usage:\n ruby #{__FILE__} [isoforms.fpkm_tracking]"
   exit
 end

 require 'csv'
 tair_fa = "/scratch/bio373_2020/data/TAIR10_annotation.txt"
 trans_list = {}
 genes_list = {}
 File.readlines(tair_fa).each do |line|
   if line =~ />/
     atid, *others = line.chomp.split('|')
     trans_list[atid.gsub(/>/,'').strip] = others.join(" ").gsub(/Symbols:\s+/, '')
     genes_list[atid.gsub(/>/,'').strip.split('.').first] ||= others.join(" ").gsub(/Symbols:\s+/, '')
   end
 end

 list = []
 File.readlines(isoforms_fpkm).each do |line|
   x = line.split
   gid = x[0]
   if gid =~ /\./
     list << [x[9].to_f, gid, trans_list[gid]]
   else
     list << [x[9].to_f, gid, genes_list[gid]]
   end
 end

 print "FPKM\tGeneID\tDescription\n"
 list.sort.reverse.each do |line|
   print line.join("\t"), "\n"
 end
```

**Question1**
* What is the highest expressed gene?

**Question2**
* How is the HMA4 gene expressed?

# 3. Cufflinks other samples

1. Please run cufflinks on the other samples, too.
	* w2271_1_lyr.fastq
	* w2271_2_hal.fastq
	* w2271_2_lyr.fastq
	* w2271_3_hal.fastq
	* w2271_3_lyr.fastq
2. Check the differences between them

Hint
* You need the sorted BAM files for the samples before calling *cufflinks*
* You can make a shell script

cufflinks_batsh.sh
```bash
source /usr/local/ngseq/etc/lmod_profile
module load Tools/Cufflinks/2.2.1
cufflinks -p 2 -o cufflinks_out_w2271_1_lyr -G /scratch/bio373_2020/data/athal_genome.gtf --library-type fr-unstranded star_out_w2271_1_lyr/w2271_1_lyr.sorted.bam
cufflinks -p 2 -o cufflinks_out_w2271_2_hal -G /scratch/bio373_2020/data/athal_genome.gtf --library-type fr-unstranded star_out_w2271_2_hal/w2271_2_hal.sorted.bam
cufflinks -p 2 -o cufflinks_out_w2271_2_lyr -G /scratch/bio373_2020/data/athal_genome.gtf --library-type fr-unstranded star_out_w2271_2_lyr/w2271_2_lyr.sorted.bam
cufflinks -p 2 -o cufflinks_out_w2271_3_hal -G /scratch/bio373_2020/data/athal_genome.gtf --library-type fr-unstranded star_out_w2271_3_hal/w2271_3_hal.sorted.bam
cufflinks -p 2 -o cufflinks_out_w2271_3_lyr -G /scratch/bio373_2020/data/athal_genome.gtf --library-type fr-unstranded star_out_w2271_3_lyr/w2271_3_lyr.sorted.bam
```

# 4. Make a FPKM table

* Make sure you have done all the cufflinks jobs before running the following commands
* Go to your working directory on the server
```sh
 $ cut -f 1,10 cufflinks_out_w2271_1_hal/genes.fpkm_tracking | sed 's/FPKM/H1/g' | sort -r > w2271_H1_fpkm.dat
 $ cut -f 1,10 cufflinks_out_w2271_1_lyr/genes.fpkm_tracking | sed 's/FPKM/L1/g' | sort -r | cut -f 2 > w2271_L1_fpkm.dat
 $ cut -f 1,10 cufflinks_out_w2271_2_hal/genes.fpkm_tracking | sed 's/FPKM/H2/g' | sort -r | cut -f 2 > w2271_H2_fpkm.dat
 $ cut -f 1,10 cufflinks_out_w2271_2_lyr/genes.fpkm_tracking | sed 's/FPKM/L2/g' | sort -r | cut -f 2 > w2271_L2_fpkm.dat
 $ paste w2271_*.dat > kam_fpkm.tsv
```

Note
* Let's think about what each command above does
* *paste*: merge lines of files with tab space

**Question3**
* How different are Zinc-binding protein homeologs expressed?

## Advanced exercise1 (Kallisto quantification)

* Let's try if you have time but you do not have to do

Kallisto quantification
* Login the server
* Go to your working directory
```sh
 $ source /usr/local/ngseq/etc/lmod_profile
 $ module load Aligner/kallisto/0.46.1
 $ mkdir athal_cdna
 $ kallisto index -i athal_cdna/athal_cdna.kallisto_index /scratch/bio373_2020/data/athal_cdna.fa
 $ mkdir kallisto_out
 $ kallisto quant -i references/athal_cdna.kallisto_index -o kallisto_out/H1 --single -l 100 -s 1 -t 2 --bias -b 10 kam_samples/w2271_1_hal.fastq
 $ kallisto quant -i references/athal_cdna.kallisto_index -o kallisto_out/H2 --single -l 100 -s 1 -t 2 --bias -b 10 kam_samples/w2271_2_hal.fastq
 $ kallisto quant -i references/athal_cdna.kallisto_index -o kallisto_out/L1 --single -l 100 -s 1 -t 2 --bias -b 10 kam_samples/w2271_1_lyr.fastq
 $ kallisto quant -i references/athal_cdna.kallisto_index -o kallisto_out/L2 --single -l 100 -s 1 -t 2 --bias -b 10 kam_samples/w2271_2_lyr.fastq
```

Cheating
```sh
  $ source /usr/local/ngseq/etc/lmod_profile
  $ module load Aligner/kallisto/0.46.1
  $ cp /scratch/bio373_2020/data/kallisto_batch.sh .
  $ sh kallisto_batch.sh
```

* Let's compare the calculation time with STAR+Cufflinks
* Let's compare the gene expressions with the ones from Cufflinks

Note
* *time*: calculate command execution time
  e.g. $ time kallisto quant -i references/athal_cdna.kallisto_index -o kallisto_out/H1 --single -l 100 -s 1 -t 2 --bias -b 10 kam_samples/w2271_1_hal.fastq

## Advanced exercise2 (Scatter plot and correlation matrix in R)

* Let's try if you have time but you do not have to do

Online RStudio
* https://fgcz-genomics.uzh.ch/
* ID, PASS: B-Fabric account

Note
* **Home directory is shared between RStudio and the server, but /scratch is not shared with RStudio**

Copy the count data file to your home directory on the server, 172.23.xx.xx.
```sh
 $ cd (go to your home directory)
 $ mkdir bio373_2020
 $ cp /scratch/bio373_2020/your_directory/kam_fpkm.tsv bio373_2020/
 $ cd ~/bio373_2020
 $ ls
```

Note
* This is Tab-Separated-Value, text file.
* This is the result coming from the cufflinks command for all samples.
* RSudio server at FGCZ is available once you get a B-Fabric account.

Go to RStudio server
* http://fgcz-genomics.uzh.ch/
* and start a new project
	* [File]-[New Project], select *Existing Directory*, and select the *bio373_2020/* folder.

Load data (in R environment)
```R
 > data <- read.table("kam_fpkm.tsv", header=T, row.names=1)
 > head(data)
```

Correlation matrix
```R
 > cor(data)
```

Scatter plot
```R
 > plot(data$H1, data$H2)
 > plot(data$H1, data$L1)
 > plot(data$H1, data$H2, log="xy", xlim=c(0.001,100000), ylim=c(0.001, 100000))
 > plot(data$H1, data$L1, log="xy", xlim=c(0.001,100000), ylim=c(0.001, 100000))
```

Note
* plot(x-axis data, y-axis data)
* log="xy": scaling logarithmically in x- and y-axis
* xlim: plot region of x-axis
* ylim: plot region of y-axis

Boxplot
```R
 > data <- read.table("kam_fpkm.tsv", header=T, row.names=1)
 > head(data)
 > hal <- t(data["AT5G16980",1:2])
 > lyr <- t(data["AT5G16980",3:4])
 > tb <- cbind(hal, lyr)
 > colnames(tb) <- c("H", "L")
 > boxplot(tb, ylab="FPKM", main="AT5G16980")
```

![box.png](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/box.png)

Scatter plot matrix
```R
 > pairs(data)
```

Plot scatter matrix and correlation matrix together
```R
 > library(psych)
 > pairs.panels(data)
```

![corr.png](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/corr.png)

# 5. Hierarchical Clustering

Load data (in R environment)
```R
 > data <- read.table("kam_fpkm.tsv", header=T, row.names=1)
 > head(data)
 > distance <- dist(t(data))
 > distance
 > clustering <- hclust(distance)
 > plot(clustering)
```

![clustering1.png](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/clustering1.png)

# 6. Log transform (scaling)

```R
 > data.log2 <- log2(data+1)
 > distance <- dist(t(data.log2))
 > clustering <- hclust(distance)
 > plot(clustering)
```

**Question4**
* What is the difference with and without scaling?

# 7. MDS (Multi-Dimensional Scaling) plot (by ezRun)

```R
 > library(ezRun)
 > ezMdsPlot(data, sampleColors = c("red", "blue", "green", "yellow"), main="MDS plot")
```

![mds.png](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/mds.png)

Reference
* ezRun: https://github.com/uzh/ezRun/
