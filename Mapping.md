
Contents
1. [Short URL](#1-short-url)
2. [Quality Check](#2-quality-check)
3. [Setup sample RNAseq data in your working directory](#3-setup-sample-rnaseq-data-in-your-working-directory)
4. [Mapping by STAR](#4-mapping-by-star)
5. [Check mapping result on IGV](#5-check-mapping-result-on-igv)
   * [Let's Try\!\!](#lets-try)

----
# 1. Short URL

https://git.io/JU4wN

# 2. Quality Check

Login
```
$ ssh your_account_name@172.23.xx.xx
```

Please select your assigned server:
1. 172.23.30.6 (fgcz-kl-003)
2. 172.23.30.200 (fgcz-kl-004)
3. 172.23.32.84 (fgcz-kl-005)


Change your working directory 
```sh
$ cd /scratch/bio373_2020/your_directory
```

Note
* **NOT** *your_directory*, but *the directory that you made in the previous exercise*

# 3. Setup sample RNAseq data in your working directory

Copy the following file in your working directory (/scratch/bio373_2020/*your_directory*)
```sh
$ hostname
  (make sure you are working on the server)
$ pwd
  (make sure you are working in the working directory, /scratch/bio373_2020/your_directory)
$ cp /scratch/bio373_2020/data/akam_samples2.tgz /scratch/bio373_2020/your_directory/
$ ls
```

Note (again)
* Please **NOT** type like **your_directory** but the directory name that you made

Extract (decompress) it under your directory
```sh
 $ ls
  (make sure there is certainly akam_samples2.tgz in your current working directory)
 $ tar zxvf akam_samples2.tgz
 $ ls
```

Check fastq file
```sh
 $ ls
  (make sure there is certainly akam_samples2/ directory created)
 $ cd akam_samples2/
 $ head w2271_1_hal.fastq
```

Note
* There should be 6 samples
* naming rule: 
* e.g: w2271_1_hal: Arabidopsis kamchatica accession "w2271", biological replicate No. 1, *A. halleri* subgenome 
* These samples were extracted from Leaf tissue RNA with Zinc treated after 48 hours 
* 1 million Illumina reads are randomly sampled in each sample (for the quick process in this course practice)

# 4. Mapping by STAR

References
* Alexander Dobin, et.al., Bioinformatics, 2012, http://bioinformatics.oxfordjournals.org/content/early/2012/10/25/bioinformatics.bts635
* Source Code https://github.com/alexdobin/STAR

Make an index file for STAR
```sh
 $ cd /scratch/bio373_2020/your_directory
 $ source /usr/local/ngseq/etc/lmod_profile
 $ module load Aligner/STAR/2.7.4a
 $ mkdir athal_genome
 $ cp /scratch/bio373_2020/data/athal_genome.fa athal_genome/
 $ STAR --runMode genomeGenerate --genomeDir ./athal_genome --genomeFastaFiles athal_genome/athal_genome.fa
```

Note
* It will take 5 minutes

Command options
* --runMode genomeGenerate: make an index file
* --genomeDir: output directory where a reference sequence is located
* --genomeFastaFiles: the reference genome FASTA file

Run STAR mapping
```sh
 $ mkdir star_out_w2271_1_hal
 $ STAR --genomeDir ./athal_genome --runThreadN 2 --genomeLoad NoSharedMemory --outSAMstrandField intronMotif --outFileNamePrefix star_out_w2271_1_hal/ --readFilesIn ./akam_samples2/w2271_1_hal.fastq
```

Note
* It will take 5 minutes

Command options
* --genomeDir: directory where a reference sequence is located
* --runThreadN: the number of threads (cores) to be used for the calculation
* --outFileNamePrefix: output directory
* --readFilesin: input fastq file
* --genomeLoadNoSharedMemory: do not use shared memory space. Sometimes error happens if you map a big size data without this option.
* --outSAMstrandField intronMotif: this is needed for cufflinks (quantification)

# 5. Check mapping result on IGV

Reference
* Integrated Genome Viewer IGV: http://www.broadinstitute.org/igv/
* IGV User Guide: http://software.broadinstitute.org/software/igv/userguide

Go to the STAR result directory
```sh
 $ cd star_out_w2271_1_hal
 $ ls
```

Check mapping results
```sh
 $ cat Log.final.out
```

Check the result SAM file
```sh
 $ head Aligned.out.sam
 $ grep -v '@SQ' Aligned.out.sam|head
```

References
* SAM format: http://samtools.sourceforge.net/SAMv1.pdf
* BAM is a compressed binary file of SAM
* samtools: http://samtools.sourceforge.net/

Convert the SAM to a sorted BAM file (it is necessary to look at the alignment on IGV)
```sh
 $ cd /scratch/bio373_2020/your_directory/star_out_w2271_1_hal/
 (You do not have to execute the command above, if you are working under the star_out_w2271_1_hal/ directory)
 $ module load Tools/samtools/1.10
 $ samtools sort -O BAM -o w2271_1_hal.sorted.bam Aligned.out.sam
```

Command options
* -O: output file format
* -o: output file name
 
Make BAM index file (it is necessary to look at the alignment on IGV)
```sh
 $ samtools index w2271_1_hal.sorted.bam
 $ ls 
 (check if w2271_1_hal.sorted.bam.bai is created)
```

Download
* **Logout from server**

```sh
 $ hostname
   (make sure you are working on Mac)
 $ cd ~/Desktop
 $ scp your_account_name@172.23.xx.xx:/scratch/bio373_2020/your_directory/star_out_w2271_1_hal/*sorted* ./
 $ ls (check if *.bam and *.bai files are downloaded)
```

Note
* **Please confirm first that you are logged-out from the server**
* The symbol * is available for matching any characters (it is called *wildcard*). 
* In this case, for example, this command copies all the files matching with *sorted* in the directory.

View Mapping Result
* start IGV (igv.jar)
* Select "A. thaliana (TAIR10)" from the left top selector

![igv_select_genome.png](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/igv_select_genome.png)

* [File]-[Load from File] -> select a sorted.bam

![igv2.png](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/igv2.png)

* for example: AT1G67090.1

![igv3.png](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/igv3.png)

Download (if IGV is not installed)
* https://software.broadinstitute.org/software/igv/download

![IGV_down.png](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/IGV_down.png)

Note
* The '+' and '-' button at the top right makes the view zoom-in and -out
* Drag&drop in a region of the top space makes also zoom-in
* Drag&drop in the alignment region makes the view moving to left and right
* You can also use *samtools tview* to see the alignment result on command line

Check mapping quality
```sh
 (Log-in the server)
 (go to your working directory)
 $ source /usr/local/ngseq/etc/lmod_profile
 $ module load QC/SAMStat/1.5.1
 $ module load Tools/samtools/1.10
 $ samstat star_out_w2271_1_hal/w2271_1_hal.sorted.bam
 (Log out from the server)
 (Download the sorted.bam.samstat.html)
 (Look at it using web-browser)
```

![samstat.png](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/samstat.png)

**Question 1**
* What kind of files are generated after STAR mapping?

## Let's Try!!

Mapping

1. Map the following reads to A. thaliana genome using STAR
	* w2271_1_lyr.fastq
	* w2271_2_hal.fastq
	* w2271_2_lyr.fastq
	* w2271_3_hal.fastq
	* w2271_3_lyr.fastq
2. Check the mapping result using IGV

Cheating
```sh
  $ cp /scratch/bio373_2020/data/star_batch.sh .
  $ bash star_batch.sh
```

star_batch.sh
```bash
$ cat star_batsh.sh
source /usr/local/ngseq/etc/lmod_profile
module load Aligner/STAR/2.7.4a
mkdir -p star_out_w2271_1_lyr
mkdir -p star_out_w2271_2_hal
mkdir -p star_out_w2271_2_lyr
mkdir -p star_out_w2271_3_hal
mkdir -p star_out_w2271_3_lyr
STAR --genomeDir ./athal_genome --runThreadN 2 --genomeLoad NoSharedMemory --outSAMstrandField intronMotif --outFileNamePrefix star_out_w2271_1_lyr/ --readFilesIn ./akam_samples2/w2271_1_lyr.fastq
STAR --genomeDir ./athal_genome --runThreadN 2 --genomeLoad NoSharedMemory --outSAMstrandField intronMotif --outFileNamePrefix star_out_w2271_2_hal/ --readFilesIn ./akam_samples2/w2271_2_hal.fastq
STAR --genomeDir ./athal_genome --runThreadN 2 --genomeLoad NoSharedMemory --outSAMstrandField intronMotif --outFileNamePrefix star_out_w2271_2_lyr/ --readFilesIn ./akam_samples2/w2271_2_lyr.fastq
STAR --genomeDir ./athal_genome --runThreadN 2 --genomeLoad NoSharedMemory --outSAMstrandField intronMotif --outFileNamePrefix star_out_w2271_3_hal/ --readFilesIn ./akam_samples2/w2271_3_hal.fastq
STAR --genomeDir ./athal_genome --runThreadN 2 --genomeLoad NoSharedMemory --outSAMstrandField intronMotif --outFileNamePrefix star_out_w2271_3_lyr/ --readFilesIn ./akam_samples2/w2271_3_lyr.fastq

module load Tools/samtools/1.10
samtools sort -O BAM -o star_out_w2271_1_lyr/w2271_1_lyr.sorted.bam star_out_w2271_1_lyr/Aligned.out.sam
samtools sort -O BAM -o star_out_w2271_2_hal/w2271_2_hal.sorted.bam star_out_w2271_2_hal/Aligned.out.sam
samtools sort -O BAM -o star_out_w2271_2_lyr/w2271_2_lyr.sorted.bam star_out_w2271_2_lyr/Aligned.out.sam
samtools sort -O BAM -o star_out_w2271_3_hal/w2271_3_hal.sorted.bam star_out_w2271_3_hal/Aligned.out.sam
samtools sort -O BAM -o star_out_w2271_3_lyr/w2271_3_lyr.sorted.bam star_out_w2271_3_lyr/Aligned.out.sam
```

Note
* It will take about 30 minutes.
* When you run several commands at once, you can use *shell script*.
* Now you can also run *samstat* for all the BAM files by making a shell script.

**Question 2**
* How many reads are mapped with >= MAPQ 30 in each sample?
	* w2271_1_lyr.fastq
	* w2271_2_hal.fastq
	* w2271_2_lyr.fastq
	* w2271_3_hal.fastq
	* w2271_3_lyr.fastq
