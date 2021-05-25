Contents
1. [URL](#1-short-url)
2. [B-Fabric registration](#2-b-fabric-registration)
3. [First task: login the server](#3-first-task-login-the-server)
4. [Check the server specs](#4-check-the-server-specs)
5. [Make your own working directory](#5-make-your-own-working-directory)
6. [Calculate G-C content](#6-calculate-g-c-content)
7. [FastQC](#7-fastqc)
8. [Advanced practice](#8-advanced-practice)
9. [For Unix experienced students](#9-for-unix-experiencd-students)
10. [Appendix: Text Editor](#10-apendix-text-editor)

----
# 1. URL

* https://gist.github.com/masaomi/38be7b693b63a51ed431b3f79be724b1

# 2. B-Fabric registration

Link
* https://fgcz-bfabric.uzh.ch

**FGCZ servers**
![fgcz_servers.png](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/fgcz_servers.png)

* B-Fabric: Project management, user management
   * https://fgcz-bfabric.uzh.ch
* gStore: File server, all results are stored here
   * https://fgcz-gstore.uzh.ch/projects/p3662/
* SUSHI: Application server, typical NGS applications can run
   *  https://fgcz-sushi.uzh.ch

Note
* All projects at FGCZ are managed under B-Fabric system
* You can access the NGS data after the sequencing through B-Fabric, gStore, or SUSHI
* If you have your account at FGCZ, you do not have to register again

ToDo
1. Go to https://fgcz-bfabric.uzh.ch/
2. Click [Register](https://fgcz-bfabric.uzh.ch/bfabric/common/edit-user.html)
3. Follow the instructions and input your information

After the registration,
* Please send your account name to masaomi.hatakeyama@ieu.uzh.ch
* Please wait for a while until Masa sets-up your configuration
* Login B-Fabric again, and please make sure you have a new project (p1535) by clicking [Projects]-[List My]


# 3. First task: login the server

1. Start and open a terminal

![terminal](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/terminal.png)

2. Login the server
```bash
 $ ssh your_BFabric_account_name@172.23.xx.xx
 your_BFabric_account_name@172.23.xx.xx's password:
```

Server list
1. 172.23.30.6 (fgcz-kl-003)
2. 172.23.30.200 (fgcz-kl-004)
3. 172.23.32.84 (fgcz-kl-005)




Note
* The password characters are not shown in the terminal, but the system can recognize your password. Please press enter-key after typing your password
* If you fail in log-in, please try it again. If you fail several times, please confirm your account name and password by logging-in B-Fabric.

Reference
* *ssh*: connect to a remote computer with Secure SHell Protocol
* *scp*: Secure CoPy, file copy via ssh protocol between local computer and a remote server

**What is SSH?**

![ssh](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/ssh.png)

Secure Shell: (Wikipedia) http://en.wikipedia.org/wiki/Secure_Shell
```
 Secure Shell (SSH) is a cryptographic network protocol for secure data communication, 
 remote command-line login, remote command execution, and other secure network s
 ervices between two networked computers that connects, via a secure channel over 
 an insecure network, a server and a client (running SSH server and SSH client programs, respectively)
```

```bash
 $ hostname
   (make sure you are working on the server fgcz-kl-00x, not on Mac)
 $ pwd
```

# 4. Check the server specs

Check the system (OS) information
```bash
 $ uname -a
```

Check CPU and RAM information of the server (How many CPUs are there? and How many bytes are the total RAM capacity?)
```bash
 $ cat /proc/cpuinfo
 $ cat /proc/meminfo
```

Check the current running process
```bash
 $ htop
 (in order to exit from htop, type 'q')
```

Reference
* *pwd*: show the current working directory path
* *ls*: show files and sub directories in the current directory
* *uname*: show OS information
* *cat*: show text file contents
* */proc/cpuinf, meminfo*: they contain computer hard ware information
* *htop*: show running processes, memory usage, and so on.

# 5. Make your own working directory

Please make your directory under /scratch as follows:

```bash
 $ cd /scratch/bio373_2020
 $ mkdir your_name
```

* *cd*: Change Working directory
* *mkdir*: MaKe a new DIRectory

Note
* Please DO NOT type just **your_name** but **YOUR NAME**, such as *masaomi*, *tim*, to identify who uses it  (e.g. first name, family name, nick name, whatever, but do not type just like *your_name*)
* You should use this directory for your command practice during the course, rather than home directory, because no enough disk space is allocated to a user home directory.
* All data should be saved in this directory, otherwise your data might be deleted by another user or you might delete other user's data by accident

Copy a sample file to your working directory

```bash
 $ cp /scratch/bio373_2020/data/TAIR10_chr1.fa.gz /scratch/bio373_2020/your_name/
 $ ls
```

Note
* NOT just *your_name* but **YOUR NAME** again.

Change the current working directory to your directory

```bash
 $ cd /scratch/bio373_2020/your_name
 $ ls
```
* Make sure that there is **TAIR10_chr1.fa.gz file**.

Decompress the compressed file & confirm the decompressed file
```bash
 $ gunzip -c TAIR10_chr1.fa.gz > TAIR10_chr1.fa
 $ ls
```

Note
* *gunzip*: decompress a compressed file by *gzip* command
* *-c*: keep the original file and make an output in standard output
* *>*: (redirection) change the standard output to a file

Check if the file is really decompressed
```bash
 $ ls
 $ file TAIR10_chr1.fa.gz
 $ file TAIR10_chr1.fa
```

Note
* *file*: check the file type

Check and compare the file sizes between compressed and decompressed files with different options
```bash
 $ ls -l
 $ ls -lh
```

Note
* *-l*: show detail information of files
* *-lh*: file size is shown with unit prefix (k, M, G, ...)

Show the first and last some lines of the FASTA file
```bash
 (Check the differences of the commands below)
 $ head TAIR10_chr1.fa
 $ tail TAIR10_chr1.fa
 $ head -n 100 TAIR10_chr1.fa
 $ head -n 100 TAIR10_chr1.fa|less
 $ less TAIR10_chr1.fa
 (To exit less command mode, type 'q', to go to the next page, type space or enter key)
```

Note
* *head*: show first 10 lines
* *tail*: show last 10 lines
* | (pipe): command output is passed to the next command as an input
* *less*: show text file data one screen by screen, typing 'h' key shows sub-commands

List all the gene annotation
```bash
 $ grep '>' TAIR10_chr1.fa
 $ grep '>' TAIR10_chr1.fa | less
 (to exit less condition, type 'q', to go to the next page, type space)
```

Count the total number of lines, words, and characters
```bash
 $ wc TAIR10_chr1.fa
```

Count only the total number of lines
```bash
 $ wc -l TAIR10_chr1.fa
```

Count only the number of characters
```bash
 $ wc -c TAIR10_chr1.fa
```

# 6. Calculate G-C content

Count how many genes are defined in the FASTA file (TAIR10_chr1.fa, **A. thaliana** Chromosome1)
 (Count only the gene annotation lines)

```bash
 $ grep '>' TAIR10_chr1.fa | wc -l
```

Count the total number of nucleotide bases
 (Do not count the gene annotation line)

 Hint: **grep** with **-v** option can skip a specific line. Check **grep --help** or search by Google!

```bash
 $ grep -v '>' TAIR10_chr1.fa | wc -c
 $ grep -v '>' TAIR10_chr1.fa | wc -l
```

**Question1**
* How much is the GC content in the FASTA file (TAIR10_chr1.fa, **A. thaliana** chromosome1)?

Hints
* Hint1: sed command can replace characters in a line
* Hint2: wc can count lines and characters

```bash
 $ echo 'AAATTTGGGCCC' | sed 's/A/Z/g'
 ZZZTTTGGGCCC
```

Reference/Hint
* GC-Content (wikipedia):http://en.wikipedia.org/wiki/GC-content
* *sed* command: replace strings, -e "s/[target character]/[replace character]/g"
* regular expression: character pattern, wikipedia http://en.wikipedia.org/wiki/Regular_expression
* line break is also counted as one character in wc command

Compress TAIR10_chr1.fa with different options (What is the difference?)
```bash
 $ gzip -h
 (for example)
 $ time gzip --fast -c TAIR10_chr1.fa > TAIR10_chr1.fa.fast.gz
 $ time gzip --best -c TAIR10_chr1.fa > TAIR10_chr1.fa.best.gz
 $ ls -lh
```

**Question2**
* How much (How many bytes) is the total data size (after decompression)?

# 7. FastQC

Copy the following file in your working directory (/scratch/bio373_2020/your_name)
```bash
  $ hostname
    (make sure you are working on the server fgcz-kl-00x)
  $ pwd
    (make sure you are working in the working directory, /scratch/bio373_2020/your_name)
  $ cp /scratch/bio373_2020/data/akam_samples.tgz ./
    (do not forget the dot at the end, which means the current working directory)
  $ ls
```

Extract (decompress) it under your directory

```bash
 $ ls
  (make sure there is certainly akam_samples.tgz in your current working directory)
 $ tar zxvf akam_samples.tgz
 $ ls
```

Check fastq file
```bash
 $ ls
  (make sure there is certainly akam_samples/ directory)
 $ cd akam_samples
 $ cat akam_sample1.fastq | more
 (to exit, type 'q')
```

Note
* There are 4 samples data.
* sample1 and 2 are in control (normal condition)
* sample3 and 4 are in insecticide treated
* sample1 and 3 libraries were prepared by a FGCZ staff
* sample2 and 4 libraries were prepared by students

FastQC (Quality Control)
```bash
 $ source /usr/local/ngseq/etc/lmod_profile
 $ module load QC/FastQC/0.11.9
 $ ls
   (make sure there is certainly akam_sample1.fastq)
 $ fastqc akam_sample1.fastq
 $ ls
   (check the generated files)
```

Note
* module: setup environmental variables appropriately for **fastqc** command. This command is configurated and available only in fgcz-kl-00x servers as the example above

**Question3**
* What kind of files are generated by fastqc command?

Logout from fgcz-kl-00x
```bash
 $ exit
```

Download the result file from the server %red% DO NOT forget the dot (it means the current directory)!!
```bash
  $ hostname
     (You are sure it is on Mac, not on fgcz-kl-00x)
  $ cd ~/Desktop
  $ scp your_BFabric_accout_name@172.23.xx.xx:/scratch/bio373_2020/your_name/akam_samples/akam_sample1_fastqc.html ./
```

Note
* Be careful of the small letter **l** (L) and the number **1** (one)

* Open the downloaded html file on a web-browser

![fastqc_report](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/fastqc_report.png)

**Check it out!!**
* What point is good/bad with the QC report?

# Let's Try!!

* Check the other samples

```bash
 akam_sample2.fastq
 akam_sample3.fastq
 akam_sample4.fastq
```

**Question4**
* How different are they?

# 8. Advanced practice

* MultiQC can gather all the fastqc results into one file nicely.
* After making all the FastQC results of 1-4 samples, in the **akam_examples** directory try to run  as follows:
```bash
 $ source /usr/local/ngseq/etc/lmod_profile
 $ module load Dev/Python/3.8.3
 $ multiqc .
```

Note
* Please do not forget the last dot after one space.
* You will get **multiqc_report.html** file at the end and let's download it to your local computer by **scp** command, and open it on a web-browser in your local computer

![multiqc](https://raw.githubusercontent.com/masaomi/Bio373_2020/master/png/multiqc.png)

MultiQC
* http://multiqc.info/
* https://doi.org/10.1093/bioinformatics/btw354


# 9. For Unix experienced students

Try to make a shell script using a text editor and run several **fastqc** commands at once.

fastqc_batch.sh
```bash
 source /usr/local/ngseq/etc/lmod_profile
 module load QC/FastQC/0.11.9
 for i in `seq 1 3`
 do
   fastqc akam_samples/akam_sample${i}.fastq
 done
```

To execute it
```bash
 $ bash fastqc_batch.sh
```

Note
* Please make the shell script using a text editor (it is fine to make it locally and upload to the server, but for the text editor practice, let's try to use one of the CLI text editors)
* In this example, the shell script file name is **fastqc_batch.sh** but it does not matter whatever the file name is.
* *seq*: make a list of numbers
* *for*: iterate a process between **do** and **done** with assigning each element in the variable **i**
* The iterate variable can be referred to with **$** symbol
* To execute the shell script, you can call by **bash**
* In the shell script, backquotations `` returns the command result. For the application of this, the following script has the same function as the shell script above

Additional samples
* /scratch/bio373_2020/data/akam_samples2.tgz

* Copy it in your working directory
* Extract the archived files (please use tar command)
* Execute FastQC (by the following shell script)

fastqc_batch2.sh
```bash
 source /usr/local/ngseq/etc/lmod_profile
 module load QC/FastQC/0.11.9
 for file in `ls -1 akam_samples2/w2271_*.fastq`
 do
   fastqc $file
 done
```


# Appendix: Text Editor

The followings are famous text editors available on this server.

* vi
* emacs
* nano

* I recommend for you to use **nano** if you have no experience to use a text editor in command line

How to use
* make a new text file
```bash
 $ nano text_file_name.txt
```

* save and exit
```bash
 Ctrl + x, y
```

Note
* The symbol **^** means pressing *Control* Key
