**Question1**
* How much is the GC content in the FASTA file (TAIR10_chr1.fa, **A. thaliana** chromosome1)?

Hints
* Hint1: sed command can replace characters in a line
* Hint2: wc can count lines and characters

```bash
 $ echo 'AAATTTGGGCCC' | sed 's/A/Z/g'
 ZZZTTTGGGCCC
```

Example answer
```bash
$ grep -v '>' TAIR10_chr1.fa | wc -c
25441459
$ grep -v '>' TAIR10_chr1.fa | wc -l
257293
$ grep -v '>' TAIR10_chr1.fa| sed 's/G//g'|sed 's/C//g' | wc -c
15586736
$ echo "scale=3;(25441459-15586736)/(25441459-257293)" | bc
```

**Question2**
* How much (How many bytes) is the total data size (after decompression)?

Answer
```bash
 $ ls -l TAIR10_chr1.fa
 $ ls -l TAIR10_chr1.fa.gz
```
