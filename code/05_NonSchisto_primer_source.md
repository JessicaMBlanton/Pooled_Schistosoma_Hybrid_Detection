# Primer set contribution to non-Schistosoma reads

Jessica Blanton

Last updated 2026-10-08

---

### **Purpose**

There is a lot of human sequence in the urine field samples.  Used primer sequences to find out which PCR rxns contribute most.
These libraries were prepared using the rapid barcode kit chemistry with random fragmentation.  Do not expect primer sequences at all ends, so this analysis is an approximation.
Evaluate which primer sets (18S, ITS, and COX1 genes) contribute the most unmapped reads.

The remaining` unmapped_primercat.txt` and `mapped_primercat.txt` files are input for calculations and plotting in R.

### **Steps**

1. Separate read data into mapped and unmapped (to target regions of *Schistosoma* species)
1. Search for primer sequences in all data
1. Get names for read linked to each primer set
1. Get names all reads mapped/unmapped to Sh refs
1. Tally all mapped/unmapped reads attributed to each primer set
 
### **Programs**

- SeqKit v2.10.0
- Minimap2 v2.30-r1287
- Cutadapt v5.2

### **Inputs**

Sample list

* `~/SHyb_2025/sample_IDS.txt `

Reference marker sequences

* `Sh_hyb_refs.fasta`

Trimmed reads

* `~/SHyb_2025/trim_combo/combo_*_trim.fastq.gz`


Primer sequences to look for

* OR-18S: GCACCAGACTTGCCCTCCAATTGGTCC
* OF-18S: GCATTTATTAGAACAGAACCAAYCGGGCG
* OR-ITS: TCGTGCGTATTACACACACCATCGGTACAAACC
* OF-ITS: GCATGCAAATCCGCCCCGTTATTGTTCCT
* Asmit.1: TTTTTTGGTCATCCTGAGGTGTAT
* Asmit.2: TAAAGAAAGAACATAATGAAAATG

### **Outputs**

Mapping read counts (also available from [this repository](../data/))


- `unmapped_primercat.txt`
- `mapped_primercat.txt`


Logfiles per sample

- cutadapt_*.log

## Commands

Nanopore data: use cutadapt options:

`-b ADAPTER, --anywhere ADAPTER
                      Sequence of an adapter that may be ligated to the 5' or 3' end...`
                        
`--match-read-wildcards Interpret IUPAC wildcards in reads`

```bash
mkdir ~/SHyb_2025/unmap_primer/
cd ~/SHyb_2025/unmap_primer/

# Set up files to record results from looping
echo 'sample;chr;unmapped_reads' | tr ";" "\t" > unmapped_primercat.txt
echo 'sample;chr;mapped_reads' | tr ";" "\t" > mapped_primercat.txt

# Loop over all samples with minimap followed cutadapt per primer set

for i in `cat ~/SHyb_2025/sample_IDS.txt` ; do

	# mapping
	
	minimap2 -ax map-ont -t 10 -k10 -w5 -sr --secondary=no -O 8,24 -E 4,2 \
	~/SHyb_2025/Sh_hyb_refs.fasta \
	~/SHyb_2025/trim_combo/combo_${i}_trim.fastq.gz |
	samtools view -f 4 --threads 10 | samtools fasta > reads_unmapped.fasta
	
	minimap2 -ax map-ont -t 10 -k10 -w5 -sr --secondary=no -O 8,24 -E 4,2 \
	~/SHyb_2025/Sh_hyb_refs.fasta \
	~/SHyb_2025/trim_combo/combo_${i}_trim.fastq.gz |
	samtools view -F 4 --threads 10 -h | samtools fasta > reads_mapped.fasta

	############################################################
	# 18S
	
	cutadapt -j 10 \
	-b GCACCAGACTTGCCCTCCAATTGGTCC \
	-b GCATTTATTAGAACAGAACCAAYCGGGCG \
	--match-read-wildcards --error-rate 0.2 --overlap 18 \
	--report=minimal --action=none --fasta -o trimmed.fa --untrimmed-output noprim.fa \
	~/SHyb_2025/trim_combo/combo_${i}_trim.fastq.gz | column -t 1>> cutadapt_${i}.log 
	
	seqkit seq trimmed.fa --name --only-id | sort > 18S_reads_primerset.txt
	
	seqkit grep -f 18S_reads_primerset.txt reads_mapped.fasta |
	seqkit seq --name --only-id | wc -l  | 
	sed 's/^/set_18S\t/' | sed "s/^/${i}\t/" >> mapped_primercat.txt
	
	seqkit grep -f 18S_reads_primerset.txt reads_unmapped.fasta |
	seqkit seq --name --only-id | wc -l  | 
	sed 's/^/set_18S\t/' | sed "s/^/${i}\t/" >> unmapped_primercat.txt
	
	############################################################	# ITS
	
	cutadapt -j 10 \
	-b TCGTGCGTATTACACACACCATCGGTACAAACC \
	-b GCATGCAAATCCGCCCCGTTATTGTTCCT \
	--match-read-wildcards --error-rate 0.2 --overlap 18 \
	--report=minimal --action=none --fasta -o trimmed.fa --untrimmed-output noprim.fa \
	~/SHyb_2025/trim_combo/combo_${i}_trim.fastq.gz | column -t 1>> cutadapt_${i}.log 
	
	seqkit seq trimmed.fa --name --only-id | sort > ITS_reads_primerset.txt
	
	seqkit grep -f ITS_reads_primerset.txt reads_mapped.fasta |
	seqkit seq --name --only-id | wc -l  | 
	sed 's/^/set_ITS\t/' | sed "s/^/${i}\t/" >> mapped_primercat.txt
	
	seqkit grep -f ITS_reads_primerset.txt reads_unmapped.fasta |
	seqkit seq --name --only-id | wc -l  | 
	sed 's/^/set_ITS\t/' | sed "s/^/${i}\t/" >> unmapped_primercat.txt
	
	############################################################	# COX1
	
	cutadapt -j 10 \
	-b TTTTTTGGTCATCCTGAGGTGTAT \
	-b TAAAGAAAGAACATAATGAAAATG \
	--match-read-wildcards --error-rate 0.2 --overlap 15 \
	--report=minimal --action=none --fasta -o trimmed.fa --untrimmed-output noprim.fa \
	~/SHyb_2025/trim_combo/combo_${i}_trim.fastq.gz | column -t 1>> cutadapt_${i}.log 
	
	seqkit seq trimmed.fa --name --only-id | sort > COX1_reads_primerset.txt
	
	seqkit grep -f COX1_reads_primerset.txt reads_mapped.fasta |
	seqkit seq --name --only-id | wc -l  | 
	sed 's/^/set_COX1\t/' | sed "s/^/${i}\t/" >> mapped_primercat.txt
	
	seqkit grep -f COX1_reads_primerset.txt reads_unmapped.fasta |
	seqkit seq --name --only-id | wc -l  | 
	sed 's/^/set_COX1\t/' | sed "s/^/${i}\t/" >> unmapped_primercat.txt
	
	###### aggregate readnames to check overlap later: ######
	
	cat 18S_reads_primerset.txt >> all_18S_reads_primerset.txt
	cat ITS_reads_primerset.txt >> all_ITS_reads_primerset.txt
	cat COX1_reads_primerset.txt >> all_COX1_reads_primerset.txt

done
```
Check that reads are unique to primer sets across all samples. No results for each comparison indicates no duplicates

```bash
cat all*reads_primerset.txt | sort | uniq -c | grep "^[ ]*1" -v

```
If there are duplicates (cutadapt settings too permissive), see where. No results for each comparison indicates no duplicates

```bash
echo 'ribo'
cat all_18S_reads_primerset.txt all_ITS_reads_primerset.txt | sort | uniq -c | grep "^[ ]*1" -v

echo '18SvCOX'
cat all_18S_reads_primerset.txt all_COX1_reads_primerset.txt | sort | uniq -c | grep "^[ ]*1" -v

echo 'ITSSvCOX'
cat all_ITS_reads_primerset.txt all_COX1_reads_primerset.txt | sort | uniq -c | grep "^[ ]*1" -v
```
Cleanup intermediate files

```bash
rm reads_unmapped.fasta reads_mapped.fasta trimmed.fa noprim.fa *reads_primerset.txt *mapped.fasta

```















