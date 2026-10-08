# Read Processing Nanopore raw reads

Jessica Blanton

Last updated 2025-10-04

---

## **Purpose**
Process raw data from Oxford Nanopore sequencing of datasets containing pooled 18S, ITS, and COX1 marker amplifications

- Libraries were prepared with Rapid Sequencing V14 kit SQK-RBK114.24 and sequenced on FLO-MN114 flow cells
- Basecalled with dorado's "super sccurate" model `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` via minknow desktop GUI
- Barcodes were assigned per sample, basecalling combined markers into one datafile during 
- This workflow uses looping. Alternatively could easily switch GNU-parallel in case of larger datasets.

-
### **Steps**

1. Trim Barcodes
1. Length filter
2. Record length stats
 
### **Required programs**

- Porechop v0.2.4
- SeqKit v2.10.0

### **Inputs**

Raw reads

* `~/SHyb_2025/reads_combo_100425_SUP05/combo_*_raw.fastq.gz `

Sample Info

* `~/SHyb_2025/run_metadata.txt`

### **Outputs**

Trimmed reads per sample

* `~/SHyb_2025/trim_combo/combo_*_trim.fastq.gz`

Read file stats

* `combo_raw_stats.txt `
* `combo_trim_stats.txt`

-
### Commands

Get list of samples from metadata 

```bash
cd ~/SHyb_2025/
cut -f1 run_metadata.txt > sample_IDS.txt
```
#### Remove barcodes, length filter to > 200 bp

```bash
mkdir ~/SHyb_2025/trim_combo

for i in `cat ~/SHyb_2025/sample_IDS.txt` ; do
	porechop \
	-i ~/SHyb_2025/reads_combo_100425_SUP05/combo_${i}_raw.fastq.gz \
	-o trim_combo/porechop_${i}.fastq.gz \
	--verbosity 1 \
	-t 20 \
	1> porechop_${i}.log
	
	seqkit seq -j 20 --min-len 200 --max-len 1000 trim_combo/porechop_${i}.fastq.gz -o trim_combo/combo_${i}_trim.fastq.gz

	rm trim_combo/porechop_${i}.fastq.gz
done

```
Get read statistics

```bash
cd ~/SHyb_2025/reads_combo_100425_SUP05
seqkit stats -j 20 --tabular combo_*_raw.fastq.gz ../combo_raw_stats.txt

cd ~/SHyb_2025/trim_combo
seqkit stats -j 20 --tabular combo_*_trim.fastq.gz ../combo_trim_stats.txt
```

```
cd ~/SHyb_2025/
wc -l combo_raw_stats.txt combo_trim_stats.txt sample_IDS.txt

```
