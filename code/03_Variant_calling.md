# Variant calling with FreeBayes

Jessica Blanton

Last updated 2026-10-08

---

### **Purpose**
Identify SNVs present at a minimum of 1% of reads with the Bayesian genetic variant detector FreeBayes v1.3.10.  Parameters are chosen to accomodate pooled samples.

The resulting `mapped.txt` `combo_*_bc.tsv` and `freebayes_all_SNVs.txt` files are input for calculations and plotting in R.

### **Steps**

1. Trim Barcodes
1. Length filter
2. Record length stats
 
### **Programs**

- SeqKit v2.10.0
- Minimap2 v2.30-r1287
- samtools v 1.22.1 
- bam-readcount v1.0.1
- `brc-parser.py` (script downloaded from https://github.com/sridhar0605/brc-parser, saved in `~/opt/`)
- FreeBayes v1.3.10
- vt v0.57721


### **Inputs**

Trimmed reads: 

* `~/SHyb_2025/trim_combo/combo_*_trim.fastq.gz`

Reference marker sequences

* `Sh_hyb_refs.fasta`

Sample list

* `~/SHyb_2025/sample_IDS.txt `

Amplicon positions on *S. haematobium* reference sequences

|Gene|Sh amplicon positions|
|----|---------------------|
| 18S | 232-557 |
| ITS | 274-878 |
| COX1 | 719-1113 |

### **Outputs**

Mapping files, mapping count tallies, and mapping stats- per sample

* `~/SHyb_2025/mapping/combo_${i}_Sh.bam` (one file per sample)
* `~/SHyb_2025/mapping/mapped.txt`
* `~/SHyb_2025/basecounts/combo_*_bc.tsv` (one file per sample)

Variant calls

* `freebayes_all_SNVs.txt`

## Commands

### Map trimmed reads to S. haematobium reference sequences

```bash
cd ~/SHyb_2025/

mkdir ~/SHyb_2025/mapping

for i in `cat ~/SHyb_2025/sample_IDS.txt` ; do
	minimap2 \
	-ax map-ont -t 20 \
	-k10 -w5 -sr --secondary=no \
	-O 8,24 \
	-E 4,2 \
	~/SHyb_2025/Sh_hyb_refs.fasta \
	~/SHyb_2025/trim_combo/combo_${i}_trim.fastq.gz |
	samtools view -F 4 --threads 9 -b | samtools sort --threads 9 -O BAM --write-index -o  ~/SHyb_2025/mapping/combo_${i}_Sh.bam -
done
```
Summarize number of reads mapped

```
echo 'sample;GU257398.1_mapped;Z11976.1_mapped;NC_008074.1_mapped' | tr ";" "\t" > mapped.txt

for i in `cat ~/SHyb_2025/sample_IDS.txt` ; do
	samtools idxstats ~/SHyb_2025/mapping/combo_${i}_Sh.bam |
	grep "*" -v | cut -f3 | tr '\n' '\t' | sed 's/\t$/\n/' |
	sed "s/^/combo_${i}\t/" >> mapped.txt
done

```
Check mapping

```
head mapped.txt | column -t 
```

		sample     GU257398.1_mapped  Z11976.1_mapped  NC_008074.1_mapped
		combo_001  923                1222             26343
		combo_002  1996               1959             9979
		combo_003  8443               6934             1902
		combo_004  1322               1594             12929
		combo_005  11860              6645             5158
		combo_006  18066              9697             8378
		combo_007  8837               5016             7511
		combo_008  18841              11600            6183
		combo_009  35747              11408            1743

Get single nucleotide statistics

```bash
mkdir ~/SHyb_2025/basecounts
cd ~/SHyb_2025/basecounts

# Specify region range on reference to count basecoverage
echo 'GU257398_1;274;878
Z11976_1;232;557
NC_008074_1;719;1113' | tr ";" "\t" > sites

for i in `cat ~/SHyb_2025/sample_IDS.txt` ; do

	bam-readcount \
	--reference-fasta ~/SHyb_2025/Sh_hyb_refs.fasta \
	~/SHyb_2025/mapping/combo_${i}_Sh.bam \
	--min-mapping-quality 10 \
	--site-list sites > combo_${i}_bc.tsv

	# Convert to long format output_parsed.csv
	python ~/opt/brc-parser.py combo_${i}_bc.tsv 

done
```
### Variant calling 

Run freebayes followed by normalizing: Parsimony and left alignment decomposition of bi-allelic block substitutions

```bash
mkdir ~/SHyb_2025/freebayes_vcf_01

for i in `cat ~/SHyb_2025 study_IDS.txt` ; do

	cd ~/SHyb_2025/freebayes_vcf_01
	
	freebayes \
	-f ~/SHyb_2025/Sh_hyb_refs.fasta \
	-b mapping/kmkhyb_combo_${i}_Sh.bam \
	--pooled-continuous \
	--haplotype-length 0 \
	--min-alternate-fraction 0.01 \
	--min-alternate-count 20 \
	--min-coverage 200 \
	--min-mapping-quality 40 \
	--limit-coverage 10000 \
	--mismatch-base-quality-threshold 20 \
	--read-indel-limit 20 \
	--no-partial-observations |
	vt decompose -s  - | # Decomposition
	vt normalize -n -r ~/SHyb_2025/Sh_hyb_refs.fasta - | 
	vt decompose_blocksub - -o combo_${i}_Sh.vcf

done

# Check that all samples have completed successfully 
ls *vcf | wc -l

```
Concatenate to one file for analysis and parse for import to R

```bash
grep "#" ~/SHyb_2025/freebayes_vcf_01/combo_*_Sh.vcf  -v |
cut -f2 -d"/" | sed 's/_Sh.vcf:/\t/g'| tr ";" "\t" \
> ~/SHyb_2025/freebayes_all_SNVs.txt
```




























