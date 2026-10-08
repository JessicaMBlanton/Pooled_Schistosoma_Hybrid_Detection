# Mapping non-Schistosome/off-target amplification

Jessica Blanton

Last updated 2026-10-08

---

### **Purpose**

Many samples show high proportions of reads not mapping to the *Schistosoma* spp. 18S, ITS, and COX1 gene regions. To evaluate potential of-target amplification host species in urine and laboratory samples, reads were mapped to a library made up of the target regions for the 18S, ITS, and COX of reference *Schistosoma* spp., as well as human, mouse, Syrian hamster, and cattle genomes (representing the field or lab definitive hosts). Mapping parameters are set to map reads competitvely.

The resulting `genome_idxstats.txt` file is input for calculations and plotting in R.

### **Steps**

1. Make single file of all reference sequences and genomes
2. Competitively map reads
3. Get mapping counts
4. Record contig names per genome to associate hits with hosts

### **Programs**

- SeqKit v2.10.0
- Minimap2 v2.30-r1287
- samtools v 1.22.1 

### **Input**

Trimmed reads

* `~/SHyb_2025/trim_combo/combo_*_trim.fastq.gz`

Sample list

* `~/SHyb_2025/sample_IDS.txt `

18S, ITS, and COX1 amplification regions of *Schistosoma* reference sequences
	
* `Ref_18S_seqs_TARMS.fasta`
* `Ref_ITS_seqs_TARMS.fasta`
* `Ref_COX1_seqs_asmit.fasta`

Host genomes, which include MT genome

* Human (*Homo saipiens*, GCF_000001405.26)
* Mouse (*mus musculus*, GCF_000001635.27)
* Syrian hamster (*Mesocricetus auratus*, GCF_017639785.1)
* Domesic Cow (*Bos Taurus*, 	GCF_002263795.3)

### **Output**

Mapping read counts

* `genome_idxstats.txt`

Contig names for whole-genomes

* `Hs_GRCh38_ctg_nms.txt`
* `Mmus_ctg_nms.txt`
* `Maur_ctg_nms.txt`
* `Btaur_ctg_nms.txt`

## Commands

Minimap2 parameters were set to recruit reads competitively between these references (-ax map-ont -t 10 --secondary=no). Cutadapt v5.2 was used to identify read proportions containing primer sets from each of the three marker genes.

Aggregate *Schistosoma* and host reference sequences

```bash
mkdir -p ~/SHyb_2025/map_GRCh38_Mmus_Maur_Srefs/bam_genome
cd ~/SHyb_2025/map_GRCh38_Mmus_Maur_Srefs/

# Examine sequences of marker gene regions for all Schistosoma species

cat Ref_18S_seqs_TARMS.fasta Ref_ITS_seqs_TARMS.fasta Ref_COX1_seqs_asmit.fasta > Ref_seqs_cat.fasta

grep ">" Ref_seqs_cat.fasta 
	>Sb_18S_AY157238_1
	>Sc_18S_AY157236_1
	>Sg_18S_OX103898_1
	>Sh_18S_Z11976_1
	>Sj_18S_AY157226_1
	>Sm_18S_X53047_1
	>Sb_ITS_MT580950_1
	>Sc_ITS_MT580946_1
	>Sg_ITS_Z21717_1
	>Sh_ITS_GU257398_1
	>Sj_ITS_FJ852437_1
	>Sm_ITS_KX011041_1
	>Sb_cox1_MF919409_1
	>Sc_cox1_AY157210_1
	>Sg_cox1_AJ519523_1
	>Sh_cox1_NC_008074_1
	>Sj_cox1_EU340360_1
	>Sm_cox1_MF919424_1

#### Combine all references to map reads between

gzcat GCF_000001405.26_GRCh38_genomic.fna.gz \
GCF_000001635.27_GRCm39_genomic.fna.gz \
GCF_017639785.1_BCM_Maur_2.0_genomic.fna.gz |
GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz \
cat - Ref_seqs_cat.fasta \
> GRCh38_Mmus_Maur_Btaur_ShRefs.fna.gz
```
Check sequence stats of concatenated genomes file

```
seqkit stats GRCh38_Mmus_Maur_Btaur_ShRefs.fna.gz

file                                  format  type  num_seqs  sum_len         min_len  avg_len      max_len
GRCh38_Mmus_Maur_Btaur_ShRefs.fna.gz  FASTA   DNA   8,839     11,165,281,644  205      1,263,183.8  248,956,422
```
Mapping

```bash 
cd ~/SHyb_2025/map_GRCh38_Mmus_Maur_Srefs/

####
# Index reference file
minimap2 -t 10 GRCh38_Mmus_Maur_ShRefs.fna -d GRCh38_rodent_Srefs

# Map reads by sample to combined reference sequences
for i in `cat ~/SHyb_2025/study_IDS.txt` ; do
	minimap2 \
	-ax map-ont -t 10 --secondary=no --split-prefix temp_name \
	GRCh38_rodent_Srefs \
	~/SHyb_2025/trim_combo/combo_${i}_trim.fastq.gz |
	samtools view -F 4 --threads 10 -b | 
	samtools sort -O BAM --write-index -o bam_genome/comp_${i}_GRCh38_rodent_Srefs.bam -
done
```
Extract mapping stats

```bash
for i in `cat ~/SHyb_2025/study_IDS.txt` ; do
samtools idxstats bam_genome/comp_${i}_GRCh38_rodent_Srefs.bam | 
grep "*" -v | sed "s/^/sample_${i};/" | 
tr ";" "\t" > bam_genome/idxstats_${i}.txt
done

echo "sample;chr;chr_length;mapped;unmapped" | tr ";" "\t" > genome_idxstats.txt
cat bam_genome/idxstats_*.txt >> genome_idxstats.txt

# Cleanup- remove large genome file and index
rm GRCh38_rodent_Srefs GRCh38_Mmus_Maur_ShRefs.fna

```
The file `genome_idxstats.txt` contains genome sequences (contigs & chromosomes) spread across multiple fasta entries .  

Record contig names per host genome for association of host source with mapping counts in downstream analysis.

```bash
zgrep ">"  GCF_000001405.26_GRCh38_genomic.fna.gz | tr -d ">" | sed 's/ /\t/' > Hs_GRCh38_ctg_nms.txt

zgrep ">" GCF_000001635.27_GRCm39_genomic.fna.gz | tr -d ">" | sed 's/ /\t/' > Mmus_ctg_nms.txt

zgrep ">" GCF_017639785.1_BCM_Maur_2.0_genomic.fna.gz | tr -d ">" | sed 's/ /\t/' > Maur_ctg_nms.txt

zgrep ">" GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz | tr -d ">" | sed 's/ /\t/' > Btaur_ctg_nms.txt

```









