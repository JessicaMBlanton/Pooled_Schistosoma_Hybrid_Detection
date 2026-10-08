# Determination of reference variants
Jessica Blanton

Last updated 2026-10-08

---

### **Purpose**

Publicly deposited sequences for the 18S, ITS and COX1 genes from *S. bovis*, *S. curassoni*, *S. guineensis*, *S. mansoni*, and *S. japonicum* were used to determine species-indicative variant positions versus the *S. haematobium* sequences.  

These species were selected based on geographic overlap in Senegal and Gabon and prevalence on the African continent in humans and livestock. The human-infecting *S. japonicum* was included to add phylogenetic breadth.

Thus, the species-indicative SNVs identified here are relative to the taxa being compared, and may not be indicative across all known *Schistosoma* species.

### **Steps**

1. Retrieve reference sequences from 7 *Schistosoma* species. 
1. Sequences were aligned and trimmed to start positions of the *S. haematobium* sequence. 
1. Polymorphic sites differing from *S. haematobium* were identified within each amplicon region using R packages.
 
### **Programs**

Run from terminal

- `muscle v 5.3`
- `SeqKit v2.10.0`

Run in R

- `tidyvers v 2.0.0`
- `seqinr v 4.2.36`

### **Inputs**

Data

* `Ref_18S_seqs_seqs.fasta`
* `Ref_ITS_seqs.fasta`
* `Ref_COX1_seqs.fasta`
	
* Amplicon positions on *S. haematobium* reference sequences
	
	|Gene|Sh amplicon positions|
	|----|---------------------|
	| 18S | 232-557 |
	| ITS | 274-878 |
	| COX1 | 719-1113 |


### **Outputs**

* `18S_muscle.fasta`
* `ITS_muscle.fasta`
* `COX1_muscle.fasta`
* `ref_snps.RData`
* `expected_SNVs_extended.txt`

---

NCBI nr sequences were downloaded where available, or were extracted from NCBI genome assemblies (*)

Accessions:

| Species | 18S | ITS | COX1 |
|-------- | --- | --- | ---- |
| *S. haematobium* | Z11976.1 | GU257398.1 | NC_008074.1: 6333-7874* |
| *S. bovis* | AY157238.1 | MT580950.1 | MF919409.1 - short |
| *S. mansoni* | X53047.1 | KX011041.1 | NC_002545.1:924-2456 |
| *S. curassoni* | AY157236.1 | MT580946.1 | OX104147.1: 573-2113* |
| *S. guineensis* | OX103898.1 | Z21717.1 | OX103896.1: 553-2093 |
| *S. japonicum* | AY157226.1 | FJ852557.1 | KU196409.1: 8928-10571* |
| *S. intercalatum* | CALYCO020000237.1: 16010-17998* | CALYCO020000237.1: 17980-18935* | OX103731.1: 574-2102* |
\* positions extraction from genome

## Commands

#### Prepare multi-sequence files

```bash
mkdir ~/SHyb_2025/Ref_gene_seqs/
cd ~/SHyb_2025/Ref_gene_seqs/

```
Collect sequences into multifasta files by gene and save in directory:

* `Ref_18S_seqs_seqs.fasta`
* `Ref_ITS_seqs.fasta`
* `Ref_COX1_seqs.fasta`

#### Align with muscle

```bash
muscle -in Ref_18S_seqs_seqs.fasta -out 18S_muscle.afa 
muscle -align Ref_ITS_seqs.fasta -output ITS_muscle.afa 
muscle -in Ref_COX1_seqs.fasta -out COX1_muscle.afa 

# Unwrap lines
seqkit sort 18S_muscle.afa --line-width 0 -o 18S_muscle_free.fasta
seqkit sort ITS_muscle.afa --line-width 0 -o ITS_muscle_free.fasta
seqkit sort COX1_muscle.afa --line-width 0 -o COX1_muscle_free.fasta

```
#### Manually trim alignments to *S. haematobium* ends, as we will define the reference position for SNVs relative to the *S. haematobium* sequences. Could be done programmatically, but these are just 3 small files.

Name final files as

* `18S_muscle.fasta`
* `ITS_muscle.fasta`
* `COX1_muscle.fasta`

-

### Identify expected SNV positions for each representative species reference

*The following was done in R*

Import alignments

```R
# Load libraries
library(tidyverse)
library(seqinr)

setwd("~/SHyb_2025/Ref_gene_seqs/")

# Setup comparison pairs
index_pairs <- data.frame(ref = c("Sh", "Sh", "Sh", "Sh", "Sh","Sh"),
                          query = c("Sb", "Sc", "Sg", "Si", "Sj", "Sm"))

# Import aligned sequences from FASTA
Ref_18S_seqs <- read.alignment("~/SHyb_2025/Ref_gene_seqs/18S_muscle.fasta", 
  format = "fasta", forceToLower = FALSE)
seqs_18S <- Ref_18S_seqs$seq
names(seqs_18S) <- Ref_18S_seqs$nam

Ref_ITS_seqs <- read.alignment("~/SHyb_2025/Ref_gene_seqs/ITS_muscle.fasta", 
  format = "fasta", forceToLower = FALSE)
seqs_ITS <- Ref_ITS_seqs$seq
names(seqs_ITS) <- Ref_ITS_seqs$nam

Ref_cox1_seqs <- read.alignment("~/SHyb_2025/Ref_gene_seqs/COX1_muscle.fasta", 
  format = "fasta", forceToLower = FALSE)
seqs_cox1 <- Ref_cox1_seqs$seq
names(seqs_cox1) <- Ref_cox1_seqs$nam
```

### Extract differences from Sh ref for each query sequences, tag amplicon regions

18S

```R
# Initialize new df
ref_snp_18S <- data.frame()

# Find positions where aligned sequences differ, excluding gaps
# Calculate the ungapped coordinates.
for (j in seq_len(nrow(index_pairs))) {
  REF   <- index_pairs[j, 1]
  QUERY <- index_pairs[j, 2]
  
  s1 <- strsplit(unlist(seqs_18S[grepl(paste0(REF, "_"), names(seqs_18S))]), "")[[1]]
  s2 <- strsplit(unlist(seqs_18S[grepl(paste0(QUERY, "_"), names(seqs_18S))]), "")[[1]]
  
  pos1 <- cumsum(s1 != "-")
  pos2 <- cumsum(s2 != "-")
  
  snps <- s1 != s2 & s1 != "-" & s2 != "-"
  
  ref_snp_18S <- rbind(
    ref_snp_18S,
    data.frame(
      Alignment_Pos = which(snps),
      SeqRef_Char   = s1[snps],
      SeqQuery_Char = s2[snps],
      Pos_Ref       = pos1[snps],
      Pos_Query     = pos2[snps],
      Ref           = REF,
      Query         = QUERY
    )
  )
}

# Add label specifying Sh chromosome accession 
ref_snp_18S$chr <- "Z11976_1"

# label positions that fall within amplicon region
ref_snp_18S[ref_snp_18S$Pos_Ref %in% 232:557,"target_desig"] <- "PCR_region"

# Add unique index term for each SNP 
ref_snp_18S <- 
  ref_snp_18S %>%
  unite("index", c("chr", "Pos_Ref", "SeqQuery_Char"), remove = F) 
  
```
ITS - repeat as above

```R
# Initialize new DF
ref_snp_ITS <- data.frame()

for (j in seq_len(nrow(index_pairs))) {
  REF   <- index_pairs[j, 1]
  QUERY <- index_pairs[j, 2]
  
  s1 <- strsplit(unlist(seqs_ITS[grepl(paste0(REF, "_"), names(seqs_ITS))]), "")[[1]]
  s2 <- strsplit(unlist(seqs_ITS[grepl(paste0(QUERY, "_"), names(seqs_ITS))]), "")[[1]]
  
  pos1 <- cumsum(s1 != "-")
  pos2 <- cumsum(s2 != "-")
  
  snps <- s1 != s2 & s1 != "-" & s2 != "-"
  
  ref_snp_ITS <- rbind(
    ref_snp_ITS,
    data.frame(
      Alignment_Pos = which(snps),
      SeqRef_Char   = s1[snps],
      SeqQuery_Char = s2[snps],
      Pos_Ref       = pos1[snps],
      Pos_Query     = pos2[snps],
      Ref           = REF,
      Query         = QUERY
    )
  )
}

ref_snp_ITS$chr <- "GU257398_1"
ref_snp_ITS[ref_snp_ITS$Pos_Ref %in% 274:878, "target_desig" ] <- "PCR_region"

ref_snp_ITS <- 
  ref_snp_ITS %>%
  unite("index", c("chr", "Pos_Ref", "SeqQuery_Char"), remove = F) 
  
```

COX1 - repeat as above

```R
Initialize new DF
ref_snp_cox1 <- data.frame()

for (j in seq_len(nrow(index_pairs))) {
  REF   <- index_pairs[j, 1]
  QUERY <- index_pairs[j, 2]
  
  s1 <- strsplit(unlist(seqs_cox1[grepl(paste0(REF, "_"), names(seqs_cox1))]), "")[[1]]
  s2 <- strsplit(unlist(seqs_cox1[grepl(paste0(QUERY, "_"), names(seqs_cox1))]), "")[[1]]
  
  pos1 <- cumsum(s1 != "-")
  pos2 <- cumsum(s2 != "-")
  
  snps <- s1 != s2 & s1 != "-" & s2 != "-"
  
  ref_snp_cox1 <- rbind(
    ref_snp_cox1,
    data.frame(
      Alignment_Pos = which(snps),
      SeqRef_Char   = s1[snps],
      SeqQuery_Char = s2[snps],
      Pos_Ref       = pos1[snps],
      Pos_Query     = pos2[snps],
      Ref           = REF,
      Query         = QUERY
    )
  )
}

ref_snp_cox1$chr <- "NC_008074_1"
ref_snp_cox1[ref_snp_cox1$Pos_Ref %in% 719:1113,"target_desig"] <- "PCR_region"

ref_snp_cox1 <- 
  ref_snp_cox1%>%
  unite("index", c("chr", "Pos_Ref", "SeqQuery_Char"), remove = F) 
  
```
#### Create master list of reference SNVs

- NOTE: One SNV was discovered between two S. haematobium strains in the course of this study- removed as this is not a species-indicative position.

```R
# Combine SNVs for all markers
ref_snps <- bind_rows(ref_snp_18S,
                      ref_snp_ITS,
                      ref_snp_cox1)

# Remove real variant from in S. haematobium reference Egypt strain vs Mali Ref:
ref_snps <- ref_snps[!grepl("NC_008074_1_786", ref_snps$index),]
```

Examine Reference SNV df 

```R
head(ref_snps)

	  Alignment_Pos SeqRef_Char          index SeqQuery_Char Pos_Ref Pos_Query Ref Query      chr target_desig
	1            87           C  Z11976_1_87_T             T      87         1  Sh    Sb Z11976_1         <NA>
	2            93           T  Z11976_1_93_A             A      93         7  Sh    Sb Z11976_1         <NA>
	3           225           T Z11976_1_225_C             C     225       139  Sh    Sb Z11976_1         <NA>
	4           250           C Z11976_1_250_T             T     250       164  Sh    Sb Z11976_1   PCR_region
	5           297           T Z11976_1_297_C             C     297       211  Sh    Sb Z11976_1   PCR_region
	6           687           T Z11976_1_685_C             C     685       599  Sh    Sb Z11976_1         <NA>
```

1532 SNVs found across all positions vs S. haematobium references

```R
dim(ref_snps)
[1] 1532   10
```

454 SNVs found within PCR amplicons in this study 

```R
dim(ref_snps[ref_snps$target_desig %in% "PCR_region",])
[1] 454  10
```

```R
# Save object for downstream analyses
save(ref_snps, file = "R_analysis/Rdata/ref_snps.RData")

# Export dataframe for reference
ref_snps%>%.[.$target_desig %in% "PCR_region",] %>%
  group_by(chr, Ref, Query) %>%
  summarise(count=n()) %>%
  pivot_wider(., names_from = chr, values_from = count, values_fill = 0) %>%
  write.table(file = "~/SHyb_2025/Ref_gene_seqs/expected_SNVs_pcr.txt", 
              quote = F, sep = "\t", row.names = F)

```



















