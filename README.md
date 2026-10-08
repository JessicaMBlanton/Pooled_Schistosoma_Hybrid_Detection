# Pooled_Schistosoma_Hybrid_Detection

This repository contains workflows documenting programs and commands for:

- [Reference species-indicative SNVs](code/01_Determination_of_reference_variants.md)  
- [Processing of sequence reads](code/02_Read_processing.md)
- [Variant calling with FreeBayes](code/03_Variant_calling.md)
- [Identification of off-target host amplification](code/04_Non_Schistosoma_read_mapping.md)
- [Estimate primer set contributions to off-target amplification](code/05_NonSchisto_primer_source.md)


Data files:

01

* Ref_18S_seqs.fasta
* Ref_ITS_seqs.fasta
* Ref_COX1_seqs.fasta
* 18S_muscle.fasta
* ITS_muscle.fasta
* COX1_muscle.fasta
* ref_snps.RData
* expected_SNVs_extended.txt

03

* freebayes_all_SNVs.txt

04

* genome_idxstats.txt
* Btaur_ctg_nms.txt
* Hs_GRCh38_ctg_nms.txt
* Maur_ctg_nms.txt
* Mmus_ctg_nms.txt

05

* mapped_primercat.txt
* unmapped_primercat.txt



To do
These are available on NCBI, but there were some extractions from genomes...

- provide output files in gen.  This documentation could focus on the first steps to get away from any sequence handling (raw or ref)
	
- bc files for frequency counts from mapping
- Metadata - sample_IDS.txt and sample_metadata.txt

Analyses (time permitting)
- _Walk through figure by table to identify other code that should be posted_
