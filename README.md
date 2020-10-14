# ![CNV2SV](images/cnv2sv.png) 

## Goal
  
The root cause of the copy number variations is from the underlying genome structure changes. Most detection methods / technologies only take care about whether there are additional copies or missing at one particular locus.  In the light of better DNA detection / sequencing technology, it is possible to construct a full picture of a genome. We would like to address the missing links from CNV detection to the full genome sequence information. In the future, this may help to understand whether the detailed information, e.g. breaking points, of the additional or missing copies are indeed important as makers for pathogenic effects.

## Usage

Required input for CNV2SV linking:

* CNV calls in BED or VCF formats (from Parliament2, Control-FREEC)
* Putative SVs from genome-genome alignment (vcf format, from dipcall)
* .fa files for both assembly versions

Required Python (3.8.\*) packages for CNV2SV linking:

* intervaltree
* mappy
* pyfaidx

Additionally required for visualization of the linking results:

* R, circlize package
* matplotlib
* seaborn
* pandas

Additional dependencies:

* Control-FREEC (provides CNV calls from short read data, requires per chromosome reference FASTA)

Recommended usage (via DNAnexus):

* Read alignments in BAM format for processing with Parliament2
* Genome assembly in FASTA format for processing with dipcall

Quickstart using example data (CHM13 vs GRCh38):

```bash
#Download reference
curl -O -J https://s3.amazonaws.com/nanopore-human-wgs/chm13/assemblies/chm13.draft_v1.0.fasta.gz
curl -O -J https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip chm13.draft_v1.0.fasta.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

#Download CNV calls from T2T CH13 short reads aligned against GRCh38
curl -O -J https://github.com/collaborativebioinformatics/CNV2SV/blob/main/cnvlink/input_cnv_calls/parliament2/PCRfree.cnvnator.vcf

#Download SV calls from genome-genome alignment of T2T CHM13 vs GRCh38
curl -O -J https://github.com/collaborativebioinformatics/CNV2SV/blob/main/cnvlink/chm13_alignment/chm13_grch38.pair.vcf

#Run cnvlink
#Links CNV calls from PCRfree.cnvnator.vcf to SVs in chm13_grch38.pair.vcf (assumes all calls were made against GRCh38.no_alt_analysis_set.fa as ref).
#Output is saved to cnvlink_out_parliament_cnvnator
python cnvlink.py PCRfree.cnvnator.vcf chm13_grch38.pair.vcf GRCh38.no_alt_analysis_set.fa cnvlink_out_parliament_cnvnator


#Visualize CNV-SV linkage statistics
python visualize.py cnvlink_out_parliament_cnvnator

#Circlize plot connecting position of CNV calls to linked SV insertion sites
python cnv_tsv_to_bed.py cnvlink_out_parliament_cnvnator/cnvlink_out.tsv
R circ.R

```

## Output format

cnvlink.py produces linkage information for each CNV as .tsv files. Each row links a CNV to its best matching SV. Information about the additional links with lower scores is included as well.

| Field index | Field name | Description |
| --------------- | --------------- | --------------- |
| 0| cnv_id | The id for the CNV, assigned by the CNV caller |
| 1| status | Summarizes the SV linkage status of this CNV (no_match, match_adj, perfect_adj, match_far, perfect_far) |
| 2| reason | Reason why the CNV could not be linked to an adjacent SV, if any. |
| 3| svtype | DUP or DEL, dependeing on whether the CNV was called as copy number gain or loss |
| 4| good_matches_adj | Total number of valid adjacent SV links discovered |
| 5| good_matches_far | Total number of valid distant SV links discovered |
| 6| alignment_matches | Number of matching characters in the CNV-SV event alignment (mlen value from mappy) |
| 7| cnv_chr | Chromosome the CNV was called on |
| 8| cnv_start | Start position of the CNV |
| 9| cnv_length | Length of the CNV |
| 10| asm_chr | Chromosome of the best linked SV |
| 11| asm_start | Start position of the best linked SV |
| 12| asm_length | Length of the best linked SV |
| 13| all_good_matches_adj | List of all valid adjacent SV links, separated by commas, in the format chr:start-length |
| 14| all_good_matches_far | List of all valid distant SV links, separated by commas, in the format chr:start-length |
| 15| aln_cigar | CIGAR string of the alignment of the CNV sequence to the one of the best linked SV |
| 16| aln_length | Length of the alignment of the CNV sequence to the one of the best linked SV |
| 17| aln_NM | Number of mismatches in the alignment of the CNV sequence to the one of the best linked SV |
| 18| aln_mapq | Mapping quality of the alignment of the CNV sequence to the one of the best linked SV |
| 19| cnv_seq | Full sequence of the called CNV |
| 20| asm_seq | Full sequence of the best linked SV |


## Limitations
* Currently, analysis only includes duplication events

## Aims

We will focus on a couple of genomes that have data from various sequencing technologies where the best genome assemblies and short reads are easy to get. The first target is the T2T2 CHM13. While it is a haploid genome and it will not reflect the "real world" use cases, it will simplify the "equation" so we can develop methods before jumping into a jungle where the information is too complicated to analyze initially.  

The general idea is straight-forward. We will detect CNV using short read data. The alignments and SV (not CNV) calls are straight forward, there are off-shelf solutions. (Yih-chii, Chai and I from DNAnexus has done some preliminary work in the last couple days, so we have some jump-start first without dealing with those long alignment computation tasks.) 

### Aim1

figure out the right way to call CNV from a germline only WGS data. It seems to me most CNV code in the market is using normal-tumor pairs for calling CNV or only focusing on exomes. We need to find useful tools for germline WGS or come up with some quick way to reanalyze the SVs or variant calls (e.g. het-variant call from CHM13) to identify CNV candidates.

### Aim2 

With the CNV candidates called on GRCh38, we will need to find a mapping / liftover, or quick ad-hoc alignment using existing tools to get the sequences from a CHM13 assembly.

### Aim3 

Make some visualization of the identified CNV, potentially focus on pathogenic alleles.

### Aim4

Explore the infrastructure for automation of these processes and make a gallery all related to CNV / SV for demonstration.

## Methods

<!--- ## Awesome Logo -->

Overview Diagram

![Data Processing](/images/CNV2SV.svg)

<!---
# Software Workflow Diagram
-->

### CNV/SV calling from short read data

We are using Parliament for SV/CNV calls from the short read data. We rely on both CNVnator individual calls and combined calls to get the potential duplication locations. Additionally we also run Control-FREEC to estimate the copy number for different regions of the genome based on the short read data. The Control-FREEC output is converted into BED and VCF files for downstream processing.

### Assembly alignment and SV calling

We use the draft assembly of CHM13 T2T genome. We align it to the reference GRCh38 human genome using dipcall, by building a simulated diploid genome via using two copies of CHM13. We then construct dot plots of the alignment to identify potential regions of interest, and extarct a VCF for downstream processing.

### CNV/SV linking

Moritz Smolka has developed a Python script for merging short read and assembly based CNV/SV calls to locate regions in which the duplication events overlap. We are currently working on additional visualization scripts for the final product.

![cnvlink.py](/images/cnvlink_workflow.png)

## Results

### CNVlink raw output

Here, an example for the raw output (CNV - SNV links) produced by cnvlink.py is shown. Upcoming work includes visualizing key elements from these raw data.
![First results](/images/data_example_parliament_new.png)

### CNV/SV linkage information plot

The figures below shows linked CNV calls (from CNVnator) and respective SV insertion calls. The site of the CNV call corresponds to the wider end of a chord, and the insertion SV corresponds to the narrower end. First figure shows the best match used for linking, while the second figure shows all potential SV matches that have the alignment identiity of at least 80%.
![Best CNV2SV link](/images/cnvlink_out_parliament_cnvnator_best_match_resized.pdf)
![All CNV2SV links](/images/cnvlink_out_parliament_cnvnator_all_matches_resized.pdf)
