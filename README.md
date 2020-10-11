# CNV2SV

## Goal
  
The root cause of the copy number variations is from the underlying genome structure changes. Most detection methods / technologies only take care about whether there are additional copies or missing at one particular locus.  In the light of better DNA detection / sequencing technology, it is possible to construct a full picture of a genome. We would like to address the missing links from CNV detection to the full genome sequence information. In the future, this may help to understand whether the detailed information, e.g. breaking points, of the additional or missing copies are indeed important as makers for pathogenic effects.

## Aims and Methods

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



<!--- ## Awesome Logo -->

Overview Diagram

![Data Processing](/images/CNV2SV.svg)

<!---
# Software Workflow Diagram
-->
