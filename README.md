# NEMs ETL Pipeline

This repo is the start of a refactored NEMs pipeline. The primary purpose of this pipeline is to convert raw sequencing data (fastq files) into NEMs. In order, to do this the pipeline is going to be able to conduct alignment, feature count, differential expression, transformation of differential expression to S-Gene Matric, and NEMs network reconstruction. 

## Pipeline Steps:
	1. Alignment (needs to be implmented)
	2. Feature Counting 
	3. Differential Expression Analysis (DEA)
	4. Transform DEA Results to S-Gene Matrix
	5. Network Reconstruction with NEMs

## TODO List:
- Assess Reconstructed Networks stability and quality
- Enable Differential Topological Analysis of Networks from multiple conditions
- Implement Alignment step of pipeline

