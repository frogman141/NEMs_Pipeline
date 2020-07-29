#!/home/baker02/miniconda3/envs/Renv/bin/Rscript

# execute the analysis pipeline for the ER data
# each step will write intermediate output files to disk

#! Give your job a name
#SBATCH -J 003_diffexp.NEMpipeline
#! How many cores per task?
#SBATCH --cpus-per-task=6
#! How much memory do you need?
#SBATCH --mem=12G
#! How much wallclock time will be required?
#SBATCH --time=01-00:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL
#! Specify your email address here otherwise you won't recieve emails!
#SBATCH --mail-user=alexander.baker@cruk.cam.ac.uk
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
#SBATCH --no-requeue
#! General partition
#SBATCH -p general
library(dplyr)
library(ggplot2)
library(pheatmap)
library(BiocParallel)

multicore_param = MulticoreParam(6)
register(multicore_param, default=TRUE)

# source("002_lfc.R")
source("003_dea.R")
source("004_etl.R")
#source("050_nems.R")

run_E2V2_pipeline <- function() {

    ################################################################################
    #
    # Definitions
    #
    ################################################################################

    # file locations
    base_dir = "/scratcha/fmlab/baker02/ER_Model/E2V2/Initial_ER_Interactome_Network_Model/data"
    
    # output locations
    source_dir = file.path(base_dir, "000_source")
    feature_count_dir = file.path(base_dir, "002_lfc")
    diffexp_dir = file.path(base_dir, "003_diffexp")
    prepared_dir = file.path(base_dir, "004_prepared")
    nems_dir = file.path(base_dir, "005_nems")
#     consensus_dir = file.path(base_dir, "/data/060_consensus")
#     plots_dir = file.path(base_dir, "/data/070_plots")
#     benchmark_dir = file.path(base_dir, "/data/075_benchmark")
#     egenes_dir = file.path(base_dir, "/data/080_egenes")
    
    dirs_to_check = c(base_dir, feature_count_dir, diffexp_dir, prepared_dir, nems_dir)
    
    # create directory is the specified directory doesn't exists
    for (output_dir in dirs_to_check) {
        if( !dir.exists(output_dir)) {
            dir.create(output_dir)
        }
    }
    
    # aligners: hisat2, bowtie
    aligner = "hisat2"
    cell_lines = c('MCF7', 'T47D')
    project = 'Modeling ER Interactome of MCF7 and T47D'
    
    feature_count_fp = file.path(feature_count_dir, 'cleaned_feature_counts.Rds')
    metadata_fp = file.path(base_dir, "metadata/clean_metadata.csv")

# prep_methods: binary, lfc, pvalue
#     binary - discretize expression data to 0 or 1 (there is an effect or not)
#     pvalue - use pvalue of differential expression
    
    prep_method = "binary"
    
# filter_method: How do you want to select E-Genes: logfc_aprior, fdr_aprior, or custom
#     fdr_aprior - use FDR (Adjusted P-Values) to select to ranked E-Genes
#     logfc_aprior - use logfc (Adjusted P-Values) to select to ranked E-Genes
#     custom - Custom Gene set
    
    filter_method = 'fdr_aprior'
    genes_to_select = c()
    
# Which NEM methods should be applied?
# The way the data is prepared determines which NEM methods are compatible with which prepared
    
# Traditional NEM Methods: search, triples, pairwise, ModuleNetwork.orig, ModuleNetwork
#    search -- Markowetz 2005 "Non-transcriptional pathway features reconstructed from secondary effects of RNA interference"
#    triples -- Markowetz 2007 "Nested effects models for high-dimensional phenotyping screens"
#    pairwise -- Markowetz 2007 "Nested effects models for high-dimensional phenotyping screens"
#    nem.greedy -- H. Fröhlich, A. Tresch, T. Beißbarth, Nested Effects Models for Learning Signaling Networks from Perturbation Data
#    ModuleNetwork.orig -- Frohlich 2007 "Large scale statistical inference of signaling pathways from RNAi and microarray data"
#    ModuleNetwork -- Frohlich 2008, "Estimating large-scale signaling networks through nested effect models with intervention effects from microarray data"
    nem_method_compat = list("binary" = c("triples", "pairwise", "nem.greedy", "ModuleNetwork", "ModuleNetwork.orig"),
                             "pvalue" = c("triples", "pairwise", "nem.greedy", "ModuleNetwork", "ModuleNetwork.orig"))

    nem_method = "triples"
    
    # logfc and the adjusted p-values cut off points
    logfc_cutoff = 0.25
    padj_cutoff = 0.05

    
    report_attached_egenes = FALSE

    ################################################################################
    #
    # Pipeline
    #
    ################################################################################

    print("Running ER pipeline")
    print(paste("project: ", project))
    print(paste("aligner: ", aligner))

    # turn bam files into log fold change
    # step_002_lfc(project, aligner, experiment_definitions, base_input_dir, lfc_dir)
    
#     for (cell_line in cell_lines) {
#         # calculate differential expression per cell line
#         step_003_dea(feature_count_fp, metadata_fp, cell_line, diffexp_dir, run_parallel=TRUE)
#     }
   
    for (cell_line in cell_lines) {
        # prepare data in correct format to use with NEMs
        results_fp = file.path(diffexp_dir, paste0('deseq2_results_', cell_line, ".Rds"))
        step_004_etl(cell_line, results_fp, prepared_dir, prep_method, filter_method, 
                     pval_threshold=padj_cutoff, logfc_threshold=logfc_cutoff, genes_to_keep=genes_to_select)
    }
    
#     for (cell_line in cell_lines) {
#         # create filepath for egene matrix and bootstrap list
#         egene_mat_fp = file.path(prepared_dir, paste(cell_line, prep_method, 'gene_mat.Rds', sep="_"))
#         boot_list_fp = file.path(prepared_dir, paste(cell_line, prep_method, 'boot_list.Rds', sep="_"))
        
#         step_005_nems(egene_mat_fp, prep_method, cell_line, nem_method, nem_method_compat,  nems_dir)
#         step_005_nems(boot_list_fp, prep_method, cell_line, nem_method, nem_method_compat,  nems_dir)
#     }
}

run_E2V2_pipeline()
