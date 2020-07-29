step_003_dea = function(feature_count_fp, metadata_fp, cell_line, output_dir, cell_line_col='Cell.line', 
                        condition_col='Gene', sample_col='sample_names', ctrl_label='CTRL', run_parallel=FALSE,
                        logfc_col='log2FoldChange', pval_col='pvalue', fdr_col='padj') {
    
    # load metadata and feature count matrix.
    metadata = load_metadata(metadata_fp, cell_line, ctrl_label, cell_line_col, condition_col, sample_col)
    feature_counts = load_feature_count(feature_count_fp, metadata, sample_col)
    
    # conducting Differential Expression analyis and assessing Perturbation Quality.
    print (paste0("working on cell line: ", cell_line))
    
    results = diff_expr_analysis(feature_counts, metadata, cell_line, output_dir, condition_col, ctrl_label, run_parallel, logfc_col, pval_col, fdr_col)
    assess_ko_quality(results, cell_line, metadata, names(results), output_dir)
    
}

load_feature_count = function(feature_count_fp, metadata, sample_col){
    # check feature count matrix file paths is valid.
    stopifnot(file.exists(feature_count_fp))
    
    # load Feature Count Matrix
    feature_counts = readRDS(feature_count_fp)
    feature_counts = feature_counts$counts
    
    # selecting samples that are only in the metadata. then make sure the matrix isn't empty
    feature_counts = feature_counts[ , colnames(feature_counts) %in% metadata[[sample_col]]]
    stopifnot(ncol(feature_counts) >= 2)
    feature_counts = feature_counts[!duplicated(rownames(feature_counts)), ]
    
    return(feature_counts)
}

load_metadata = function(metadata_fp, cell_line, ctrl_label, cell_line_col, condition_col, sample_col){
    
    # check file path if valid load and check metadata
    stopifnot(file.exists(metadata_fp))
    metadata = read.csv(metadata_fp)
    metadata_cols = colnames(metadata)
    
    # checking if the user provided columns are in the metadata.
    stopifnot(cell_line_col %in% metadata_cols)
    stopifnot(condition_col %in% metadata_cols)
    stopifnot(sample_col %in% metadata_cols)
    
    # checking that cell line and ctrl condition labels are in the metadata.
    stopifnot(cell_line %in% metadata[[cell_line_col]])
    stopifnot(ctrl_label %in% metadata[[condition_col]])
    
    # extracting metadata related to the cell_line being analyzed
    metadata = metadata[metadata[[cell_line_col]] == cell_line, ]
    
    # resetting the factor for the condition column 
    cond = factor(metadata[[condition_col]])
    cond = relevel(cond, ctrl_label)
    metadata[[condition_col]] = cond
    
    return (metadata)
}

################################## Utility Function for Differential Expression Analysis #####################################

diff_expr_analysis = function(feature_counts, metadata, cell_line, output_dir, condition_col, ctrl_label, run_parallel, logfc_col, pval_col, fdr_col) {
    # Differential Expression Analysis with DESeq2 is a two step process. First estimate gene expression distributions and then second extract 
    # the results from DEA. 
    dds_filename = paste0('deseq2_model_', cell_line, '.Rds')
    dds_fp = file.path(output_dir, dds_filename)

    results_filename = paste0('deseq2_results_', cell_line, ".Rds")
    results_fp = file.path(output_dir, results_filename)

    if (!file.exists(dds_fp)){
        # a deseq 2 object doesn't exists conducting DEA with deseq 2 now
        dds = running_deseq(feature_counts, metadata, dds_fp, run_parallel)
    } else {
        # detected a deseq 2 object loading it. no point rerunning deseq 2
        dds = readRDS(dds_fp)
    }

    if (!file.exists(results_fp)){
        # Extract DEA results, set internal columns names, and then save the results object
        results = extract_results(dds, condition_col, ctrl_label, results_fp, run_parallel)
        results = set_internal_columns(results, logfc_col, pval_col, fdr_col)
        saveRDS(results, results_fp)
    
    } else {
        results = readRDS(results_fp)
    }
    
    return (results)
}

running_deseq = function(feature_counts, metadata, dds_fp, run_parallel=FALSE) {
    # Estimating the Expression Distributions of samples
    # create DESeq2 object from cell line count and meta data matrix.
    print ("Creating DESeq2 object...")
    dds = DESeq2::DESeqDataSetFromMatrix(countData=feature_counts,
                                          colData=metadata,
                                          design= ~ Gene)

    # estimate the size factors and run DESeq2
    print ("Estimating Size Factors...")
    dds = DESeq2::estimateSizeFactors(dds)
    
    print ("Running DESeq2...")
    dds = DESeq2::DESeq(dds, parallel=run_parallel)
    
    saveRDS(dds, dds_fp)
    return (dds)
}

extract_results = function(dds, condition_col, ctrl_label, results_fp, run_parallel=FALSE){
    
    results = list()
    conditions = dds[[condition_col]]
    selected_genes = unique(conditions[which(conditions != ctrl_label)])
    
    # fetching the results for all shRNA and applying shrinkage to there log2 fold change estiamtes.
    for(gene in selected_genes){
        gene_contrast = c(condition_col, gene, ctrl_label)
        gene_coefficient = paste0(condition_col, "_", gene, "_vs_", ctrl_label)

        print(paste0("Fetching the Results for ", gene))
        res_unshrunken = DESeq2::results(dds, alpha = 0.05, contrast = gene_contrast, parallel=run_parallel)
        results[[gene]] = DESeq2::lfcShrink(dds, gene_coefficient, res=res_unshrunken, type="normal", parallel=run_parallel)
    }
    
    return (results)
}

set_internal_columns = function(results, logfc_col, pval_col, fdr_col){
    # setting internal column names for downstream analysis. By creating these columns now.
    for (gene in names(results)){
        # checking that the same logfc, p-value, FDR (adjusted P-Value) columns are in each target
        # genes result dataframe.
        res = results[[gene]]
        
        stopifnot(logfc_col %in% colnames(res))
        stopifnot(pval_col %in% colnames(res))
        stopifnot(fdr_col %in% colnames(res))
        
        res$logfc = res[[logfc_col]]
        res$pval = res[[pval_col]]
        res$fdr = res[[fdr_col]]
        
        results[[gene]] = res
    }
    
    return (results)
}

################################## Visualizing Knockout Effect of Conditions (i.e. Target Gene) #####################################

assess_ko_quality = function(results, cell_line, metadata, selected_genes, output_dir){
    # plotting adjusted p-value
    padj = ko_matrix(results, cell_line, metadata, selected_genes, 'padj')
    plot_ko_matrix(padj, 'padj', cell_line, output_dir)

    # plotting log 2 fold change
    logfc = ko_matrix(results, cell_line, metadata, selected_genes, 'log2FoldChange')
    plot_ko_matrix(logfc, 'logFold2Change', cell_line, output_dir)
}

ko_matrix = function(results, cell_line, metadata, selected_genes, stat='padj') {
    # setting global parameters
    df = data.frame()
    selected_genes = toupper(selected_genes)
    
    # creating heatmap matrix
    for (gene in selected_genes) {
        df = rbind(df, results[[gene]][selected_genes, stat])
    }
    
    # naming columns and rows after selected genes
    colnames(df) = selected_genes
    rownames(df) = selected_genes
    
    #group genes and guides
    f = df[colnames(df), colnames(df)]
    
    if (stat == 'padj') {
        # if an adjusted p-value remove 0.0 values and perform a log transformation
        min_pval = min(f[f > 0.0])
        f[f == 0.0] = min_pval
        
        f = -log(f, 10)                                
    }
    
    return (f)
}

plot_ko_matrix = function(ko_mat, stat, cell_line, output_dir) {
    # file path to save heatmap figure
    output_file_name = paste0(cell_line, "_", stat, "_heatmap", '.pdf')
    output_fp = file.path(output_dir, output_file_name)
    
    # creating a KO assessment heatmap
    plot_meta = ko_plot_metadata(stat)
    ko_heatmap = pheatmap(ko_mat, cluster_rows=FALSE, cluster_cols=FALSE, color=plot_meta$colors, breaks=plot_meta$breaks)
    ggplot2::ggsave(output_fp, ko_heatmap)
}

ko_plot_metadata = function(stat){
    if (stat == 'padj') {
        # colors for the adjusted p-values
        breaksup = seq(2, 5, 0.5)
        breaksdn = seq(0, 1.5, 0.5)
        colordn = rep("#FFFFFF", length(breaksdn))
        colorup = colorRampPalette(RColorBrewer::brewer.pal(n=7, name="Blues"))(length(breaksup))

        breaks = c(breaksdn, breaksup)
        colors = c(colordn, colorup)
    } else {
        # colors for the logfc change
        breaksdn = seq(-1.5, -0.5, 0.5)
        breaksup = seq(0.5, 3, 0.5)
        breaksdn = seq(-3, -0.5, 0.5)
        breaksup = seq(0.5, 2, 0.5)
        colordn = rev(colorRampPalette(RColorBrewer::brewer.pal(n=7, name="Reds"))(length(breaksdn)))
        colorup = colorRampPalette(RColorBrewer::brewer.pal(n=7, name="Blues"))(length(breaksup))

        breaks = c(breaksdn, 0, breaksup)
        colors = c(colordn, "#FFFFFF", colorup)
    }
    
    return (list(breaks=breaks, colors=colors))
}