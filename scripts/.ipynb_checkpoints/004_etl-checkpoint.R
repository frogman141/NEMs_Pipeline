# vector of nem methods currently implemented. 
# NEM Method used determines how I transform the data.
implemented_methods = c('binary', 'pvalue')
required_columns = c('pval', 'fdr', 'logfc')
egene_selection_methods = c('custom', 'fdr_aprior', 'logfc_aprior')

step_004_etl = function(cell_line, results_fp, output_dir, prep_method, filter_method='fdr_aprior', 
                        pval_threshold=0.05, logfc_threshold=1.5, genes_to_keep=c(), ntop=100) {
    # main function for step 005 Extraction, Transformation, and Loading
    results = load_dea(results_fp)
    qc_params(results, output_dir, prep_method, filter_method)
    
    # load, clean, and filter degs
    mats = transform_data(results, prep_method, pval_threshold, output_dir)
    mats = egenes_selection(mats, filter_method, genes_to_keep, ntop)
    
    egene_mat_fp = file.path(output_dir, paste(cell_line, prep_method, filter_method, 'gene_mat.Rds', sep="_"))
    saveRDS(mats, egene_mat_fp)
}

load_dea = function(results_fp){
    stopifnot(file.exists(results_fp))
    results = readRDS(results_fp)
    
    return (results)
}

qc_params = function(results, output_dir, prep_method, filter_method){
    # checking if directories exist, if user is filtering degs, if the prep_method is implemented
    stopifnot(dir.exists(output_dir))
    stopifnot(prep_method %in% implemented_methods)
    stopifnot(filter_method %in% egene_selection_methods)
    
    for (res in results){
        stopifnot(required_columns %in% colnames(res))
    }
}

################# Transform Differential Expression Analysis results list to a Matrix with E-Genes x S-Genes #################

transform_data = function(results, prep_method, pval_threshold, output_dir){
    mats = list()
    egenes = extract_egenes(results)
    
    mats$fdr = build_matrix(results, egenes, 'fdr')
    mats$pval = build_matrix(results, egenes, 'pval')
    mats$logfc = build_matrix(results, egenes, 'logfc')
    mats$data = transform_matrix(mats, prep_method, pval_threshold, output_dir)
    
    return (mats)
}

extract_egenes = function(results){
    # checking that the same logfc, p-value, FDR (adjusted P-Value) columns are in each target genes result dataframe.
    egenes = c()

    for (sgene in names(results)){
        res = results[[sgene]]
        res = clean_outliers(res)
        degs = rownames(res)

        if (length(egenes) == 0){
            egenes = degs
        } else {
            egenes = intersect(egenes, degs)
        }
    }
    
    return (egenes)
}

build_matrix = function(results, egenes, col_used){
    mat = empty_matrix(egenes, names(results))
    mat = fill_matrix(mat, results, egenes, col_used)
    
    return (mat)
}

transform_matrix = function(mats, prep_method, pval_threshold, output_dir){
    # transform egene matrix into either a binarized or log density version of itself.
    if (prep_method == 'pvalue') {
        # if prep method is pvalue then i need to generate a log density matrix from those p-values
        qc_dir = file.path(output_dir, "pvalue_diagnostic_plots")
        dir.create(qc_dir, showWarnings=FALSE)
        data_mat = logdensity_matrix(mats, qc_dir)
    } else {
        # binarize inputted data
        data_mat = binarize_matrix(mats$fdr, pval_threshold)
    }
    
    return (data_mat)
}

clean_outliers = function (res){
    # Remove null values and set fdr/pvals that are at 0 to the min pval or fdr.
    # I am trying to remove all unnecessary values and perserve the overall structure of pval and fdr distribution
    res = res[!is.na(res$fdr), ]
    num_of_ones = nrow(res[res$fdr == 1, ])
    num_of_zeros = nrow(res[res$fdr == 0,])
    
    if (num_of_ones != 0) {
        res[res$fdr == 1, ]$fdr = max(res$fdr[res$fdr != 1])
        res[res$pval == 1, ]$pval = max(res$pval[res$pval != 1])
    }
    
    if (num_of_zeros != 0) {
        res[res$fdr == 0, ]$fdr = min(res$fdr[res$fdr != 0])
        res[res$pval == 0, ]$pval = min(res$pval[res$pval != 0])
    }
    
    return (res)
}

empty_matrix = function(degs, target_genes){
    # create an empty matrix based on the number of degs x number of target gene
    mat = matrix(0L, nrow=length(degs), ncol=length(target_genes))
    
    rownames(mat) = degs
    colnames(mat) = target_genes
    
    return (mat)
}

fill_matrix = function(mat, results, egenes, col_used){
    
    for (sgene in names(results)){
        res = results[[sgene]]
        mat[egenes, sgene] = res[egenes, ][[col_used]]
    }
    
    return (mat)
}

logdensity_matrix = function(mats, qc_dir){
    # removing extreme P-Values (defined as 1 and 0 P-Values)
    mat = mats$pval
    
    for (sgene in colnames(mat)){
        egenes = mat[, sgene]
        egenes[egenes == 0] = min(egenes[egenes != 0])
        egenes[egenes == 1] = max(egenes[egenes != 1])
        mat[, sgene] = egenes
    }
    
    mat = nem::getDensityMatrix(mat, dirname=qc_dir)
    return (mat)
}

binarize_matrix = function(mat, pval_threshold) {
    
    for (sgene in colnames(mat)){
        egenes = mat[, sgene]
        egenes[egenes < pval_threshold] = 1
        egenes[egenes != 1] = 0
        mat[, sgene] = egenes
    }
    
    return (mat)
}

################# Filtering E-Genes from E-Genes x S-Genes Matrix #################
egenes_selection = function(mats, filter_method, genes_to_keep, ntop){
    
    if (filter_method == 'logfc_aprior') {
        mats = logfc_aprior_selection(mats, ntop)
    } else if (filter_method == 'fdr_aprior'){
        # use nem packages filterEGenes to select Statitically Signficant E-Genes (uses FDR to select top n genes)
        mats = fdr_aprior_selection(mats, ntop)
    } else if (filter_method == 'custom') {
        mats = custom_selection(mats, genes_to_keep)
    } 
    
    return (mats)
}

fdr_aprior_selection = function(mats, ntop=100){
    # This function implements the ad hoc a prior filtering of E-Genes using FDR Frohlich et al. 2008
    # 1. For each S-Gene rank order it's E-Genes by FDR (Adjusted P-Value).
    # 2. Select top n genes (default 100) and save to top_ranked_egenes vector
    # 3. Select all E-Genes saved to top_ranked_egenes in transformed E-Genes x S-Gene Matrix
    selected = filterEGenes(mats$pval, mats$data, mats$fdr, ntop=ntop)
    mats$data = selected$dat
    mats$misc = selected
    
    return (mats)
}

logfc_aprior_selection = function(mats, logfc_threshold, ntop=100) {
    # This function implements the ad hoc a prior filtering of Frohlich et al. 2008
    # 1. Filter logfc matrix using user specified logfc_threshold (default 1.5 log fold change)
    # 2. For each S-Gene rank order it's E-Genes by Log Fold Change.
    # 3. Select top n genes (default 100) and save to top_ranked_egenes vector
    # 4. Select all E-Genes saved to top_ranked_egenes in transformed E-Genes x S-Gene Matrix
    logfc_mat = mats$logfc
    top_ranked_egenes = c()
    
    for (sgene in colnames(logfc_mat)){
        egenes = logfc_mat[, sgene]
        egenes = egenes[egenes < abs(logfc_threshold)]
        top_egenes = order(egenes, decreasing=TRUE)[1:ntop]
        top_ranked_egenes = union(top_ranked_egenes, top_egenes)
    }
    
    mats$data = mats$data[top_ranked_egenes, ]
    return (mats)
}

custom_selection = function(mats, genes_to_keep){
    # Custome Gene Selection based on whatever the user defines as genes of interest.
    # Purpose of this is to allow the user to define a E-Genes nems use to reconstruct networks.
    degs = rownames(mats$data)
    mats$data = mats$data[degs %in% genes_to_keep, ]
    
    return (mats)
}

######################## UTILITY FUNCTION FOR E-GENE SELECTION ########################

filterEGenes = function(Porig, D, Padj=NULL, ntop=100, fpr=0.05, adjmethod="bonferroni", cutoff=0.05){
    # filterEGenes was originally from the nems package. Unfortunately, the nems package version of this function contains a bug when it 
    # filting our patterns with that have 0 across all columns. In this version we supply usie the ! to invert the boolean values to select rows
    # in the pattern matrix. 
    
    if(is.null(Padj)){
      Padj = apply(Porig, 2, p.adjust, method="fdr")
    }

    n = ncol(Porig)
    ntop = min(nrow(Porig), ntop)
    I1 = apply(Porig, 2, function(x) order(x)[1:ntop])
    I1 = unique(as.vector(I1))

    print(paste("Selecting top",ntop," genes from each list -->",length(I1),"genes total"))

    disc = (Padj[I1,] <= fpr) * 1
    N = nrow(disc)
    nsig = colSums(Padj[I1,] < fpr)

    patterns = unique(disc)
    patterns = patterns[!apply(patterns, 1, function(r) all(r == 0)), ] # bug is right here in the nem package. it uses -which() this erases the matrix

    if(nrow(patterns) < 1){i
        stop("No patterns found!")
    }

    idx = apply(patterns,1, function(p){
        cl = which(apply(disc,1, function(r) all(r == p)))
    })

    nobserved = sapply(idx, length)
    patterns = patterns[nobserved > 0,]
    idx = idx[nobserved > 0]
    nobserved = nobserved[nobserved > 0]

    cat("Testing ", nrow(patterns), " patterns\n")

    p.values = sapply(1:nrow(patterns), function(j){
        p = patterns[j,]
        pexpected = max(0,prod((fpr * p * nsig + (1 - fpr) * (1 - p) * (N - nsig))/ N))
        p.value = binom.test(nobserved[j], N, pexpected, alternative="greater")$p.value
        cat("pattern ", p, ": (#observed = ", nobserved[j], ", #expected = ", floor(pexpected*N), ", raw p-value = ", p.value,")\n")
        p.value
    })

    cat("\n")

    p.values = p.adjust(p.values,method=adjmethod)
    if(!any(p.values < cutoff)) {
        stop("No significant patterns found!\n")
    }

    patterns = patterns[p.values < cutoff,]	
    idx = idx[p.values < cutoff]
    nobserved = nobserved[p.values < cutoff]
    p.values = p.values[p.values < cutoff]

    I = I1[unlist(idx)]	
    cat(length(p.values), " significant patterns -->", length(I), "E-genes in total\n")
    D = D[I,]	

    return (list(selected=I, dat=D, patterns=patterns, nobserved=nobserved, p.values=p.values))
}