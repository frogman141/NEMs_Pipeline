#!/home/baker02/miniconda3/envs/Renv/bin/Rscript
# filepath way above uses my R virtual environment to Run this script. Recommendation is to use a virtual environment
# select differentially expressed genes from the parsed bam lfc data

process_with_deseq <- function(counts, design, selected.genes, diffexp_dir, output_file_name) {
    
    deseq2_results <- list()
    rownames(design) <- colnames(counts)
    cell.lines <- unique(design$Cell.line)

    guides <- unique(design$shRNA)
    guides <- gsub('-', '.', guides)
    guides <- as.character(guides[ -which(guides == "Renilla") ])
    
    # loop through all of the unique samplse conditions (in current example unique sample conditions are cell lines).
    for (i in 1:length(cell.lines)){
      # selecting metadata and rows of the cell line currently being analyzed
        result <- list()
        cell.line <- as.character(cell.lines[[i]])
        print (paste0("working on cell line: ", cell.line))
        meta <- design[design$Cell.line %in% c(cell.line), ]
        cell.line.counts <- counts[ , colnames(counts) %in% rownames(meta)]
        meta$Gene <- factor(meta$Gene, levels=c("CTRL", levels(meta$Gene)[-which(levels(meta$Gene) == "CTRL")]))
        
        # create DESeq2 object from cell line count and meta data matrix.
        dds <- DESeq2::DESeqDataSetFromMatrix(countData=cell.line.counts,
                                              colData=meta,
                                              design= ~ shRNA + Gene)

        # estimate the size factors and run DESeq2
        # ddf = DESeq2::collapseReplicates(dds, meta$Sample.name)
        dds <- DESeq2::estimateSizeFactors(dds)
        dds <- DESeq2::DESeq(dds, parallel=TRUE)

        # cache raw deseq2 object of the cell line for further use
        dds_object_name <- paste0('raw_deseq2_', cell.line, '_genes.RDS')
        saveRDS(dds, file.path(diffexp_dir, dds_object_name))
        # dds <- readRDS(file.path(diffexp_dir, dds_object_name))

        # fetching the results for all shRNA and applying shrinkage to there log2 foldchange estiamtes.
        for(gene in selected.genes){
            gene_contrast <- c("Gene", gene, "CTRL")
            gene_coefficient <- paste0("Gene_", gene, "_vs_CTRL")

            print(paste0("Fetching the Results for ", gene))
            res_unshrunken <- DESeq2::results(dds, alpha = 0.05, contrast = gene_contrast, parallel=TRUE)
            result[[gene]] <- DESeq2::lfcShrink(dds, gene_coefficient, res=res_unshrunken, type="apeglm", parallel=TRUE)
        }
        names(result) <- selected.genes
        #save cell line name and results to deseq_results list
        dds_results_name <- paste0('raw_DESeq_', cell.line, "_results_genes.Rds")
        # result <- readRDS(file.path(diffexp_dir, dds_results_name))
        deseq2_results[[i]] <- list(cell.line, result)
        saveRDS(result, file.path(diffexp_dir, dds_results_name))
    }
    
    return(deseq2_results)
}

guide_heatmap <- function(diffexp, diffexp_method, input_file_name, selected.genes, guides, cellline, diffexp_dir, experiment_definitions) {
    stat <- "log2FoldChange"
    selected.genes <- toupper(selected.genes)
    
    d <- data.frame()
    for (g in 1:length(selected.genes)) {
        d <- rbind(d, diffexp[[g]][selected.genes, stat])
    }

    colnames(d) <- selected.genes
    rownames(d) <- selected.genes

    e <- experiment_definitions %>% filter(Gene != "CTRL") %>% group_by(Gene, shRNA) %>% summarize() %>% print(n=200)
    e$shRNA <- gsub("-", ".", e$shRNA)

    #group genes and guides
    f <- d[as.character(e$Gene), sort(colnames(d))]
    # log2foldchange
    breaksdn <- seq(-1.5, -0.5, 0.5)
    breaksup <- seq(0.5, 3, 0.5)
    breaksdn <- seq(-3, -0.5, 0.5)
    breaksup <- seq(0.5, 2, 0.5)
    colordn <- rev(colorRampPalette(RColorBrewer::brewer.pal(n=7, name="Reds"))(length(breaksdn)))
    colorup <- colorRampPalette(RColorBrewer::brewer.pal(n=7, name="Blues"))(length(breaksup))

    ifn <- gsub("Rds", "pdf", input_file_name)
    output_file_name <- paste("gene_heatmap", cellline, diffexp_method, ifn, sep=".")
    pdf(file.path(diffexp_dir, output_file_name))
    pheatmap::pheatmap(f, cluster_rows=FALSE, cluster_cols=FALSE, color=c(colordn, "#FFFFFF", colorup), breaks=c(breaksdn, 0, breaksup))
    dev.off()
}

pval_heatmap <- function(diffexp, diffexp_method, input_file_name, selected.genes, guides, cellline, diffexp_dir, experiment_definitions) {
    # padj
    stat <- "padj"
    selected.genes <- toupper(selected.genes)
     
    d <- data.frame()
    for (g in 1:length(selected.genes)) {
        d <- rbind(d, diffexp[[g]][selected.genes, stat])
    }
    colnames(d) <- selected.genes
    rownames(d) <- selected.genes

    e <- experiment_definitions %>% filter(Gene != "CTRL") %>% group_by(Gene,shRNA) %>% summarize() %>% print(n=200)
    e$shRNA <- gsub("-", ".", e$shRNA)

    #group genes and guides
    f <- d[as.character(e$Gene), sort(colnames(d))]

    breaksup <- seq(2.5, 5, 0.5)
    breaksdn <- seq(0, 2, 0.5)
    colordn <- rep("#FFFFFF", length(breaksdn))
    colorup <- colorRampPalette(RColorBrewer::brewer.pal(n=7, name="Blues"))(length(breaksup))

    ifn <- gsub("Rds", "pdf", input_file_name)
    output_file_name <- paste("gene_pval_heatmap", cellline, diffexp_method, ifn, sep=".")
    pdf(file.path(diffexp_dir, output_file_name))
    pheatmap::pheatmap(-log(f,10), cluster_rows=FALSE, cluster_cols=FALSE, color=c(colordn, colorup), breaks=c(breaksdn, breaksup))
    dev.off()
    # rows are knockdowns
    # cols are lfc
}


# main

step_030_diffexp <- function(project, aligner, lfc_dir, diffexp_dir, experiment_definitions, expr.cutoff, samples, selected.genes) {
    print("030_diffexp")
    
    
    input_file_name <- paste(aligner, project, "Rds", sep=".")
    fc <- readRDS(file.path(lfc_dir, input_file_name))
    foo <- preprocess(fc, experiment_definitions, expr.cutoff, samples)
    counts <- foo[[1]] #x$counts
    experiment_definitions = foo[[2]]

    # extracting metadata from experiment definitions
    design <- experiment_definitions[,c(1,4,6,7,8,9)]
    
    if (diffexp_method == "DESeq") {
        # diffexp is now a list
        print ("starting anaylze data with DESeq2")
        # diffexp <- process_with_deseq(counts, design, selected.genes, diffexp_dir, "raw_deseq2_object")
        mcf7 = readRDS('../../NEMpipelineE2V2output/030_diffexp/raw_DESeq_MCF7_results.Rds')
      	t47d = readRDS('../../NEMpipelineE2V2output/030_diffexp/raw_DESeq_T47D_results.Rds')
      	diffexp = list(list('MCF7', mcf7), list('T47D', t47d))
      	print ("finished analyzing data with DESeq2")
    } else if (diffexp_method == "edgeR") {
        #diffexp <- process_with_edger(sel, controls, selected.genes)
    }

    # generate sanity check heatmap
    guides <- unique(design$shRNA)
    guides <- gsub('-', '.', guides)
    guides <- as.character(guides[ -which(guides == "Renilla") ])

    for (idx in 1:length(diffexp)) {
        cell.line.data <- diffexp[[idx]]
        cell.line <- cell.line.data[1]
        diffexp_data <- cell.line.data[2][[1]]
        
        print (paste0("Generating output for cell line ", cell.line))
        
        output_file_name <- paste(diffexp_method, cell.line, input_file_name, sep=".")
        guide_heatmap(diffexp_data, diffexp_method, input_file_name, selected.genes, guides,cell.line, diffexp_dir, experiment_definitions)
        pval_heatmap(diffexp_data, diffexp_method, input_file_name, selected.genes, guides, cell.line, diffexp_dir, experiment_definitions)
    }

    print ("Finished Conducting Differentianl Gene Expression Analysis")
}