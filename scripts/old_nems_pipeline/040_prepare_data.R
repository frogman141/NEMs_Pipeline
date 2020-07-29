# pipeline
# prepare the differential expression data for input to nems in a variety of ways
library(dplyr)
library(rlist)
source('ko_quality_control.R')

# functions
decompose <- function(shrink.list, k, adjusted_pvalue_cutoff, regulon) {
    # has to work on binary data because we use original nem method on 4 node subsets
    prepared <- binary(shrink.list, adjusted_pvalue_cutoff)
    filtered <- filter_regulon(prepared, regulon)
    n <- ncol(filtered)
    combinations <- combn(selected.genes, k)
    decomposed_list <- list()
    for (comb in 1:ncol(combinations)) {
        decomposed_list[[comb]] <- filtered[,combinations[,comb]]
    }
    return(decomposed_list)
}

binary <- function(shrink.list, adjusted_pvalue_cutoff) {
    # generate discretized matrix
    # (will be needed for downstream modeling of the NEMs)
    sig.names <- c()
    
    for(i in 1:length(shrink.list)){
        sig.names <- c(sig.names, rownames(shrink.list[[i]])[which(shrink.list[[i]]$padj < adjusted_pvalue_cutoff)])
    }
    print ("after first for loop")
    sig.names <- unique(sig.names)
    # discrete
    DESeq.disc.mat <- matrix(0L, 
                             nrow = length(sig.names), 
                             ncol= length(names(shrink.list))
                             )
    
    colnames(DESeq.disc.mat) <- names(shrink.list)
    rownames(DESeq.disc.mat) <- sig.names
    
    for(i in 1:length(colnames(DESeq.disc.mat))){
        DESeq.disc.mat[rownames(shrink.list[[i]])[which(shrink.list[[i]]$padj < adjusted_pvalue_cutoff)], i] <- 1
    }
    return(DESeq.disc.mat)
}

cluster_bin <- function(d, k) {
    # cluster the input genes to try to reduce noise
    # use binary data that has already been filtered for the regulon
    d_toplot <- unique(d)
    if (cluster_method == "hierarchical") {
        library(dynamicTreeCut)
        library(lsa)
        c <- cosine(t(d))
        dissimilarity <- 0.5*(1-c)
        dista <- as.dist(dissimilarity)
        # cluster, and cut the tree 
        h <- hclust(dista, method="average")
        clusters <- cutreeDynamic(h, distM = as.matrix(dista), minClusterSize = 1, method = "tree", deepSplit = 1, verbose = 0)
        # get points as the centroid of the cluster
        listy <- list()
        for (idx in 1:max(clusters)) {
            listy[[idx]] <- colMeans(d[which(clusters == idx),])
        }
        clustered_data <- do.call(rbind, listy)
        # we now have no row names because the 'egenes' are centroids
        # it's kind of odd to take the centroids because they are no longer binary
        # TODO print the dendrogram and the cluster definitions

        k <- 0

        return(clustered_data)
    } else if (cluster_method == "kmeans") {
        # attempted to use spherical kmeans, but this doesn't improve the silhouette measure of how appropriately clustered things are
        #library(skmeans)
        #kmeans_out <- skmeans(d_toplot,k)
        kmeans_out <- kmeans(d_toplot,k)
        clustered_data <- kmeans_out$centers
        clusters <- kmeans_out$cluster
    }
    library(factoextra)
    fviz_cluster(kmeans_out, data=as.data.frame(d_toplot))
    library(cluster)
    # use silhouette to measure how well clustered the data is
    s <- silhouette(clusters, dist(d_toplot))
    fviz_nbclust(d_toplot, FUNcluster=skmeans, method="silhouette", k.max=100, nboot=10)


    colours <- clusters[rownames(d_toplot)]
    perplexities <- c(3, 10, 30)
    for (perp in perplexities) {
        #tsne_plot(d_toplot, cluster_method, k, perp, colours)
    }

    find_representatives <- TRUE
    if (find_representatives) {
        # find real data points that are closest to the centroids, and return these
        representatives <- list()
        idx <- 1
        for (row in 1:nrow(clustered_data)) {
            centroid <- clustered_data[row,]
            l <- sort(apply(d_toplot, 1, function(x) {cosine(centroid, x)})) 
            representative <- names(l)[[length(l)]]
            representatives[[idx]] <- d_toplot[representative,]
            idx <- idx + 1
        }
        return(do.call(rbind, representatives))
    }
    return(clustered_data)
}

tsne_plot <- function(d_toplot, cluster_method, k, perplexity, colours) {
    library(Rtsne)
    library(ggplot2)
    set.seed(42)
    tsne_out <- Rtsne(d_toplot, pca=FALSE, perplexity=perplexity, theta=0)
    output_file_name <- file.path(prepared_dir, paste("tsne", cluster_method, k, perplexity, prep_method, project, "pdf", sep="."))
    p <- as.data.frame(tsne_out$Y)
    colnames(p) <- c("x", "y")
    p$colour <- factor(unname(colours))
    ggplot(p, aes(x=x, y=y, col=colour)) + geom_point() + theme_bw() + coord_fixed() + theme(legend.position = "none") 
    ggsave(output_file_name)
}

bootstrap <- function(shrink.list, adjusted_pvalue_cutoff, regulon) {
    # bootstrap from top egenes determined by adjusted_pvalue_cutoff, take half each bootstrap, ten times
    # then compare these graphs
    # does the filtering by regulon, so this does not need to happen again
    # TODO: refactor other functions so the filtering happens before

    sig.names <- c()
    for(i in 1:length(shrink.list)){
        sig.names <- c(sig.names,rownames(shrink.list[[i]])[which(shrink.list[[i]]$padj < adjusted_pvalue_cutoff)])
    }
    sig.names <- unique(sig.names)
    # TODO make this dynamic
    
    top_200_genes <- matrix(0L, 
                             nrow = length(sig.names), 
                             ncol= length(names(shrink.list))
                             )
    colnames(top_200_genes) <- names(shrink.list)
    rownames(top_200_genes) <- sig.names

    # significantly DEGs with padj < 0.05 are set to 1
    for(i in 1:length(colnames(top_200_genes))){
        top_200_genes[rownames(shrink.list[[i]])[which(shrink.list[[i]]$padj < adjusted_pvalue_cutoff)],i] <- 1
    }
    
    print(paste0("genes available before filtering: ", nrow(top_200_genes)))
    top_200_genes <- filter_regulon(top_200_genes, regulon)

    #bootstraps
    n_bootstrap <- 50
    #n_rows <- nrow(top_200_genes) / 2  # number of rows to bootstrap is half (???) of the available rows
    n_rows <- nrow(top_200_genes) # number of rows to bootstrap is all of the available rows, sampled with replacement
    print(paste0("number of genes sampled: ", n_rows, " from ", nrow(top_200_genes)))
    bootstrap_list <- list()
    for(i in 1:n_bootstrap) {
        rows <- sample.int(n_rows, replace=TRUE)
        bootstrap_list[[i]] <- top_200_genes[rows,]

    }
    return(bootstrap_list)
}

random <- function(shrink.list) {
    # generate random data to compare the bootstraps against
    # then draw bootstraps from the randomly generated data
    rando <- sample(c(0,1), 1042*10, replace=TRUE) # TODO make size dynamic based on size of shrink.list after filtering, so that we get the right number
    rando <- matrix(rando, nrow=1042, ncol=10)

    n_bootstrap <- 50
    nrows <- 1042
    bootstrap_list <- list()
    for (i in 1:n_bootstrap) {
        # bootstrap all randomly generated data, sampled with replacement
        rows <- sample.int(1042, replace=TRUE)
        sampled <- rando[rows,]
        colnames(sampled) <- selected.genes
        rownames(sampled) <- paste("row", 1:1042, sep="")
        bootstrap_list[[i]] <- sampled
    }
    return(bootstrap_list)
}

progressive <- function(shrink.list, adjusted_pvalue_cutoff) {
    # take top 50, 100, 150 and so forth egenes
    
    # 
    all_pvalues <- list()
    for (i in 1:length(shrink.list)) {
        all_pvalues[[i]] <- sort(shrink.list[[i]]$padj)
    }
    p <- sort(unlist(all_pvalues))

    gap <- 200
    prev_length <- 0

    bootstrap_list <- list()
    # don't take any genes that we don't think are significant, and don't fall off the end of the p list
    #for (prog in seq(50, min(nrow(shrink.list[[1]])/length(names(shrink.list)), max(which(p < adjusted_pvalue_cutoff))), 100)) {
    for (prog in seq(gap, max(which(p < adjusted_pvalue_cutoff)), gap)) {
        test <- prog * length(names(shrink.list))
        if(test > length(p)) {
            next
        }
        # guess at an appropriate pvalue cutoff to give us this number of genes across all of the experiments
        pval_cutoff <- p[[test]]

        sig.names <- c()
        for(i in 1:length(shrink.list)){
            sig.names <- c(sig.names,rownames(shrink.list[[i]])[which(shrink.list[[i]]$padj < pval_cutoff)])
        }
        sig.names <- unique(sig.names)

        top_n_genes <- matrix(0L, 
                                 nrow = length(sig.names), 
                                 ncol= length(names(shrink.list))
                                 )
        colnames(top_n_genes) <- names(shrink.list)
        rownames(top_n_genes) <- sig.names


        for(i in 1:length(colnames(top_n_genes))){
            top_n_genes[rownames(shrink.list[[i]])[which(shrink.list[[i]]$padj < pval_cutoff)],i] <- 1
        }

        top_n_genes <- filter_regulon(top_n_genes, regulon)

        if (nrow(top_n_genes) - prev_length >= gap) {
            bootstrap_list[[as.character(nrow(top_n_genes))]] <- top_n_genes
            prev_length <- nrow(top_n_genes)
        }
    }
    return(bootstrap_list)
}

lfc <- function(shrink.list, adjusted_pvalue_cutoff) {
    # use log fold change instead of pvalue because we are more interested in the size of the effect than whether it is significant
    sig.names <- c()
    for(i in 1:length(shrink.list)){
        sig.names <- c(sig.names,rownames(shrink.list[[i]])[which(shrink.list[[i]]$padj < adjusted_pvalue_cutoff)])
    }
    sig.names <- unique(sig.names)

    DESeq.padj.mat <- matrix(0L, 
                             nrow = length(sig.names), 
                             ncol= length(names(shrink.list))
                             )
    colnames(DESeq.padj.mat) <- names(shrink.list)
    rownames(DESeq.padj.mat) <- sig.names
    # significantly DEGs with padj < 0.05 are set to padj
    for(i in 1:length(colnames(DESeq.padj.mat))){
        padj <- shrink.list[[i]][which(shrink.list[[i]]$padj < adjusted_pvalue_cutoff),]$log2FoldChange
        DESeq.padj.mat[rownames(shrink.list[[i]])[which(shrink.list[[i]]$padj < adjusted_pvalue_cutoff)],i] <- padj
    }
    return(DESeq.padj.mat)

}

pvalue <- function(shrink.list, adjusted_pvalue_cutoff) {
    # padj for mc.eminem
    sig.names <- c()
    for(i in 1:length(shrink.list)){
        sig.names <- c(sig.names,rownames(shrink.list[[i]])[which(shrink.list[[i]]$padj < adjusted_pvalue_cutoff)])
    }
    sig.names <- unique(sig.names)

    DESeq.padj.mat <- matrix(0L, 
                             nrow = length(sig.names), 
                             ncol= length(names(shrink.list))
                             )
    colnames(DESeq.padj.mat) <- names(shrink.list)
    rownames(DESeq.padj.mat) <- sig.names
    # significantly DEGs with padj < 0.05 are set to padj
    for(i in 1:length(colnames(DESeq.padj.mat))){
        padj <- shrink.list[[i]][which(shrink.list[[i]]$padj < adjusted_pvalue_cutoff),]$padj
        DESeq.padj.mat[rownames(shrink.list[[i]])[which(shrink.list[[i]]$padj < adjusted_pvalue_cutoff)],i] <- padj


                                                }
    return(DESeq.padj.mat)
    
}
    
geneset <- function(diffexp) {
    sigpathways <- c()
    keep <- list()
    for(sgene_idx in 1:length(diffexp)){
        d <- diffexp[[sgene_idx]] 
        d <- d[order(d$log2FoldChange),]
        # remove any rows that are named NA
        d <- d[!grepl("^NA\\.", rownames(d)),]
        # convert gene names to entrez ids
        entrez_genes <- unlist(mget(x=rownames(d),envir=org.Hs.egALIAS2EG))
        # only rows that are in entrez_genes
        d <- d[rownames(d) %in% names(entrez_genes),]

        # TODO there is a problem with duplicates somewhere
        lfc <- d[,'log2FoldChange']
        names(lfc) <- entrez_genes[unique(rownames(d))]

        # get pathways associated with these genes
        pathways <- reactomePathways(entrez_genes)

        fgseaRes <- fgsea(pathways = pathways, 
                          stats = lfc,
                          minSize=15,
                          maxSize=500,
                          nperm=1000)
        plotEnrichment(pathways[[fgseaRes[order(pval), ][1,]$pathway]], lfc) + labs(title=fgseaRes[order(pval), ][1,]$pathway)
        geneset_pval_threshold <- 0.005
        # get list of all significant pathways
        sigpathways <- unique(unlist(c(sigpathways, fgseaRes[fgseaRes$pval < geneset_pval_threshold,]$pathway)))
        keep[[sgene_idx]] <- fgseaRes
    }

    # create matrix for nems
    gsea_mat <- matrix(0L, 
                       nrow = length(sigpathways), 
                       ncol= length(names(diffexp))
                      )
    for(i in 1:ncol(gsea_mat)){
        fgseaRes <- keep[[i]]
        for (j in 1:length(sigpathways)) {
            gsea_mat[j,i] = fgseaRes[fgseaRes$pathway == sigpathways[j],]$pval
        }
    }
    rownames(gsea_mat) <- sigpathways
    colnames(gsea_mat) <- names(diffexp)
    return(gsea_mat)
}

filter_regulon <- function(prepared, regulon) {
    #filter for only genes that are found to be regulated by E2
    return(prepared[which(rownames(prepared) %in% regulon),])
}

# main
step_040_prepare_data <- function(project, aligner, diffexp_method, prep_method, diffexp_dir, prepared_dir, adjusted_pvalue_cutoff, regulon, exp_def) {
    matching <- dir(diffexp_dir, pattern=diffexp_method)
    perturb_by_gene <- groupPerturbationByGene(exp_def)
    
    # filter for project name
    matching <- Filter(function(x) grepl('raw_DESeq', x), matching)
    
    for (input_file_name in matching) {
        diffexp <- readRDS(file.path(diffexp_dir, input_file_name))
        
        ######### REMOVE THESE FUNCTIONS #########
        good_ko <- checkQualityOfKnockOut(perturb_by_gene, diffexp, adjusted_pvalue_cutoff)
        cleaned_diffexp <- cleaningDiffexpData(exp_def, good_ko, diffexp, adjusted_pvalue_cutoff)
        
        if (prep_method == "binary") {
            
            prepared <- binary(cleaned_diffexp, adjusted_pvalue_cutoff)    
            filtered <- filter_regulon(prepared, regulon)
            output_file_name <- paste(prep_method, 'raw.diffexp', input_file_name, sep=".")
            saveRDS(filtered, file.path(prepared_dir, output_file_name))
        }
        if (prep_method == "cluster_bin") {
            prepared <- binary(diffexp, adjusted_pvalue_cutoff)
            filtered <- filter_regulon(prepared, regulon)
            ks <- c(5, 15, 30, 100)
            for (k in ks) {
                clustered <- cluster_bin(filtered, k)
                output_file_name <- paste(prep_method, "rep", cluster_method, k, input_file_name, sep=".")
                saveRDS(clustered, file.path(prepared_dir, output_file_name))
            }
        }
        if (prep_method == "lfc") {
            prepared <- lfc(diffexp, adjusted_pvalue_cutoff)
            filtered <- filter_regulon(prepared, regulon)
            output_file_name <- paste(prep_method, input_file_name, sep=".")
            saveRDS(filtered, file.path(prepared_dir, output_file_name))
        }
        if (prep_method == "pvalue") {
            prepared <- pvalue(diffexp, adjusted_pvalue_cutoff)
            filtered <- filter_regulon(prepared, regulon)
            output_file_name <- paste(prep_method, input_file_name, sep=".")
            saveRDS(filtered, file.path(prepared_dir, output_file_name))
        }
        
        if (prep_method == "geneset") {
            # no filter
            prepared <- geneset(diffexp)
            output_file_name <- paste(prep_method, input_file_name, sep=".")
            saveRDS(prepared, file.path(prepared_dir, output_file_name))
        }
        if (prep_method == "decompose") {
            # filters
            prepared_list <- decompose(diffexp, decompose_k_nodes, adjusted_pvalue_cutoff, regulon)
            acc <- 1
            for (i in 1:length(prepared_list)) {
                output_file_name <- paste(paste(prep_method, decompose_k_nodes, acc, sep="_"), input_file_name, sep=".")
                saveRDS(prepared_list[[i]], file.path(prepared_dir, output_file_name))
                acc <- acc + 1
            }
        }
        
        ################ THESE METHODS CAN BE AND SHOULD BE RAN ON ALL DATA PROCESSING METHODS ###############
        if (prep_method == "bootstrap") {
            # filters
            prepared_list <- bootstrap(diffexp, adjusted_pvalue_cutoff, regulon)
            acc <- 1
            for (i in 1:length(prepared_list)) {
                output_file_name <- paste(paste(prep_method, acc, sep="_"), input_file_name, sep=".")
                saveRDS(prepared_list[[i]], file.path(prepared_dir, output_file_name))
                acc <- acc + 1
            }
        }
        if (prep_method == "random") {
            prepared_list <- random(diffexp)
            acc <- 1
            for (i in 1:length(prepared_list)) {
                output_file_name <- paste(paste(prep_method, acc, sep="_"), input_file_name, sep=".")
                saveRDS(prepared_list[[i]], file.path(prepared_dir, output_file_name))
                acc <- acc + 1
            }
        }
        if (prep_method == "progressive") {
            # filters
            prepared_list <- progressive(diffexp, adjusted_pvalue_cutoff)
            acc <- names(prepared_list)
            for (i in 1:length(prepared_list)) {
                output_file_name <- paste(paste(prep_method, acc[[i]], sep="_"), input_file_name, sep=".")
                saveRDS(prepared_list[[i]], file.path(prepared_dir, output_file_name))
            }
        }
        
    }
}










