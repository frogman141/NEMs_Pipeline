step_005_nems = function(egene_mat_fp, prep_method, filter_metod, cell_line, nem_method, nems_dir, bootstrap=FALSE, report_attached_egenes=FALSE) {
    # main function of the nems wrapper script
    
    if (bootstrap){
        nem_model = modelling_bootstrap(egene_mat_fp, prep_method, nem_method)
        model_fp = file.path(nems_dir, paste(cell_line, nem_method, 'boot_nem_models.Rds', sep="_"))
        
    } else {
        egene_mat = load_egene_mat(egene_mat_fp)
        nem_type = egene_mat_datatype(prep_method)
        nem_model = run_nems(egene_mat, nem_method, nem_type)
        
        if (report_attached_egenes) {
            nem_model = attach_egenes(nem_model, egenes_mat)
        }
        
        model_fp = file.path(nems_dir, paste(cell_line, nem_method, 'nem_models.Rds', sep="_"))
    }
    
    saveRDS(nem_model, model_fp)
}

modelling_bootstrap = function(egene_mat_fp, prep_method, nem_method) {
    # looping through all of the bootstraps generated egene matrix
    bootstrap_nems = list()
    bootstrap_egene_list = load_egene_mat(egene_mat_fp)
    
    for (boot_name in names(bootstrap_egene_list)) {
        egene_mat = bootstrap_egene_list[[boot_name]]
        
        egene_mat = load_egene_mat(egene_mat_fp)
        nem_type = egene_mat_datatype(prep_method)
        nem_model = run_nems(egene_mat, nem_method, nem_type)
        
        bootstrap_nems[boot_name] = nem_model
    }
    
    return (bootstrap_nems)
}

load_egene_mat = function(egene_mat_fp){
    # load the egene matrix or bootstrap list of egene matrix
    stopifnot(file.exists(egene_mat_fp))
    egene_mat = readRDS(egene_mat_fp)
    
    return (egene_mat)
}

egene_mat_datatype = function(prep_method){
    # based on how the E-Gene matrix was prepared will determine the inference method NEMs package needs to use.
    # for binarized matrices we need to mLL for log densities use CONTmLLBayers.
    
    if (prep_method == 'binary') {
        type = "mLL" 
    } else {
        type = "CONTmLLBayes"
    }
    
    return (type)
}

run_nems = function(egene_mat, nem_method, type, nboot=1) {
    # running NEMs
    print (paste0("constructing network models using ", nem_method))
    
    contr = c(0.15,0.05)
    hyper = nem::set.default.parameters(colnames(egene_mat), 
                               para=contr, 
                               type=type)
    
    b = nem::nem.bootstrap(egene_mat,
        inference=nem_method,
        control=hyper, 
        verbose=F,
        nboot=nboot)
    
    return(b)
}

# pipeline
# do nem calculations
attach_egenes = function(nem_model, egenes) {
    # we are going to get $graph, an adjacency matrix and $mappos, a list of lists of egenes
    # and want to turn this into a single matrix

    m = as(nem_model$graph, "matrix")
    sgenes = rownames(m)
    egenes = unique(unlist(nem_model$mappos))
    genes = append(sgenes, egenes)
    
    bigger = matrix(0L, nrow = length(genes), ncol= length(genes))
    colnames(bigger) = genes
    rownames(bigger) = genes
    bigger[sgenes, sgenes] = m[sgenes, sgenes]
    
    for (sg in sgenes) {
        at = nem_model$mappos[[sg]]
        bigger[sg, at] = 1
        bigger[at, sg] = 1
    }
    nem_model$graph = bigger
    return(nem_model)
}