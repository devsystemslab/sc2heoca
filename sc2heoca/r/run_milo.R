library(igraph)
library(data.table)

run_milo <- function(t_milo, adata_query, group_by){
    set.seed(0)
    o_milo <- Milo(adata_query)
    o_milo <- buildGraph(o_milo, k=30, d=2, reduced.dim="UMAP")
    o_milo <- makeNhoods(o_milo, prop=0.05, k=30, d=2, refined=T, reduced_dims="UMAP")
    o_milo <- buildNhoodGraph(o_milo)
    assay(o_milo, "logcounts") <- assays(o_milo)[['X']]

    out <- scrabbitr::calcNhoodSim(t_milo, o_milo, 
                                   orthologs= cbind(rownames(o_milo), rownames(o_milo)),
                                   sim_preprocessing="gene_spec", 
                                   sim_measure="spearman",
                                   verbose = FALSE
                                  )

    # Calculate maximum correlations  
    out$nhood_sim[is.na(out$nhood_sim)] <- 0

    t_maxNhoods <- scrabbitr::getMaxMappings(out$nhood_sim, 1, long_format=FALSE) # rabbit-mouse
    o_maxNhoods <- scrabbitr::getMaxMappings(out$nhood_sim, 2, long_format=FALSE) # mouse-rabbit
    df_simFilt <- rbind(t_maxNhoods, o_maxNhoods)

    graph <- miloR::nhoodGraph(t_milo)

    df <- data.frame(nhood=igraph::vertex_attr(graph)$name,
                   max_sim=t_maxNhoods$sim,
                   group=colData(t_milo)[as.numeric(vertex_attr(graph)$name), group_by])
    return(df)
}
