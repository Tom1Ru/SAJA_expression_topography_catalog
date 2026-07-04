#' Order gene IDs by apcluster label (C01, C02, …); within cluster, gene_id order.
order_genes_by_cluster <- function(cluster_df, gene_col = "gene_id", cluster_col = "cluster") {
  cluster_num <- as.integer(sub("^C0*", "", cluster_df[[cluster_col]]))
  ord <- order(cluster_num, cluster_df[[gene_col]])
  as.character(cluster_df[[gene_col]][ord])
}

load_DEGn3_genes_cluster_order <- function(catalog_root_dir = catalog_root()) {
  paths <- catalog_paths(catalog_root_dir)
  cluster_df <- utils::read.csv(paths$cluster, stringsAsFactors = FALSE)
  order_genes_by_cluster(cluster_df)
}
