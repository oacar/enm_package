#' given a dendrogram, take the smallest child and cut it into 4 subclusters.
#' Assign the resulting cluster numbers to genes, using given dataframe
#' 
#' @param hcl hclust object, to calculate dendrogram
#' @param df the dataframe to map clusters back to genes. rownames of the dataframe needs to include the dendrogram labels
#' @param k number of clusters for dendrogram cut, optional. default = 4 
#' @return The dataframe joined with cluster assignments for some or all genes, depending on the input df 
get.subclusters.dendrogram <- 			# 
    function
(
 hcl,df, k=4	
)
{
    library(dendextend)
    dend  <- as.dendrogram(hcl)
    if(nleaves(dend[[1]])<nleaves(dend[[2]])){
      sub_dend <- dend[[1]]
    }else{
      sub_dend <- dend[[2]]
    }
    sub_clusters <- sub_dend %>% cutree(k=k)%>%as.data.frame()%>%rownames_to_column()%>%mutate(rowname=as.numeric(rowname)-1)%>%select(clusters='.', rowname)
    df_with_clusters <- df%>%rownames_to_column()%>%mutate(rowname=as.numeric(rowname))%>%left_join(sub_clusters)
    df_with_clusters
}
