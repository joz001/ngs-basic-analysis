library(clusterProfiler)
library(ggplot2)
library(readxl)
library(xlsx)
library(org.Hs.eg.db)

# make GO of data from DESeq output
{
  GO_summary <- function(res_data_frame, type){
    signif_res = res_data_frame %>% filter(!is.na(padj) & padj < 0.05)
    if (type == "up")
      signif_res = signif_res %>% filter(log2FoldChange > 1)
    else if (type == "down")
      signif_res = signif_res %>% filter(log2FoldChange < -1)
    else
      signif_res = signif_res %>% filter(log2FoldChange > 1 | log2FoldChange < -1)
    
    ego <- enrichGO(gene = signif_res$ENSG,
                    keyType = "ENSEMBL",
                    OrgDb = org.Hs.eg.db,
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE)
    if (length(ego$ID) < 1){
      return(NULL)
    }
    df = ego %>% dplyr::select(-"ID") %>% data.frame
    return(df[order(df$p.adjust),])
  }
  
  up_go = GO_summary(res_data_frame, "up") # 196
  down_go = GO_summary(res_data_frame, "down") # 79
  all_go = GO_summary(res_data_frame, "")
}
save_as = ''
write.xlsx(x = up_go, file = save_as, sheetName = "Upregulated Genes GO", append = TRUE)
write.xlsx(x = down_go, file = save_as, sheetName = "Downregulated Genes GO", append = TRUE)
write.xlsx(x = all_go, file = save_as, sheetName = "All Genes GO", append = TRUE)


# visualize top 10 up and down pathways 
{                                                      
  merged = list(down_go[1:10,], up_go[1:10, ])
  df = data.frame('Pathway' = c(merged[[1]]$Description, merged[[2]]$Description),
                  'negLogPAdj_Value' = -log10(c(merged[[1]]$p.adjust, merged[[2]]$p.adjust)),
                  'count' = c(merged[[1]]$Count, merged[[2]]$Count),
                  'direction_discrete' = factor(c(rep("Up-Regulated", nrow(merged[[1]])),
                                                  rep("Down-Regulated", nrow(merged[[2]]))),
                                                levels = c("Up-Regulated", "Down-Regulated")))
  
  ggplot(data=df, aes(x=reorder(Pathway, negLogPAdj_Value), 
                      y=negLogPAdj_Value, 
                      fill = direction_discrete)) +
    facet_grid(direction_discrete ~ ., scales = "free", space = "free_x") +
    geom_bar(stat="identity", show.legend = FALSE) +
    geom_text(aes(label = count, y = negLogPAdj_Value), size = 3, hjust = -0.5) +
    ggtitle("GRAS1 Knockdown Top 1000 Genes GO") +
    coord_flip() +
    labs(x = element_blank(), y = expression(Enrichment~(-log[10](P-Value)))) +
    scale_y_continuous(expand = c(0,0), 
                       limits = c(0,max(1.2 * df$negLogPAdj_Value)))+
    scale_fill_manual(values = c("red", "blue")) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title.x = element_text(size = 12))
}