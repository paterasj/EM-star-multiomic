# Title: BRCA Subtype Differentially Expressed Gene Analysis 
# Description: Conduct differentially expressed gene analysis based on EM* subtyping for BRCA Patients 
# Author: Musaddiq Lodi 
# Last Update: 9/7/2023 

 # Load all necessary packages 
library(CancerSubtypes)
library(RTCGA.mRNA)
library(limma)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DOSE)
library(clusterProfiler)
library(dplyr)
library(GOSemSim)
library(ggplot2)
library(ggrepel)
library(patchwork)



setwd("INSERT DESIRED DIRECTORY HERE")

# Load in normal and diseased as .rda file (or dataframe format of RNA-seq expression matrix, patients as columns and genes as rows)
load('0802_new_geneexp_files/brca_normal.rda')
load('0802_new_geneexp_files/brca_tumor.rda')


# These groups are the subtypes for each patient informed by the EM* Analysis. Please load as a dataframe column, and follow conversion 
new_group <- list(2,4,3,4,2,4,4,3,4,4,4,4,4,4,2,2,4,4,4,4,4,2,4,4,4,4,4,4,2,4,4,4,2,4,4,1,4,3,4,4,4,4,4,1,4,4,4,4,4,2,3,4,4,4,4,4,4,4,1,2,4,4,4,3,4,4,4,2,4,4,4,4,4,4,2,4,4,1,4,3,1,3,4,4,4,4,4,1,4,4,3,4,4,4,2,4,1,4,4,1,4,4,1,1,3,4,4,1,4,2,4,4,4,4,4,4,3,2,4,4,4,4,3,4,1,4,4,4,4,3,4,4,3,2,4,2,2,2,4,4,4,1,4,4,3,2,4,4,4,2,4,1,4,1,4,2,4,4,4,2,3,3,4,2,1,4,4,3,2,1,1,4,4,2,2,4,2,3,4,4,1,4,4,4,4,2,3,4,4,2,4,4,4,2,3,4,4,2,4,4,1,4,4,2,4,2,2,2,4,4,4,3,2,2,2,2,3,4,4,1,2,2,3,2,4,2,4,1,2,2,2,2,1,2,4,1,4,2,4,1,4,2,3,3,4,1,4,4,4,4,2,1,3,2,4,4,2,1,2,1,4,1,2,2,2,2,4,1,4,3,1,2,3,2,1,4,1,4,4,4,4,4,3,4,4,2,4,4,2,4,4,3,3,2,4,4,4,2,3,2,1,1)
new_group_df <- as.data.frame(new_group)
t_new_group_df <- as.data.frame(t(new_group_df))

# Calculate DE Analysis using Limma 
result2 <- DiffExp.limma(Tumor_Data=brca_tumor,Normal_Data=brca_normal,group=t_new_group_df$V1,topk=NULL,RNAseq=TRUE)

# The top 50 genes are sorted in the result list index based on subtype. We are accessing the top 50 genes for each subtype and using them for further downstream analysis 

# TODO 9/7/2023: GENERALIZE THIS TO MAKE IT FOR VARIABLE NUMBER OF SUBTYPES 
result2[[1]]$ID <- gsub("\\..*","", result2[[1]]$ID)
result2[[2]]$ID <- gsub("\\..*","", result2[[2]]$ID)
result2[[3]]$ID <- gsub("\\..*","", result2[[3]]$ID)
result2[[4]]$ID <- gsub("\\..*","", result2[[4]]$ID)

Subtype1_gene <- as.vector(na.omit(result2[[1]]$ID[1:50]))
Subtype2_gene <- as.vector(na.omit(result2[[2]]$ID[1:50]))
Subtype3_gene <- as.vector(na.omit(result2[[3]]$ID[1:50]))
Subtype4_gene <- as.vector(na.omit(result2[[4]]$ID[1:50]))



create_volcano <- function(cancer_type, subtype_name, fc_table) {
gene_list <- fc_table$ID
annotations_orgDb <- AnnotationDbi::select(org.Hs.eg.db, # database
                                           keys = gene_list,  # data to use for retrieval
                                           columns = c("ENSEMBL", "SYMBOL","GENENAME"), # information to retreive for given data
                                           keytype = "ENSEMBL") # type of data given in 'keys' argument

# These parameters may be modified to color up/down regulated based on user preference. Consider adding this as a function argument with defaults 
# add a column of NAs
fc_table$diffexpressed <- "NO"
# if log2Foldchange > 3.5 and pvalue < 0.05, set as "UP" 
fc_table$diffexpressed[fc_table$logFC > 3.5 & fc_table$P.Value < 0.05] <- "UP"
# if log2Foldchange < -3.5 and pvalue < 0.05, set as "DOWN"
fc_table$diffexpressed[fc_table$logFC < -3.5 & fc_table$P.Value < 0.05] <- "DOWN"


fc_table$tolabel <- "NO"

fc_table$tolabel[fc_table$logFC > 6.5 & fc_table$P.Value < 0.05] <- "YES"
fc_table$tolabel[fc_table$logFC < -6.5 & fc_table$P.Value < 0.05] <- "YES"


fc_table$delabel <- NA
fc_table$delabel[fc_table$tolabel != "NO"] <- annotations_orgDb$SYMBOL[fc_table$tolabel != "NO"]

volcano_1 <- ggplot(data=fc_table, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=delabel)) + geom_point() + theme_minimal() + geom_text_repel()
volcano_1_cont <- volcano_1 + geom_vline(xintercept=c(-3.5, 3.5), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red") + ylim(0,80) + ggtitle(paste("Volcano Plot of Genes For BRCA", subtype_name))

volcano_1_fin <- volcano_1_cont + scale_color_manual(values=c("blue", "black", "red"))


ggsave(paste(cancer_type,"/DS_Analysis_Figures/Volcano",subtype_name,"/volcano.png"), plot = volcano_1_fin, device = "png", width = 10, height = 15)


return(volcano_1_fin)

}


# Will generate plots and write to appropriate directory 
subtype_1_volc <- create_volcano("08022023_BRCA", "Subtype_1", result2[[1]])
subtype_2_volc <- create_volcano("08022023_BRCA", "Subtype_2", result2[[2]])
subtype_3_volc <- create_volcano("08022023_BRCA", "Subtype_3", result2[[3]])
subtype_4_volc <- create_volcano("08022023_BRCA", "Subtype_4", result2[[4]])


combined_vin_plot_1 <- subtype_1_volc / subtype_2_volc
combined_vin_plot_2 <- subtype_3_volc / subtype_4_volc

final_combined_vio_plot <- subtype_1_volc / subtype_2_volc / subtype_3_volc / subtype_4_volc


ggsave(paste("INSERT DESIRED DIRECTORY HERE"), plot = final_combined_vio_plot, device = "png", width = 10, height = 15)

# Function to perform all DS analysis (GoSemSim and Pathway enrichment)
ds_analysis_function <- function(cancer_type, gene_list, result_table, subtype_name) {

    # Get gene list to do downstream analysis 
    annotations_orgDb <- AnnotationDbi::select(org.Hs.eg.db, # database
                                           keys = gene_list,  # data to use for retrieval
                                           columns = c("ENSEMBL", "SYMBOL","GENENAME"), # information to retreive for given data
                                           keytype = "ENSEMBL") # type of data given in 'keys' argument

    print("Finished annotations step")
    # Determine the indices for the non-duplicated genes
    # Return only the non-duplicated genes using indices
    non_duplicates_idx <- which(duplicated(annotations_orgDb$ENSEMBL) == FALSE)
    annotations_orgDb <- annotations_orgDb[non_duplicates_idx, ]
    print("Finished duplications step")
    ## Merge the annotations with the results
    res_ids <- inner_join(result_table, annotations_orgDb, by=c("ID"="ENSEMBL"))
    print("Finished merging step")
    # Determine significant overexpressed genes 
    allOE_genes <- as.character(res_ids$ENSEMBL)
    sigOE <- dplyr::filter(res_ids, adj.P.Val < 0.01)
    sigOE_genes <- as.character(na.omit(sigOE$ID))
    sigOE_id <- as.character(na.omit(sigOE$ID))
    sigOE_name <- as.character(sigOE$SYMBOL)
    print("Finished getting sig genes step")
    #print("PRINTING SIGOE_NAME")
    #print(sigOE_name)
    #print("printing sigOE_genes")
    #print(sigOE_genes)
    # GO Overenrichment analysis 
    ego_1 <- clusterProfiler::enrichGO(gene = sigOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 1,
                pvalueCutoff = 1, 
                readable = TRUE)
    #print("printing ego_1")
    #print(ego_1)
    print("Finished cluster enrichment step")
    # Create barplot 
    barplot_enrichment <- barplot(ego_1) + ggtitle(paste("Overrepresented Pathways in", subtype_name))

    # Create goplot 
    goplot_enrichment <- goplot(ego_1, showCategory = 10)

    print("Finished barplot and goplot creation step")
    ## SEMANTIC SIMILARITY ANALYSIS 
    
    # Crate matrix of semantic similarity 
    hsGO <- godata("org.Hs.eg.db", ont = "BP")
    print("Finished creating sem sim DB step")

    genes.df <- bitr(sigOE_name, "SYMBOL", "ENTREZID", "org.Hs.eg.db")
    genes <- genes.df[,1]
    eg <- genes.df[,2]
    names(genes) <- eg


    sim <- mgeneSim(eg, semData = hsGO, measure = "Lin", drop = "IEA", combine = "BMA")
    rownames(sim) <- genes[rownames(sim)]
    colnames(sim) <- genes[colnames(sim)]

    print("Finished generating sim matrix step")


    # Create GoSemSim matrix 
    simplot_enrichment <- DOSE::simplot(sim, labs = FALSE) + ggtitle(paste("Lin Semantic Similarity Between Biomarkers of", subtype_name))

    print("Finished simplot enrichment step")

    ## SAVE PLOTS 

    # Save barplot 
    ggsave(paste(cancer_type,"/DS_Analysis_Figures/",deparse(substitute(gene_list)),"/barplot.png"), plot = barplot_enrichment, device = "png", width = 10, height = 10)

    # Save goplot 
    ggsave(paste(cancer_type,"/DS_Analysis_Figures/",deparse(substitute(gene_list)),"/goplot.png"), plot = goplot_enrichment, device = "png", width = 10, height = 10)

    # Save gosemsim matrix 
    ggsave(paste(cancer_type,"/DS_Analysis_Figures/",deparse(substitute(gene_list)),"/simplot.png"), plot = simplot_enrichment, device = "png", width = 15, height = 15)

    print("Finished saving plot step")

    return(list(barplot_enrichment, goplot_enrichment, simplot_enrichment))
}

subtype_1_plots <- ds_analysis_function("08022023_BRCA ", Subtype1_gene, result2[[1]], "Subtype_1")
subtype_2_plots <- ds_analysis_function("08022023_BRCA ", Subtype2_gene, result2[[2]], "Subtype_2")
subtype_3_plots <- ds_analysis_function("08022023_BRCA ", Subtype3_gene, result2[[3]], "Subtype_3")
subtype_4_plots <- ds_analysis_function("08022023_BRCA ", Subtype4_gene, result2[[4]], "Subtype_4")


combined_barplot <- subtype_1_plots[[1]] / subtype_2_plots[[1]] / subtype_3_plots[[1]] / subtype_4_plots[[1]]
combined_sim_plot <- subtype_1_plots[[3]] / subtype_2_plots[[3]] / subtype_3_plots[[3]] / subtype_4_plots[[3]]


ggsave(paste("INSERT DESIRED DIRECTORY HERE"), plot = combined_barplot, device = "png", width = 10, height = 15)
ggsave(paste("INSERT DESIRED DIRECTORY HERE"), plot = combined_sim_plot, device = "png", width = 10, height = 25)
