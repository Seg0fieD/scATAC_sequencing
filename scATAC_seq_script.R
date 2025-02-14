#-------------------------------------------------------------------------------------------------
# Memory size allocation-------------------------------------------------------
paste0("Memory size = ", round(object.size(atacProj)/ 10^6, 3), "MB")

#--------------------------------------------------------------------------------------------------
# Library----------------------------------------------------------------------------------
suppressPackageStartupMessages({	
	library(ArchR)
	library(dplyr)
	library(ggplot2)
	library(patchwork)
	library(SingleCellExperiment)
	library(Seurat)
	library(SeuratWrappers)
	library(tidyverse)
	library(BiocManager)
	library(BSgenome)
	library(BSgenome.Hsapiens.UCSC.hg19)
	library(writexl)
	library(WriteXLS)
   })
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# set working directory  ----------------------------------------------------------------------
setwd("/Atac_project/")

# ArchR reference set-up ----------------------------------------------------------------------
addArchRThreads(threads = 4)
addArchRGenome("hg19")

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# List and Naming file----------------------------------------------------------------------------
atac1 <- list.files("/scATAC/datas",pattern = "*.tsv.gz", full.names = TRUE)
names(atac1) <- c("scATAC_BMMC_D5T1","scATAC_CD34_D7T1","scATAC_CD34_D8T1")

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Arrow file creation---------------------------------------------------------------------------------
Arrowx <- createArrowFiles(inputFiles = atac1,
                           sampleNames = names(atac1),
                           minTSS = 4,
                           minFrags = 1000,
                           addTileMat = TRUE,
                           addGeneScoreMat = TRUE)

Arrowx
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Doublet inference_Score -----------------------------------------------------------------------------------------
doubScores <- addDoubletScores(input = Arrowx,
                               k = 10, 
                               knnMethod = "UMAP", 
                               LSIMethod = 1)


#ArchR project Creation
atacProj <- ArchR::ArchRProject(ArrowFiles = Arrowx, 
                                outputDirectory = "/scATAC/Outputs",
                                copyArrows = TRUE, #save arrow files as arrow
                                showLogo = FALSE  )
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Saving the ArchR project file-----------------------------------------------------
saveArchRProject(ArchRProj = atacProj, outputDirectory = "/scATAC",load = FALSE)
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Load the ArchR project file-------------------------------------------------------  
            
atacProj <- ArchR::loadArchRProject(path = "/scATAC",showLogo = FALSE)
atacProj
Tab<-data.frame(getCellColData(atacProj))
write_xlsx(Tab, "MetaData_added.xlsx")
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Doublet Removal--------------------------------------------------------------
atacProj <- filterDoublets(ArchRProj = atacProj)
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Creating joint data structure of samples- adding meta-data-----------------------------------------
## Add Sample list 
Samples <- gsub("scATAC_","",atacProj$Sample)
atacProj$Samples <- Samples

## Add Type list
Type <-gsub("_D5T1|_D7T1|_D8T1","",atacProj$Samples)
atacProj$Type <- Type

## Add Donor list
Donor <- gsub("BMMC_D5T1","D5", 
              gsub("CD34_D7T1","D6",
                   gsub("CD34_D8T1","D7",atacProj$Samples)))
atacProj$Donor <- Donor

## Add sex list
sex <- gsub("D5","F", 
            gsub("D6","M",
                 gsub("D7","M",atacProj$Donor)))
atacProj$Sex <- sex
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Quality Control-------------------------------------------------------------------
## How many cells do you have for each sample 
cell_num <- data.frame(table(atacProj$Sample))

## plot of number of Fragments Vs TSS Enrichment 
df <- data.frame(getCellColData(ArchRProj = atacProj, 
                                select = c("log10(nFrags)", "TSSEnrichment")))
df

pl_fragvtss <- ggPoint( x = df[,1], y = df[,2],
                        colorDensity = TRUE,
                        continuousSet = "sambaNight",
                        xlabel = "Log10 Unique Fragments", 
                        ylabel = "TSS Enrichment",
                        xlim = c(log10(500), quantile(df[,1],probs = 0.99)),
                        ylim = c(0, quantile(df[,2],probs = 0.99)))+ 
  geom_hline(yintercept = 4, lty = "dashed") + 
  geom_vline(xintercept =3, lty ="dashed")
pl_fragvtss

plotPDF(pl_fragvtss, name = "TSS-vs-Frags_log(500).pdf", 
        ArchRProj = atacProj, addDOC = FALSE)

pl_nfragsvtss <- plotGroups(ArchRProj =  atacProj,
                            groupBy = "Sample",
                            colorBy = "CellColData",
                            name = c("nFrags","TSSEnrichment"),
                            plotAs = "")
plotPDF(pl_nfragsvtss, name = "TSS-vs-nFrags.pdf", ArchRProj = atacProj, addDOC = FALSE)

## Ridge plot for each sample- TSS Enrichment 
pl_each_tss <- plotGroups( ArchRProj = atacProj,
                           groupBy = "Sample",
                           colorBy = "cellColData",
                           name = "TSSEnrichment",
                           plotAs = "ridges",
)

pl_each_tss
plotPDF(pl_each_tss, name = "TSS-Enrishcment-scores-ridge-plot-for-each-sample.pdf", 
        ArchRProj = atacProj, addDOC = FALSE)

pl_each_frag <- plotGroups(ArchRProj = atacProj,
                           groupBy = "Sample",
                           colorBy = "CellColData",
                           name = "nFrags",
                           plotAs = "ridges")

plotPDF(pl_each_frag, name = "nFrags-ridge-plot-for-each-sample.pdf", 
        ArchRProj = atacProj, addDOC = FALSE)

pl_each_ndifrag <- plotGroups(ArchRProj = atacProj,
                              groupBy = "Sample",
                              colorBy = "CellColData",
                              name = "nDiFrags",
                              plotAs = "ridges")

plotPDF(pl_each_ndifrag, name = "nDiFrags-ridge-plot-for-each-sample.pdf", 
        ArchRProj = atacProj, addDOC = FALSE)

## Plotting fragment size distribution and TSS Enrichment profiles
pl_frag <- plotFragmentSizes(ArchRProj = atacProj)
pl_frag

pl_tss <- plotTSSEnrichment(ArchRProj = atacProj)
pl_tss

plotPDF(pl_frag,pl_tss, name = "QC-Sample-FragmentSizes-TSSProfile.pdf", 
        ArchRProj = atacProj,
        addDOC = FALSE, width = 8, height = 8)
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Dimensionality Reduction(LSI)-------------------------------------------------------
atacProj <- addIterativeLSI(ArchR = atacProj, useMatrix = "TileMatrix", 
                            name= "Iterative-LSI-1.0",
                            iterations = 2,
                            clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start= 10),
                            varFeatures = 25000,
                            dimsToUse= 1:30) 
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Batch Correction by Harmony--------------------------------------------------
atacProj <- addHarmony(ArchRProj = atacProj,
                       reducedDims = "Iterative-LSI-1.0", 
                       name = "Harmony",
                       groupBy =  "Sample")
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Clustering---------------------------------------------------------------------
## Adding clusters
atacProj <- addClusters(input= atacProj, reducedDims = "Iterative-LSI-1.0",
                        method ="Seurat", 
                        name = "Clusters", 
                        resolution = 0.2, force = TRUE)  

## accessing the clusters 
head(atacProj$Clusters)
cluster_tab <- data.frame(table(atacProj$Clusters))

## creating confusion matrix
conMat<- data.frame(confusionMatrix(paste0(atacProj$Clusters), paste0(atacProj$Samples)))

## Confusion matrix to heatmap
library(pheatmap)
cM <- confusionMatrix(paste0(atacProj$Clusters), paste0(atacProj$Samples))
cM <- cM / Matrix :: rowSums(cM)

pHM <- pheatmap::pheatmap(mat = as.matrix(cM),
                          color = paletteContinuous("solarExtra"), 
                          border_color = "black")

plotPDF(pHM, name = "pHeatMap-for-Confucison Matrix.pdf",
        ArchRProj = atacProj, 
        addDOC = FALSE, 
        width = 7, 
        height = 7)
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Visualization(UMAP)-------------------------------------------------------------------
## UMAP Addition
atacProj <- addUMAP(ArchR = atacProj, reducedDims = "Iterative-LSI-1.0", name ="UMAP", 
                    nNeighbors  = 30, minDist = 0.5,
                    metric = "cosine")

plot <- plotEmbedding(ArchR = atacProj, colorBy = "cellcolData", 
                      name = c("Sample", "TSSEnrichment", "nFrags"),
                      embedding = "UMAP")
plot2 <- plotEmbedding(ArchRProj = atacProj, colorBy = "cellColData",
                       name = "Clusters",
                       embedding = "UMAP")

plotPDF(plot,
        name = "Plot-UMAP-Sample-TSS-nFrags.pdf", 
        ArchRProj = atacProj, 
        addDOC = FALSE, 
        width = 7, 
        height = 7
        )
plotPDF(plot2,
        name = "Plot-UMAP-Cells-colored-by-clusters.pdf", 
        ArchRProj = atacProj, 
        addDOC = FALSE, 
        width = 7, 
        height = 7
        )
        
## Dimensionality reduction after Harmony 
atacProj <- addUMAP(ArchRProj =  atacProj,
                    reducedDims =  "Harmony",
                    name = "UMAP-Harmony",
                    nNeighbors = 30,
                    minDist = 0.5, 
                    metric ="cosine")
                    
                    
## UMAP with embedding as Harmony
pl_H <- plotEmbedding(ArchRProj = atacProj,
                      colorBy = "cellColData",
                      name = c("Sample", "TSSEnrichment", "nFrags"),
                      embedding =  "UMAP-Harmony")
pl_h <- plotEmbedding(ArchRProj = atacProj,
                      colorBy = "cellColData",
                      name = "Clusters",
                      embedding = "UMAP-Harmony")
plotPDF(pl_H,
        name = "Plot-UMAP-Harmony-Sample-TSS-nFrags.pdf", 
        ArchRProj = atacProj, 
        addDOC = FALSE, 
        width = 7, 
        height = 7)
        
plotPDF(pl_h,
        name = "Plot-UMAP-Harmony-Clusters-coloredbycluster.pdf", 
        ArchRProj = atacProj, 
        addDOC = FALSE, 
        width = 7, 
        height = 7)
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Pseudo bulk replication----------------------------------------------

atacProj <- addGroupCoverages( ArchRProj = atacProj,
                               groupBy = "Clusters")
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Peak Calling with Tile Matrix----------------------------------------------------------
addArchRThreads(threads = 1)

## Calling peak with Tile Matrix
atacProj <- addReproduciblePeakSet(ArchRProj = atacProj,
                                   groupBy = "Clusters",
                                   peakMethod = "Tiles",
                                   method =  "p")

## Examine the merged peak set 
getPeakSet(atacProj)
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#Adding peak matrix----------------------------
atacProj <- addPeakMatrix(atacProj)

getAvailableMatrices(atacProj)
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Identifying Marker Peaks-------------------------------------
table(atacProj$Clusters)

mrkpeak <- getMarkerFeatures(ArchRProj =  atacProj, 
                             useMatrix = "PeakMatrix", 
                             groupBy = "Clusters",
                             bias = c("TSSEnrichment", "nFrags"),
                             testMethod = "wilcoxon")
mrkpeak

## GRangesList 
mrkList <- getMarkers(mrkpeak, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Plotting Marker Peak--------------------------------------
## Heat map 
HM_Peaks <- plotMarkerHeatmap(seMarker = mrkpeak,
                          cutOff = "FDR <= 0.01 & Log2FC >= 1",
                          transpose = TRUE)

plotPDF(HM_Peaks, name = "Peak-Marker-HeatMap.pdf", width = 8 , height=7, 
        ArchRProj = atacProj, addDOC = FALSE)


#


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Browser Track Containing Marker peaks-------------------
pl_bt <- plotBrowserTrack(ArchRProj =  atacProj,
                           groupBy = "Clusters",
                           geneSymbol = c("CD19", "CD34", "CD3D", "CD69"),
                           features =  getMarkers(mrkpeak, 
                                                  cutOff = "FDR <= 0.01 & Log2FC >= 1",
                                                  returnGR = TRUE),
                           upstream =  50000,
                           downstream =  50000) 
grid::grid.newpage()
grid::grid.draw(pl_bt$CD3D) # $"CD19", "CD34", "CD3D", "CD69"

plotPDF(pl_bt, name =  "Plot-Marker-with-Feature.pdf",
        ArchRProj = atacProj,
        width = 7,
        height = 6,
        addDOC = FALSE)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Gene score and Marker Features--------------------------------------
Gmrk <- getMarkerFeatures(ArchRProj = atacProj,
                          useMatrix = "GeneScoreMatrix",
                          groupBy =  "Clusters",
                          bias =  c("TSSEnrichment", "nFrags"),
                          testMethod =  "wilcoxon")

## Listing the markers for each cluster complex
mrkls <- getMarkers(Gmrk, cutOff =  "FDR<=0.01 & Log2FC >= 1") 
mrkls

## Adding impute weights
atacProj <- addImputeWeights(ArchRProj =  atacProj,
                               reducedDims = c("Iterative-LSI-1.0"),
                               dimsToUse = 1:30,
                               )

## Marker gene set
markerGenes <- c( "CD34", 
                  "GATA1", 
                  "PAX5", "MS4A1", "MME", "CD19",
                  "CD3D", "CD8A", "CD4", "TBX21",  
                  "CD14", "MPO", "CEBPB", 
                  "TBX21",
                  "CD38", "CD52", "CSF3R", 
                  "ANXA1","APLP2",
                  "CD19", "CD34" 
                  )                                                          
mrkgen_plot <- plotEmbedding(ArchRProj = atacProj, 
                             colorBy =  "GeneScoreMatrix",
                             name = markerGenes,
                             embedding = "UMAP",
                             quantCut = c(0.01,0.95),
                             imputeWeights =  NULL)

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Visualization of marker genes on Embedding _ cow plot(without MAGIC)--------------------------------------------------- 
plotPDF(poltList= cw_plot,
        name =  "Plot_Marker-Genes-WO-Imputation.pdf",
        ArchRProj = atacProj,
        addDOC = FALSE,
        width= 7,
        height = 6)

## Visual Marker gene set in form of heat map 
heatmapGS <- plotMarkerHeatmap( seMarker = Gmrk,
                            cutOff = "FDR<=0.01 & Log2FC >= 1",
                            labelMarkers = c( "CD34","GATA1","PAX5", "MS4A1", "MME", 
                                              "CD19", "CD3D", "CD8A", "CD4", "TBX21", 
                                              "CD14", "MPO","CEBPB","TBX21","CD38", 
                                              "CD52", "CSF3R","ANXA1","APLP2"), 
                            transpose = TRUE)


plotPDF(heatmapGS, name = "Gene-Score-Marker-Heatmap.pdf", ArchRProj =  atacProj, 
        addDOC = FALSE,
        width = 8,
        height = 7)

## Complex heat map 
png(filename = "Gene-Score-Marker-Heatmap", width =  1200, height = 800 )
ComplexHeatmap:: draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()
plotPDF(heatmapGS, name = "Gene-Score-Marker-Heatmap.pdf", ArchRProj =  atacProj, 
        addDOC = FALSE,
        width = 9,
        height = 8)
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Marker Gene imputation with MAGIC-------------------------------------
atacProj <- addImputeWeights(ArchRProj =  atacProj,
                             reducedDims = c("Iterative-LSI-1.0"),
                             dimsToUse = 1:30)

markerGenes <- c( "CD34", 
                  "GATA1", 
                  "PAX5", "MS4A1", "MME", "CD19",
                  "CD3D", "CD8A", "CD4", "TBX21",  
                  "CD14", "MPO", "CEBPB", 
                  "TBX21",
                  "CD38", "CD52", "CSF3R", 
                  "ANXA1","APLP2",
                  "CD19", "CD34" 
                  )
MKG_pl <- plotEmbedding(ArchRProj = atacProj, 
                        colorBy =  "GeneScoreMatrix",
                        name = markerGenes,
                        embedding = "UMAP",
                        quantCut = c(0.01,0.95),
                        imputeWeights =  getImputeWeights(atacProj))



plotPDF(plotList = MKG_pl, 
        name= "UMAP-Plot-of-Marker-genes-with Imputation.pdf",
        ArchRProj = atacProj,
        addDOC = FALSE,
        height = 8,
        width = 8)
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Tracking plot with ArchR Browser-----------------------------
pl_BrT <- plotBrowserTrack(ArchRProj = atacProj,
                           groupBy = "Clusters",
                           geneSymbol = markerGenes,
                           upstream = 50000,
                           downstream = 50000)
plotPDF(plotList = pl_BrT,
        name =  "Plot-BrowserTrack-Marker-Gene.pdf",
        ArchRProj = atacProj,
        addDOC =  FALSE,
        width = 9,
        height = 8)
        
## to plot specific gene 
grid::grid.newpage()
grid::grid.draw(pl_BrT$PAX5)
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Launching ArchR Browser--------------------------------------------------
ArchRBrowser(ArchRProj = atacProj)
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Integration with gene expression--------------------------------------

## Read the RDS file from project_1
scrna_annotation <- readRDS("/scrna_annotated_data.rds")
scrna_annotation

## Integration
atacProj <- addGeneIntegrationMatrix(ArchRProj = atacProj,
                                       useMatrix = "GeneScoreMatrix",
                                       matrixName = "GeneIntegrationMatrix",
                                       reducedDims = "Iterative-LSI-1.0",
                                       scrna_annotation = scrna_annotation,
                                       addToArrow = TRUE,
                                       force = TRUE,
                                       groupRNA = "hpca.main",
                                       nameCell = "predictedCell",
                                       nameGroup = "predictedGroup",
                                       nameScore = "predictedScore")

atacProj <- addImputeWeights(ArchRProj =  atacProj,
                             reducedDims = "Iterative-LSI-1.0")

markerHGenes <- c("CD34","GATA1","CD19","CD3D","CD14","CD38","TBX21")

pl_int <- plotEmbedding(ArchRProj = atacProj,
                        colorBy = "GeneScoreMatrix",
                        continuousSet = "horizonExtra",
                        name = markerHGenes,
                        embedding = "UMAP",
                        imputeWeights = getImputeWeights(atacProj))
plotPDF(plotList = pl_int ,
        name = "Plot-UMAP-Highlighted-gene_integeration-with-imputation.pdf",
        ArchRProj = atacProj,
        addDOC = FALSE,
        height = 6, width = 6)



palat <- paletteDiscrete(values = scrna_annotation$hpca.main)

## Linking gene expression data with unconstarined gene integration
atacProj.2 <- addGeneIntegrationMatrix(ArchRProj = atacProj,
                                       useMatrix = "GeneScoreMatrix",
                                       matrixName = "GeneIntegrationMatrix",
                                       reducedDims = "Iterative-LSI-1.0",
                                       scrna_annotation = scrna_annotation,
                                       addToArrow = FALSE,
                                       force = TRUE,
                                       groupRNA = "hpca.main",
                                       nameCell = "predictedCell_un",
                                       nameGroup = "predictedGroup_un",
                                       nameScore = "predictedScore_un"
                                       )

getAvailableMatrices(atacProj.2)

pl_un <- plotEmbedding(atacProj.2, colorBy = "cellColData",
                       name = "predictedGroup_un",
                       pal =  palat,
                       )

plotPDF( pl_un,
        name = "Plot-Unconstrained.pdf",
        ArchRProj = atacProj.2,
        addDOC = FALSE,
        height = 6, width = 6)

atacProj.2 <- addImputeWeights(ArchRProj =  atacProj.2,
                               reducedDims = c("Iterative-LSI-1.0"),
                               dimsToUse = 1:30)

pl_s1 <- plotEmbedding(ArchRProj = atacProj.2,
                       colorBy = "GeneIntegrationMatrix",
                       name = markerGenes,
                       continuousSet = "horizonExtra",
                       embedding = "UMAP",
                       imputeWeights = getImputeWeights(atacProj.2)
                       )

#pl_s2 <- plotEmbedding(ArchRProj = atacProj.2,
#                       colorBy = "GeneScoreMatrix",
#                       name = markerGenes,
#                       continuousSet = "horizonExtra",
#                       embedding = "UMAP",
#                       imputeWeights = getImputeWeights(atacProj.2)
#                      )

plotPDF(plotList = pl_s1,
        name = "Plot-UMAP-Marker-Genes-RNA-with-Imputation.pdf",
        ArchRProj = atacProj.2,
        addDOC = FALSE,
        height = 6, width = 6)

## Labelling sscATAC-seq with Clusters with scRNA-seq information
cn_atac_rna <- confusionMatrix(atacProj$Clusters, atacProj$predictedGroup)
labelOLD <- rownames(cn_atac_rna)
labelOLD

## Indetify cell type from predictedGroups
labelNEW <- colnames(cn_atac_rna)[apply(cn_atac_rna, 1 , which.max)]
labelNEW

## reclassify new clusters labels
remapC <- c("CD34" = "Early_Progenitor", 
            "GATA1" = "Erythroid", 
            "PAX5" = "B_cell", "MS4A1" = "B_cell", "MME" = "B_cell", "CD19" = "B_cell",
            "CD3D" = "T_cells", "CD8A" = "T_cells", "CD4" = "T_cells", "TBX21" = "T_cells", 
            "ILR7" = "T_cells",  
            "CD14" = "Monocyte", "MPO" = "Monocyte", "CEBPB" = "Monocyte", 
            "TBX21" = "NK",
            "CD38" = "LMPP", "CD52" = "LMPP", "CSF3R" = "LMPP", 
            "ANXA1" = "GMP","APLP2" = "GMP","AP351" = "GMP",
            "CD19" = "PRO-B", "CD34" = "PRO-B" )

remapC <- remapC[names(remapC) %in% labelNEW]

## New labele Mapping 
labelNEW2 <- mapLabels(labelNEW, oldLabels = names(remapC),
                       newLabels = remapC)
labelNEW2

## Combining old and new labels
atacProj$Clusters <- mapLabels(atacProj$Clusters,
                                 newLabels = labelNEW2,
                                 oldLabels = labelOLD)
## Plotting UMAP with new labels 
pl_labN <- plotEmbedding(atacProj, colorBy = "cellColData",
                         name = "Clusters")
plotPDF(pl_labN, name =  "Plot-UMAP-Remaped-Clusters.pdf",
        addDOC = FALSE, 
        width = 8, height = 6)
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#  Correlation Coefficient--------------------------------
segroupmotif <- getGroupSE(ArchRProj = atacProj,
                           useMatrix = "MotifMatrix",
                           groupBy = "Clusters")
segroupmotif

sez <- segroupmotif[rowData(segroupmotif)$seqnames=="z",]

rowData(sez)$maxDelta <- lapply(seq_len(ncol(sez)),function(x){
  rowMaxs(assay(sez)-assay(sez)[,x])
}
  ) %>% Reduce("cbind",.) %>% rowMaxs

corgsm_mm <- correlateMatrices(ArchRProj = atacProj,
                               useMatrix1 = "GeneScoreMatrix",
                               useMatrix2 = "MotifMatrix",
                               reducedDims = "Iterative-LSI-1.0"
                               ) 
corgsm_mm

corgsm_mm <- corgsm_mm[order(abs(corgsm_mm$cor),decreasing = TRUE),]
corgsm_mm <- corgsm_mm[which(!duplicated(gsub("\\-.*","",corgsm_mm[,"MotifMatrix_name"]))),]
corgsm_mm$TFRegulator <- "NO"
corgsm_mm$TFRegulator[which(corgsm_mm$cor > 0.5 & corgsm_mm$padj < 0.01 &  
                              corgsm_mm$maxDelta > quantile(corgsm_mm$maxDelta, 
                                                            0.75))] <- "YES"
sort(corgsm_mm[corgsm_mm$TFRegulator == "YES",1])

Corr  <- data.frame(corgsm_mm)
write_xlsx(Corr , " Correlated Matrix.xlsx")




pl_cor <- ggplot(data.frame(corgsm_mm), aes(cor,maxDelta, color = TFRegulator)) +
  geom_point() + theme_ArchR(color = "black")+ 
  geom_vline(xintercept = 0, lty = "dashed")+
  scale_color_manual(values = c("NO" = "darkgrey", "YES" = "firebrick3"))+
  xlab("Correlation To Gene Score") + ylab("Max TF Motif Delta") + 
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,max(corgsm_mm$maxDelta)*1.05)
                    )





#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Transcription Factor and Motif Activity------------------------------
## Motif Enrichment in Differential Peaks===
atacProj <- addMotifAnnotations(ArchRProj = atacProj,
                                motifSet = "cisbp",
                                name = "Motif")
## Peak Annotation Enrichment
motifsUP <- peakAnnoEnrichment( seMarker = mrkpeak, 
                                ArchRProj = atacProj,
                                peakAnnotation = "Motif",
                                cutOff = "FDR<=0.01 & Log2FC >= 1")
## Data Frame 
df <- data.frame(TF = rownames(motifsUP),
                 mlog10Padj = assay(motifsUP)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

head(df)
write_xlsx(df, "DataFrame to top Motifs.xlsx")


## Plotting gg plot 
motifggUP <- ggplot(df , aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) + 
  ggrepel::geom_label_repel(data = df[rev(seq_len(30)),], aes(x =rank, y= mlog10Padj,label =TF),
                            size = 1.5,
                            nudge_x = 2,
                            color = "black") + theme_ArchR()+
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colours = paletteContinuous(set = "comet"))

plotPDF(motifggUP,name = "GG-PLOT-For-Ranked-Sorted-TFs-Motif.pdf",
        ArchRProj = atacProj,
        addDOC = FALSE,
        width = 7, height = 6)

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Motif Deviations--------------------------------------------
if("Motif" %ni% names(atacProj@peakAnnotation)){
  atacProj <- addMotifAnnotations(ArchRProj = atacProj,
                                  motifSet = "cisbp",
                                  name= "Motif")
  }
## Adding Background peaks
atacProj <- addBgdPeaks(atacProj)
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Adding Deviation Matrix-----------
atacProj <- addDeviationsMatrix(ArchRProj = atacProj,
                                peakAnnotation = "Motif",
                                force = TRUE)
                                
# getVatDeviation
plotVarDev <- getVarDeviations(atacProj, name = "MotifMatrix", plot =  TRUE)
plotVarDev

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores.pdf",
        ArchRProj = atacProj, addDOC =  FALSE,
        width = 6, height =  6)

motifs <- c("BCL11A_194", "BCL11B_825")
mrkMotif <- getFeatures(atacProj, select = paste(motifs, collapse = "|"),
                        useMatrix = "MotifMatrix")


mrkMotif <- grep("z:", mrkMotif, value = TRUE)
mrkMotif

pl_motif <- plotGroups(ArchRProj =  atacProj,
                       groupBy = "Clusters",
                       colorBy = "MotifMatrix",
                       name = mrkMotif,
                       imputeWeights = getImputeWeights(atacProj))

plotPDF(pl_motif, name = "Plot-Groups-Devation-with-Imputation.pdf",
        ArchRProj = atacProj, addDOC =  FALSE,
        width = 6 , height = 5)

## Plot Embedding-UMAP
pl_M0tif <- plotEmbedding(ArchRProj = atacProj,
                          colorBy = "MotifMatrix",
                          name = sort(mrkMotif),
                          embedding = "UMAP",
                          imputeWeights = getImputeWeights(atacProj))
plotPDF(pl_M0tif, name = "Plot-UMAP-MotifMarker-with-Imputation.pdf",
        ArchRProj = atacProj, addDOC =  FALSE,
        width = 6 , height = 5)

                                 
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#Differentail accesibility------------------------
table(scrna_annotation$hpca.main)
table(atacProj$Clusters)
table(atacProj$predictedGroup)
markertest <- getMarkerFeatures(ArchRProj =  atacProj,
                                useMatrix = "PeakMatrix",
                                groupBy =  "Clusters",
                                testMethod = "wilcoxon",
                                bias = c("TSSEnrichment", "nFrags"),
                                useGroups = "T_cells",
                                bgdGroups = "Erythroblast"
                                )

    
# Marker peak MA & Volcano plot   

pma <- plotMarkers(seMarker = markertest, 
                   name = "T_cells",
                   cutOff = "FDR <= 0.01 & Log2FC >= 1", 
                   plotAs = "MA",
                   scaleTo = 10^4) 
pv <- plotMarkers(seMarker = markertest, 
                  name = "T_cells",
                  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
                  plotAs = "Volcano",
                  scaleTo = 10^4)

plotPDF(pma,pv, 
        name = "Marker-MA-Volcano-plot.pdf", width = 8 , height=7, 
        ArchRProj = atacProj, addDOC = FALSE)

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Foot-printing TFs-------------------------------

## to find the relevant position the motifs 
motifPosition <- getPositions(ArchRProj = atacProj, name = "Motif") 
motifPosition

## subsetting for specific motif of interest

markerMotif <- unlist(lapply(motifs, function(x)grep(x, names(motifPosition),
                                                       value = TRUE)
  
                             )
                      )
markerMotif


## Get FootPrint 
se_foot <- getFootprints(ArchRProj = atacProj,
                         positions = motifPosition[markerMotif],
                         groupBy = "Clusters")


## Normalization of Footprints for Tn5 Bias

### Subtracting the Tn5 Bias
plotFootprints(seFoot = se_foot, 
               ArchRProj = atacProj,
               normMethod = "Subtract",
               plotName = "Footprints-subtract-Bias",
               addDOC = FALSE,
               smoothWindow = 5)

### Dividing by the Tn5 Bias
plotFootprints(seFoot = se_foot,
               ArchRProj = atacProj,
               normMethod = "Divide",
               plotName = "Footprints-Divide-Bias",
               addDOC = FALSE,
               smoothWindow = 5)

## Feature Foot-printing
# TSS insertion profile

se_tss <- getFootprints(ArchRProj = atacProj,
                        positions = GRangesList(TSS = getTSS(atacProj)),
                        groupBy = "Clusters",
                        flank = 2000
                        )

# TSS insertion for each cell group
plotFootprints(seFoot = se_tss,
               ArchRProj = atacProj,
               normMethod = "None",
               plotName = "TSS-No-Normalization",
               addDOC = FALSE,
               flank = 2000,
               flankNorm = 100
               )
# Plotting Aggregate Footprintg for Tcell and Monocyte
table(scrna_annotation$hpca.main)
tcells <- scrna_annotation$hpca.main

getFootprints(ArchRProj = atacProj,
              positions =  "T_cells")


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
# Co-accessibility-------------------------------------

## Co-accessibility of peaks
###  Adding co-accessible peaks
atacProj_coap <- addCoAccessibility(ArchRProj = atacProj, 
                                  reducedDims = "Iterative-LSI-1.0")
### Retrieving the co-accessible peaks
coap_rt <- getCoAccessibility(ArchRProj = atacProj_coap,
                              corCutOff = 0.5,
                              resolution = 1,
                              returnLoops = FALSE)

coap_rt

metadata(coap_rt)[[1]]

## returnloop = TRUE
coap_rt <- getCoAccessibility(ArchRProj = atacProj_coap,
                              corCutOff = 0.5,
                              resolution = 1000,
                              returnLoops = TRUE
                              )

coap_rt[[1]]

bt_coap <- plotBrowserTrack(ArchRProj = atacProj_coap,
                            groupBy = "Clusters",
                            geneSymbol = markerGenes,
                            upstream = 50000,
                            downstream = 50000,
                            loops = getCoAccessibility(atacProj_coap))
plotPDF(plotList = bt_coap,
        name = "Plot-Track-Marker-Genes-with-CoAccessibility.pdf",
        ArchRProj = atacProj_coap,
        addDOC = FALSE,
        width = 6, height = 6)
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#Peak To gene Linkage-----------------------------------------------
atacProj_p2g <- addPeak2GeneLinks(ArchRProj = atacProj_coap,
                                  reducedDims = "Iterative-LSI-1.0")

p2g <- getPeak2GeneLinks(ArchRProj = atacProj_p2g,
                         corCutOff = 0.5,
                         resolution = 1000,
                         returnLoops = TRUE)
metadata(p2g)[[1]]
p2g[[1]]
markerGenes2 <- c("CD34","GATA1","CD19","CD3D", "CD14", "TBX21", "CD38", "AP351")

bt_p2g <-  plotBrowserTrack(ArchRProj = atacProj_p2g,
                            groupBy = "Clusters",
                            geneSymbol = markerGenes2,
                            upstream = 50000,
                            downstream = 50000,
                            loops = getPeak2GeneLinks(atacProj_p2g)
                            )

plotPDF(plotList = bt_p2g,
        name = "Plot-Track-Marker-Genes-with-Peak2GeneLinks.pdf",
        ArchRProj = atacProj_p2g,
        addDOC = FALSE,
        width = 6, height = 6)
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
