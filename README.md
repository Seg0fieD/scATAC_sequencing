# scATAC seq (Single Cell- Assay for Transposase-Accessible Chromatin) Analysis 

This project focused on analyzing single-cell chromatin accessibility data from human brain cortex cells, including progenitors and differentiated cells, using the ArchR package in R.

## Key Tasks:

1. **Preprocessing and Quality Control**:
   - The data was downloaded and filtered based on the number of fragments and TSS enrichment. Doublets were identified and removed.
   - Various quality control measures were applied, including fragment length distribution, TSS enrichment, and the number of fragments per sample.

2. **Dimensionality Reduction**:
   - Dimensionality reduction was performed using ArchR’s iterative latent semantic indexing (LSI) method.
   - UMAP plots were generated to visualize the data, with cells colored by sample, TSS enrichment, and fragment number.
   - Batch effects were identified and corrected, and the impact of batch correction was visualized.

3. **Clustering**:
   - Louvain clustering was applied to all cells, and the resulting clusters were visualized using UMAP. The sample proportions within each cluster were also analyzed.

4. **Peaks**:
   - A joint peak set was computed for the dataset, and cluster-specific marker peaks were identified. 
   - Heatmaps were created to show accessibility in these marker peaks, and specific genes such as TOP2A, MKI67, AURKA, SATB2, and VGLUT1 were analyzed.

5. **Gene Activity**:
   - Gene activity scores were computed using chromatin accessibility, and cluster marker genes were identified based on specific parameters.
   - MAGIC was applied to visualize the first five marker genes with and without smoothing.

6. **Transcription Factor (TF) Motif Activity**:
   - TF motif activity was computed using an appropriate annotation database. 
   - UMAP embeddings were generated for the top TF motifs, and their activity across different clusters was visualized.

7. **Integration with Gene Expression**:
   - The scRNA-seq data was integrated with scATAC-seq data to link gene expression with chromatin accessibility. UMAP plots were generated to visualize this integration for specific marker genes.
   - Correlation coefficients between gene expression and activity were computed.

8. **Peak-Gene Linkage**:
   - Peak-gene linkage was computed using both gene expression and chromatin accessibility data. 
   - Heatmaps were created to visualize the linkage across all cell types.

9. **Differential Accessibility**:
   - Differential peak accessibility between GluN5 and Cyc. Prog. cells was computed, and MA and volcano plots were created.
   - TF motif enrichment was analyzed for differentially accessible peaks.

10. **TF Footprinting**:
   - TF footprints were obtained for the top three motifs from the GluN5 and Cyclic progenitor cells.
   - Normalization for Tn5 bias was performed, and aggregate footprints were visualized and interpreted.

11. **Co-accessibility**:
   - Co-accessibility of peaks was computed and visualized using genome track plots. 
   - Potential enhancers for marker genes were identified by analyzing peak linkage to gene TSS.

## R Script:
The R script used for analysis is named **scATAC_seq_script**.

## Results:
The resulting plots and images can be found in the `results/img` directory of this repository.

## ArchR Documentation:
For further information on the ArchR package, please refer to the official [ArchR documentation](https://www.archrproject.com/index.html).
