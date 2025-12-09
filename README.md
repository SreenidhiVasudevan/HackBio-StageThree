# HackBio-StageThree
  # SARS-CoV-2 Infection Dynamics Single-Cell RNA-seq Analysis

This project investigates the cellular landscape changes in lung tissue during SARS-CoV-2 infection using single-cell RNA sequencing (scRNA-seq) data. The analysis tracks changes across mock (uninfected), 1 day post-infection (1 dpi), 2 dpi, and 3 dpi stages, focusing on cell type identification, gene expression patterns of key markers (ACE2, ENO2), and trajectory inference.

## How it was done and run

### 1. Environment Setup & Package Installation

All necessary Python libraries for single-cell analysis, including ⁠ scanpy ⁠, ⁠ anndata ⁠, ⁠ igraph ⁠, ⁠ celltypist ⁠, ⁠ decoupler ⁠, ⁠ fa2-modified ⁠, ⁠ louvain ⁠, and ⁠ scvelo ⁠, were installed using ⁠ pip ⁠ commands.

### 2. Data Acquisition and Loading

The scRNA-seq raw data (⁠ GSE166766_RAW.tar ⁠) was downloaded from the NCBI Gene Expression Omnibus (GEO) repository using ⁠ wget ⁠. The tar archive was then extracted, and individual 10x Genomics datasets for each infection stage (mock, 1dpi, 2dpi, 3dpi) were loaded into ⁠ AnnData ⁠ objects using ⁠ scanpy.read_10x_mtx ⁠. These ⁠ AnnData ⁠ objects were subsequently concatenated into a single ⁠ AnnData ⁠ object for combined analysis, with a 'condition' annotation added to distinguish the samples.

### 3. Quality Control (QC) & Preprocessing

*   *Gene Annotation:* Mitochondrial (MT), Ribosomal (RIBO), and Hemoglobin (HB) genes were identified for each dataset.
*   *QC Metrics Calculation:* ⁠ scanpy.pp.calculate_qc_metrics ⁠ was used to compute the number of genes by counts, total counts, and percentage of MT, RIBO, and HB gene counts.
*   *Visualization of QC Metrics:* Violin plots and scatter plots (⁠ total_counts ⁠ vs. ⁠ n_genes_by_counts ⁠ colored by ⁠ pct_counts_MT ⁠) were generated to visualize the quality distributions of cells in each sample.
*   *Normalization & Log-Transformation:* Data for each sample was normalized using ⁠ scanpy.pp.normalize_total ⁠ and log-transformed using ⁠ scanpy.pp.log1p ⁠.
*   *Feature Selection:* Highly variable genes (HVGs) were identified for each dataset using ⁠ scanpy.pp.highly_variable_genes ⁠ (top 1000 genes).

### 4. Dimensionality Reduction & Clustering

*   *Principal Component Analysis (PCA):* PCA was performed on each dataset to reduce dimensionality, and variance ratio plots were generated.
*   *Uniform Manifold Approximation and Projection (UMAP):* ⁠ scanpy.pp.neighbors ⁠ and ⁠ scanpy.tl.umap ⁠ were applied to compute neighborhood graphs and UMAP embeddings for visualization.
*   *Leiden Clustering:* The Leiden algorithm (⁠ scanpy.tl.leiden ⁠) was applied to cluster cells into distinct groups based on their similarity.

### 5. Cell Type Annotation

*   *Marker Gene Database:* The ⁠ decoupler ⁠ library was used to query the PanglaoDB marker gene database, filtered for human lung-specific cell types.
*   *Universal Lung Markers (ULM):* The ⁠ dc.mt.ulm ⁠ function was applied to assign cell type scores to individual cells based on the expression of marker genes in highly variable genes. Temporary ⁠ AnnData ⁠ objects were created for this step to use only HVGs.
*   *Differential Cell Type Expression:* ⁠ dc.tl.rankby_group ⁠ was used to identify top cell types enriched in each Leiden cluster, helping to assign a probable cell identity to each cluster.
*   *Mapping Annotations:* The most enriched cell type for each Leiden cluster was mapped back to the original ⁠ AnnData ⁠ objects as ⁠ annotated_cell_type ⁠.
*   *Visualization:* UMAP plots colored by ⁠ leiden_res_ ⁠ and specific cell type ULM scores (e.g., 'Ciliated cells', 'Clara cells') were generated to visualize cluster distribution and cell type enrichment.

### 6. Trajectory Inference

*   *Partition-based Graph Abstraction (PAGA):* ⁠ scanpy.tl.paga ⁠ was performed to infer the global topology of the single-cell data, showing relationships between cell clusters.
*   *Pseudotime Analysis (Diffusion Pseudotime - DPT):* ⁠ scanpy.tl.dpt ⁠ was run on each dataset, rooted at Leiden cluster '0', to infer developmental trajectories and pseudotime values for individual cells. UMAP plots colored by ⁠ dpt_pseudotime ⁠ and ⁠ leiden_res_ ⁠ were used for visualization.

## Observations from the Data

### Summarized Identified Cell Types per Stage:
*   *Mock (Uninfected) Stage:* Clara cells, Pulmonary alveolar type I cells, Airway goblet cells, Ionocytes.
*   *Day 1 Post-Infection:* Pulmonary alveolar type I cells, Ionocytes, Clara cells, Airway goblet cells, Pulmonary alveolar type II cells.
*   *Day 2 Post-Infection:* Pulmonary alveolar type II cells, Ciliated cells, Airway goblet cells, Ionocytes.
*   *Day 3 Post-Infection:* Ciliated cells, Ionocytes, Airway goblet cells, Alveolar macrophages.

### Correlation of Cell Types with COVID-19 Infection:
*   *Alveolar Epithelial Cells (Type I and Type II):* Crucial for gas exchange; AT2 cells express ACE2, a primary entry point for SARS-CoV-2. Damage impacts lung function.
*   *Airway Epithelial Cells (Ciliated, Goblet, Clara Cells, Ionocytes):* Line airways, involved in mucociliary clearance and protection. Disruption by viral infection impairs clearance, causes inflammation.
*   *Alveolar Macrophages:* Immune cells critical for pathogen clearance. Their activation and function are key in immune response, but dysregulation can lead to severe pathology.

### Evaluation of ACE2 as a COVID-19 Marker:
ACE2 is a strong marker for tracking COVID-19 infection. Its expression pattern showed a clear and dynamic change from low and scattered in the mock stage, to increasingly localized and intense in specific cell clusters by Day 3 post-infection, directly indicating its involvement in the infection process.

### Comparison of ENO2 and ACE2 as Biomarkers:
*   *ACE2* exhibits a dynamic and specific expression pattern, becoming progressively more localized and intense in specific cell clusters from Day 1 to Day 3 post-infection. This directly correlates with its role as the viral entry receptor, making it a strong biomarker for infection progression.
*   *ENO2* shows a broad and relatively stable expression across all stages, suggesting it acts as a general metabolic marker rather than a specific indicator of acute infection dynamics.

## Final Inference from the Results

At 3 days post-infection (3 dpi), *Ciliated cells (specifically Cluster 0)* exhibited the highest expression of ACE2. This is biologically significant as it points to ciliated cells as a primary site for viral entry and replication within the airways during the later stages of infection. This heavy involvement of ciliated cells likely contributes to compromised mucociliary clearance and localized inflammation, which are key pathological features observed in COVID-19. Understanding this cell-type-specific tropism is crucial for developing targeted therapies aimed at reducing viral load and restoring airway function. The dynamic changes in cell type composition and ACE2 expression over the infection course underscore the complex host response to SARS-CoV-2.
