# Cytof integration, batch correction and downstream analysis

# Description 
This code was developed for integrating, batch correcting, and analyzing CyTOF data acquired through multiple batches. All FCS files are normalized together using the R package premessa version 0.2.4 with the four-element calibration beads (Fluidigm, cat: 201078). After normalization, live singlets were gated. Each marker's intensities were capped at the 1st and 99th percentiles, normalized from 0 to 1, and centered at the mean. Up to 20,000 cells were randomly subsampled from each sample. Dimensional reduction was performed using the Python implementation of UMAP via the reticulate R package version 1.18. Unsupervised clustering was carried out using the PhenoGraph or ClusterX algorithm with the cytofkit54 R package version 1.4.10. The median expression of each marker in each cluster was visualized using the pheatmap R package version 1.0.12, and immune cell populations were identified based on the expression of specific markers.

# Notes
Obtain and unzip the CyTOF_renamed.zip file to access the input files

This code was originally developed for the paper: 
Lee, A.H., Sun, L., Mochizuki, A.Y. et al. Neoadjuvant PD-1 blockade induces T cell and cDC1 activation but fails to overcome the immunosuppressive tumor associated macrophages in recurrent glioblastoma. Nat Commun 12, 6938 (2021)

The code available here are the adapted codes utilized for the paper:
Everson, R.G., Hugo, W., Sun, L. et al. TLR agonists polarize interferon responses in conjunction with dendritic cell vaccination in malignant glioma: a randomized phase II Trial. Nat Commun 15, 3882 (2024). 

Requirements: R (tested with version 4.2.1).
