# MPNST LCM Explorer
## Original code and design from Haixi Yan (haixi.yan@crick.ac.uk)

Interactive Shiny app for exploring spatial copy number heterogeneity in a malignant peripheral nerve sheath tumour (MPNST) through laser capture microdissection (LCM).
Now published at https://vanloo-lab-mpnst-lcm-shiny.share.connect.posit.cloud.
 
## App features

- **LCM Spot tab** — Select a tumour region and tissue section orientation (Front/Back/Side), click individual LCM spots to view their allele-specific copy number profile and haplotype balance.
- **LCM Region tab** — View copy number profiles of all successfully sequenced spots in a region; click a genomic segment to see its copy number state distributed spatially across all three section orientations.


## App directory structure

```
App/
├── app.R                          # Shiny application
├── manifest.json                  # Posit Connect Cloud deployment manifest
├── deploy.R                       # Helper script to regenerate manifest
├── LCM_coordinates.rds            # Spot spatial coordinates and metadata
├── LCM_CN_available.rds           # Samples with successful CN profiles
├── LCM_CN_mtx.rds                 # Allele-specific copy number matrix
├── LCM_chr_probes.rds             # Chromosome probe information
├── LCM_haplotype_info.rds         # Haplotype read counts per spot
├── study_design.png               # Study design overview image
└── *.JPG / *.jpg / *.jpeg         # Tissue section images (9 sections)
```

## R packages required

`shiny`, `shinydashboard`, `ggplot2`, `dplyr`, `magrittr`, `pheatmap`, `jpeg`, `grid`

