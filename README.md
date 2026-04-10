# MPNST LCM Explorer

Interactive Shiny app for exploring spatial copy number heterogeneity in a malignant peripheral nerve sheath tumour (MPNST) through laser capture microdissection (LCM).

## App features

- **LCM Spot tab** — Select a tumour region and tissue section orientation (Front/Back/Side), click individual LCM spots to view their allele-specific copy number profile and haplotype balance.
- **LCM Region tab** — View copy number profiles of all successfully sequenced spots in a region; click a genomic segment to see its copy number state distributed spatially across all three section orientations.

## Publishing to Posit Connect Cloud (git-backed deployment)

### First-time setup

1. **Generate a fresh `manifest.json`** using your local R installation (R ≥ 4.3 recommended):

   ```r
   setwd("App/")
   source("deploy.R")   # runs rsconnect::writeManifest()
   ```

2. **Commit the updated manifest** and push to the repo:

   ```bash
   git add App/manifest.json
   git commit -m "Update manifest.json for Connect Cloud deployment"
   git push
   ```

3. **Link the repo to Posit Connect Cloud**:
   - Sign in to [connect.posit.cloud](https://connect.posit.cloud)
   - Click **Publish → Import from Git**
   - Enter the repository URL and select the branch
   - Connect will detect the `App/` directory (it contains `manifest.json`)
   - Set **Primary file** to `app.R` and click **Deploy**

### Updating the app

After changing `app.R` or any data files, re-run `deploy.R` to refresh the manifest, then commit and push. Connect Cloud will automatically redeploy on the next push (or you can trigger it manually from the **Settings → Info** panel).

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
