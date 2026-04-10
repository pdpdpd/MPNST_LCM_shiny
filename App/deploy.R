# Deploy MPNST LCM Shiny app to Posit Connect Cloud
#
# Run this script from the App/ directory before pushing to git,
# or whenever you add/remove packages or data files.
#
# Usage:
#   setwd("App/")
#   source("deploy.R")

# Install rsconnect if needed
if (!requireNamespace("rsconnect", quietly = TRUE)) {
  install.packages("rsconnect")
}

library(rsconnect)

# ------------------------------------------------------------------
# Step 1: Regenerate manifest.json
#   This captures your current R version and installed package versions.
#   Commit the updated manifest.json alongside any code changes.
# ------------------------------------------------------------------
rsconnect::writeManifest(
  appDir    = ".",
  appFiles  = c(
    "app.R",
    "LCM_coordinates.rds",
    "LCM_CN_available.rds",
    "LCM_CN_mtx.rds",
    "LCM_chr_probes.rds",
    "LCM_haplotype_info.rds",
    "study_design.png",
    "Primary_12-24POSITIONS_WRITTEN.JPG",
    "Primary_B_OV20X_cut_positions.jpg",
    "Primary_C_side_OV20x_CUT_positions.jpg",
    "T1.1_1.25X_CUT_24.JPG",
    "T1.1_B_Overview_cut_positions.jpg",
    "T1.1_side_OV_20x_CUT_positions.jpg",
    "T4.1_Front_OVERVIEW_CUT_6.3x_positions.jpeg",
    "T4.1_Back_OV_20x_cut.jpg",
    "T4.1_side_OV_cut_20x.jpg"
  )
)

message("manifest.json updated. Commit it and push to trigger a Connect Cloud redeploy.")

# ------------------------------------------------------------------
# Step 2 (optional): Push directly via rsconnect
#
# Uncomment and fill in your Connect Cloud server URL and API key
# to deploy programmatically instead of using git-backed deployment.
#
# rsconnect::addServer(
#   url  = "https://connect.posit.cloud",
#   name = "posit-connect-cloud"
# )
#
# rsconnect::connectApiUser(
#   account = "<your-account-name>",
#   server  = "posit-connect-cloud",
#   apiKey  = Sys.getenv("CONNECT_API_KEY")   # set in .Renviron or environment
# )
#
# rsconnect::deployApp(
#   appDir    = ".",
#   appName   = "MPNST_LCM_shiny",
#   server    = "posit-connect-cloud",
#   account   = "<your-account-name>",
#   forceUpdate = TRUE
# )
# ------------------------------------------------------------------
