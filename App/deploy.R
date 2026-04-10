# Deploy MPNST LCM Shiny app to Posit Connect Cloud
#
# Can be run from any working directory:
#   Rscript App/deploy.R
#   Rscript deploy.R        (from inside App/)
#   source("App/deploy.R")  (from R console)

# Locate the App/ directory regardless of where this script is called from
app_dir <- tryCatch({
  # Works when run via Rscript
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  normalizePath(dirname(sub("--file=", "", file_arg[1])))
}, error = function(e) {
  # Works when source()-ed interactively
  tryCatch(
    normalizePath(dirname(sys.frame(1)$ofile)),
    error = function(e2) getwd()
  )
})
message("App directory: ", app_dir)

# Use a local library directory if the system library is not writable
local_lib <- Sys.getenv("R_LIBS_USER", unset = file.path(Sys.getenv("HOME"), "R", "library"))
if (!dir.exists(local_lib)) dir.create(local_lib, recursive = TRUE)
if (!local_lib %in% .libPaths()) .libPaths(c(local_lib, .libPaths()))

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com"))

# Install rsconnect if needed
if (!requireNamespace("rsconnect", quietly = TRUE)) {
  install.packages("rsconnect", lib = local_lib)
}

library(rsconnect)

# ------------------------------------------------------------------
# Step 1: Regenerate manifest.json
#   This captures your current R version and installed package versions.
#   Commit the updated manifest.json alongside any code changes.
# ------------------------------------------------------------------
rsconnect::writeManifest(appDir = app_dir)

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
#   appDir    = app_dir,
#   appName   = "MPNST_LCM_shiny",
#   server    = "posit-connect-cloud",
#   account   = "<your-account-name>",
#   forceUpdate = TRUE
# )
# ------------------------------------------------------------------
