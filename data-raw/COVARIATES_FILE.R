#' List of covariates managed by Christopher Jin for the statistical analysis
config = jsonlite::fromJSON("~/config.json")
repo_local_dir = config$precovid_repo_path
COVARIATES_FILE = read.csv(file = file.path(repo_local_dir, "library", "covariates_pre_cawg.csv"))

usethis::use_data(COVARIATES_FILE, overwrite = TRUE)
