#' @title Check package installation
#'
#' @description Several secondary R packages are used and it's preferable to
#'   install them only if they are needed. Placing this check in the body of a
#'   function will provide helpful instructions for cases when a user attempts
#'   to run a function but doesn't have the required suggested package(s).
#'   (The main inspiration of this was because a lot of people had issues installing
#'   some packages needed for only specific analysis)
#'
#' @param pkg (*character*) the name of the suggested R package
#' @param fun (*character*) the name of the function that uses `pkg`
#' @param task (*character*; *optional*) the task that `fun` completes.
#'
#' @return `TRUE`, invisibly, if `pkg` is installed.
#'
#' @examples
#' \dontrun{
#' result <- MotrpacHumanPreSuspensionAnalysis::check_package_installation(pkg = 'base', fun = 'mean')
#' result
#' }

check_package_installation <- function(pkg, fun = NULL, task = NULL) {

  pkg_installed <- requireNamespace(pkg, quietly = TRUE)

  # check success - end here and return TRUE, invisibly
  if (pkg_installed) {
    library(pkg, character.only = TRUE)

    return(invisible(TRUE))
  }
  if(is.null(fun)) fun = "function unspecified"
  # chuck unsuccessful - give instructions for fix
  msg <- paste0(
    "package \"", pkg,"\" must be installed to use ", fun, "()"
  )

  if (!is.null(task)) {
    msg <- paste(msg, "for", task)
  }

  if(pkg == "MotrpacHumanPreSuspensionData"){
    msg = "This function requires the MotrpacHumanPreSuspensionData package, which is a private repository that is available upon request via the motrpac-data.org"
  }

  stop(msg, call. = FALSE)
}
