
# Modification of roxygen2 block_set_env internals
.custom_block_set_env <- function(block, env) {
  block_evaluate <- getFromNamespace("block_evaluate", "roxygen2")
  block_find_object <- getFromNamespace("block_find_object", "roxygen2")

  block <- block_evaluate(block, env)
  block <- block_find_object(block, env)

  val <- block[["object"]][["value"]]

  if (any(class(val) == "motrdat")) {
    block[["object"]][["value"]][["data"]] <- data.frame()
  } else if (is.data.frame(val)) {
    block[["object"]][["value"]] <- data.frame()
  }

  return(block)
}

# Identical to roxygen2 block_set_env internals
.block_set_env <- function(block, env) {
  block_evaluate <- getFromNamespace("block_evaluate", "roxygen2")
  block_find_object <- getFromNamespace("block_find_object", "roxygen2")

  block <- block_evaluate(block, env)
  block <- block_find_object(block, env)

  return(block)
}

# https://r-pkgs.org/code.html#sec-code-onLoad-onAttach
.onLoad <- function(libname, pkgname) {
  # Only apply the roxygen2 hack during development (when roxygen2 is loaded)
  if (!isNamespaceLoaded("roxygen2")) return(invisible())

  # Hack so that devtools::document() doesn't take hours
  environment(.custom_block_set_env) <- asNamespace("roxygen2")

  suppressWarnings({
    utils::assignInNamespace(x = "block_set_env",
                             value = .custom_block_set_env,
                             ns = "roxygen2")
  })

  invisible()
}

.onUnload <- function(libname, pkgname) {
  if (!isNamespaceLoaded("roxygen2")) return(invisible())

  environment(.block_set_env) <- asNamespace("roxygen2")

  suppressWarnings({
    utils::assignInNamespace(x = "block_set_env",
                             value = .block_set_env,
                             ns = "roxygen2")
  })

  invisible()
}
