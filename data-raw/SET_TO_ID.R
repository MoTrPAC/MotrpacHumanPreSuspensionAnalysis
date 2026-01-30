library(dplyr)
library(MotrpacHumanPreSuspensionData)

## Helper functions ----

# Split one or more strings into components separated by `split`. Then, shorten
# the string by selecting the first components that do not exceed a certain
# number of total characters, `n`.
.cutstr <- function(x, split = "", n = Inf) {
  x <- strsplit(x, split = split)
  x <- vapply(x, function(xi) {
    keep <- cumsum(nchar(xi)) + nchar(split) * (seq_along(xi) - 1L) <= n
    xi <- paste(xi[keep], collapse = split)

    return(xi)
  }, character(1L))

  return(x)
}

.cutstr <- Vectorize(.cutstr)

# Number of digits used for each ID is a function of the number of sets
n_digits <- floor(log10(sum(lengths(MOLECULAR_SIGNATURES)))) + 1L

SET_TO_ID <- data.frame(
  database = rep(names(MOLECULAR_SIGNATURES),
                 lengths(MOLECULAR_SIGNATURES)),
  set = unlist(lapply(MOLECULAR_SIGNATURES, names))
) %>%
  mutate(
    # Split PTMSigDB into categories
    database = ifelse(grepl("^PTMSIGDB", database),
                      sub("(^PTMSIGDB_[^_]+).*", "\\1", set),
                      database),
    database = factor(database,
                      levels = unique(database))
  ) %>%
  # Make sure PTMSigDB sets are last, followed by CellMarker gene sets, to avoid changing the set_id
  arrange(grepl("^CELLMARKER", database),
          grepl("^PTMSIGDB", database),
          set) %>%
  mutate(set_id = sprintf(paste0("%0", n_digits, "d"), seq_len(n())),
         # Remove database from set
         set_short = ifelse(database %in% c("MITOCARTA", "PSP",
                                            "REFMET", "CELLMARKER"),
                            sub("^[^_]+_", "", set),
                            set),
         # All CellMarker gene sets are human, so remove that keyword
         set_short = ifelse(database == "CELLMARKER",
                            sub(" Human$", "", set_short),
                            set_short),
         # Cut set_short between words so it does not exceed 40 chars
         temp = .cutstr(set_short,
                        split = ifelse(
                          database %in% c("MITOCARTA", "PSP",
                                          "REFMET", "CELLMARKER"),
                          " ", "_"), n = 50L),
         # If the description without the database is longer than the
         # shortened version with '...(uniqueID)' included (plus a few
         # characters as a buffer), use the latter; otherwise, use the
         # former.
         set_short = ifelse(nchar(set_short) >
                              # 5 = ...()    4 is a buffer
                              nchar(temp) + 5L + n_digits + 4L,
                            sprintf("%s...(%s)", temp, set_id),
                            set_short)) %>%
  select(-temp) %>%
  relocate(set_id, .before = set) %>%
  `rownames<-`(NULL) %>%
  mutate(collection = case_when(
    database %in% c("MITOCARTA", "REFMET", "PSP", "CELLMARKER") ~ database,
    grepl("^PTMSIGDB", database) ~ "PTMSIGDB",
    grepl("^GO", database) ~ "C5",
    TRUE ~ "C2"
  ),
  collection = factor(collection, levels = unique(collection))) %>%
  relocate(collection, .before = everything()) %>%
  mutate(across(.cols = everything(),
                .fns = ~ structure(.x, names = NULL)))

# str(SET_TO_ID)
#
# count(SET_TO_ID, collection, database) %>%
#   View()

# Save
usethis::use_data(SET_TO_ID, overwrite = TRUE, version = 3, compress = TRUE)
