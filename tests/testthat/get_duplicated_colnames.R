get_duplicated_colnames <- function(exdf_obj) {
    cnames <- colnames(exdf_obj)
    duplicated_cnames <- duplicated(cnames) | grepl('.1', cnames, fixed = TRUE)
    dnames <- cnames[duplicated_cnames]
    names_to_ignore <- c('TIME.1', 'time.1', 'hhmmss.1')
    dnames[!dnames %in% names_to_ignore]
}
