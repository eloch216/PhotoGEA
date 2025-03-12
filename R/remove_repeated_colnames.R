# A helping function that removes any columns of `exdf_a` that are also in
# `exdf_b`
remove_repeated_colnames <- function(exdf_a, exdf_b) {
    a_in_b <- colnames(exdf_a) %in% colnames(exdf_b)

    a_in_b_names <- colnames(exdf_a)[a_in_b]

    for (cn in a_in_b_names) {
        exdf_a[, cn] <- NULL
    }

    exdf_a
}
