# The functions in this file are for internal use only, so they are not exported
# to the package namespace.

# Helping function for identifying the points in a curve where a potential
# limiting rate is the actual limiting rate. Returns a logical vector.
A_limiting <- function(assim_exdf, an_name, a_name, tol) {
    if (!is.exdf(assim_exdf)) {
        stop('assim_exdf must be an exdf object')
    }

    rel_diff <-
        abs(assim_exdf[, a_name] - assim_exdf[, an_name]) /
            pmax(1e-10, abs(assim_exdf[, an_name]))

    return(!is.na(rel_diff) & rel_diff <= tol)
}
