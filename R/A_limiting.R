# The functions in this file are for internal use only, so they are not exported
# to the package namespace.

# Helping function for identifying the points in a curve where a potential
# limiting rate is the actual limiting rate. Returns a logical vector.
A_limiting <- function(data_table, an_name, a_name, tol) {
    rel_diff <-
        abs(data_table[, a_name] - data_table[, an_name]) /
            pmax(1e-10, abs(data_table[, an_name]))

    return(!is.na(rel_diff) & rel_diff <= tol)
}
