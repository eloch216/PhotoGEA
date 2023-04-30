# Define a small number for checking if another number is zero
quadratic_eps_zero <- 1e-15

# Finds the values of `x` where `ax^2 + bx + c = 0`; these are the roots of a
# quadratic equation. They can be found using the quadratic formula:
# `x = (-b +/- sqrt(b^2 - 4ac)) / (2a)`. Often we are only concerned with one of
# the two possible roots, so the return value can be specified as follows:
# - root_type 1: The `+` root
# - root_type 2: The `-` root
# - root_type 3: The larger root
# - root_type 4: The smaller root
#
# Note that if the root term `b^2 - 4ac` is negative, neither root exists; in
# this case, we return NA for any value of `root_type`.
#
# Another edge case is when `a = 0`; in this case, the quadratic equation
# reduces to a linear one (`bx + c = 0`) that can be solved without resorting to
# the quadratic formula: `x = -c/b`. Here there is just one root, and we return
# it for any value of `root_type`.
#
# In the code below, we refer to `a`, `b`, and `c` as `qa`, `qb`, and `qc` to
# avoid any possible confusion with the base R function `c`.
#
# Finally, note that this function expects arguments of length 1. In other
# words, it is not vectorized.
quadratic_root <- function(qa, qb, qc, root_type) {
    if (any(is.na(c(qa, qb, qc)))) {
        # one or more of the inputs is NA, so just return NA
        NA
    } else if (abs(qa) < quadratic_eps_zero) {
        # qa is zero, so just return -qc / qb
        -qc / qb
    } else {
        # qa is nonzero
        root_term <- qb^2 - 4 * qa * qc
        if (root_term < 0) {
            # there are no solutions, so just return NA
            NA
        } else {
            # find both roots and return the appropriate one
            root_plus <- (-qb + sqrt(root_term)) / (2 * qa)
            root_minus <- (-qb - sqrt(root_term)) / (2 * qa)
            if (root_type == 1) {
                root_plus
            } else if (root_type == 2) {
                root_minus
            } else if (root_type == 3) {
                max(root_plus, root_minus)
            } else {
                min(root_plus, root_minus)
            }
        }
    }
}

# The following functions are wrappers of `quadratic_root` that return
# particular roots.
quadratic_root_plus <- function(qa, qb, qc) {
    quadratic_root(qa, qb, qc, 1)
}

quadratic_root_minus <- function(qa, qb, qc) {
    quadratic_root(qa, qb, qc, 2)
}


quadratic_root_max <- function(qa, qb, qc) {
    quadratic_root(qa, qb, qc, 3)
}


quadratic_root_min <- function(qa, qb, qc) {
    quadratic_root(qa, qb, qc, 4)
}

