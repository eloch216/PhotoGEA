# This is just a helper function that identifies points in an exdf object or
# data frame that should be included in curve fits. If there is a column in the
# table called `include_when_fitting` consisting of logical (TRUE/FALSE) values,
# that column is used to determine which points are included when fitting.
# Otherwise, all points are used.
points_for_fitting <- function(curve_data) {
    if ('include_when_fitting' %in% colnames(curve_data)) {
        as.logical(curve_data[, 'include_when_fitting'])
    } else {
        rep_len(TRUE, nrow(curve_data))
    }
}
