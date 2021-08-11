calculate_gm <- function(
    licor_file
)
{
    # Add some new columns to the Licor file in preparation for calculating
    # mesophyll conductance
    variables_to_add <- data.frame(
        rbind(
            c("plant specification",  "genotype",       ""),
            c("plant specification",  "event",          ""),
            c("plant specification",  "replicate",      ""),
            c("plant specification",  "original_file",  "")
        ),
        stringsAsFactors = FALSE
    )
    colnames(variables_to_add) <- c("type", "name", "units")


}
