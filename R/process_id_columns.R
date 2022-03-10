process_id_columns <- function(
    licor_exdf,
    event_column_name,
    rep_column_name,
    unique_id_column_name
)
{
    # Make sure the required columns are defined
    required_columns <- c(
        event_column_name,
        rep_column_name
    )

    check_required_columns(licor_exdf, required_columns)

    # Make sure both columns are treated as character data and trim any
    # whitespace
    licor_exdf[,event_column_name] <-
        trimws(as.character(licor_exdf[,event_column_name]))

    licor_exdf[,rep_column_name] <-
        trimws(as.character(licor_exdf[,rep_column_name]))

    # Replace any variations of "WT" (e.g. "WT ", "wt", etc) in the event column
    # with a standardized version: "WT"
    is_wt <- tolower(trimws(licor_exdf[,event_column_name])) == "wt"
    licor_exdf[,event_column_name][is_wt] <- "WT"

    # Add a new column that uniquely identifies each measurement by its event
    # and replicate names
    licor_exdf[,unique_id_column_name] <-
        paste(licor_exdf[,event_column_name], licor_exdf[,rep_column_name])

    # Document the column that was added
    licor_exdf <- specify_variables(
        licor_exdf,
        c("process_id_columns", unique_id_column_name, "")
    )

    return(licor_exdf)
}
