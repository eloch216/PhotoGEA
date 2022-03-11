# Replace any variations of "WT" (e.g. "WT ", "wt", etc) in the event names
# with a standardized version: "WT"
fix_wt <- function(event_names) {
    is_wt <- tolower(trimws(event_names)) == "wt"
    event_names[is_wt] <- "WT"
    return(event_names)
}

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

    # Standardize the WT event name
    licor_exdf[,event_column_name] <- fix_wt(licor_exdf[,event_column_name])

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
