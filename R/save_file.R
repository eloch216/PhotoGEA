save_licor_file <- function(licor_file, output_filename) {
    # Make sure any timestamp columns written in POSIXlt form are simply printed
    # as strings
    main_data <- licor_file[['main_data']]
    for (col in colnames(main_data)) {
        if ("POSIXlt" %in% class(main_data[[col]])) {
            main_data[[col]] <- format(main_data[[col]])
        }
    }

    # Combine all the info from the Licor file into one data frame
    all_info <- licor_file[['categories']]
    all_info <- rbind(all_info, colnames(licor_file[['units']]))
    all_info <- rbind(all_info, licor_file[['units']])
    all_info <- rbind(all_info, main_data)

    # Remove column names
    colnames(all_info) <- NULL

    # Save to a csv file
    write.csv(
        all_info,
        file = output_filename,
        row.names = FALSE
    )
}
