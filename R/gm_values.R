# This script includes functions for specifying mesophyll conductance directly
# (as a single value) or by reading mesophyll conductance (gm) table files and
# adding information from these files to Licor data tables.

# Specify a value of mesophyll conductance to include as a column in the data
add_gm_to_licor_data_from_value <- function(
    licor_data,
    gm_value,
    gm_units,
    GM_COLUMN_NAME
)
{
    # Add a column for gm
    licor_data <- specify_variables(
        licor_data,
        c("gm input", GM_COLUMN_NAME, gm_units)
    )

    licor_data[['main_data']][[GM_COLUMN_NAME]] <- gm_value

    return(licor_data)
}

# Read data from a table that contains mesophyll conductance values for each
# genotype. The table should be stored as a csv file with two columns:
# GENOTYPE_COLUMN_NAME (specifying the genotype) and GM_COLUMN_NAME (specifying
# the mesophyll conductance to CO2). The first row should be the column names,
# the second row should be the units for each column, and the remaining rows
# should contain the data. The return value is an exdf object.
read_gm_table <- function(
    filename,
    GENOTYPE_COLUMN_NAME,
    GM_COLUMN_NAME
)
{
    raw_data <- utils::read.csv(filename, header = FALSE)

    table_columns <- raw_data[1,]

    if (!(GENOTYPE_COLUMN_NAME %in% table_columns) || !(GM_COLUMN_NAME %in% table_columns)) {
        stop(
            paste0(
                "The mesophyll conductance data table columns must agree with ",
                "the specified genotype column name (", GENOTYPE_COLUMN_NAME,
                ") and gm column name (", GM_COLUMN_NAME, ")"
            )
        )
    }

    table_units <- raw_data[2,]
    table_data <- raw_data[seq(3, nrow(raw_data)),]
    table_categories <- data.frame(
        matrix(
            data = "gm input",
            nrow = 1,
            ncol = ncol(raw_data)
        ),
        stringsAsFactors = FALSE
    )

    colnames(table_categories) <- table_columns
    colnames(table_units) <- table_columns
    colnames(table_data) <- table_columns

    rownames(table_units) <- NULL
    rownames(table_data) <- NULL

    # Convert the data to numeric values whenever possible
    table_data <- as.data.frame(
        lapply(table_data, try_as_numeric),
        stringsAsFactors = FALSE
    )

    return(
        exdf(
            table_data,
            table_units,
            table_categories
        )
    )
}

# Adds gm values to a list representing Licor data (as created by
# 'read_licor_file' or a similar function). The gm values should be specified
# in a list as produced by the 'read_gm_table' function.
add_gm_to_licor_data_from_table <- function(
    licor_data,
    gm_table,
    GENOTYPE_COLUMN_NAME,
    GM_COLUMN_NAME
)
{
    # Add a new column to the licor data representing gm
    licor_data <- specify_variables(
        licor_data,
        c(
            gm_table[['categories']][[GM_COLUMN_NAME]],
            GM_COLUMN_NAME,
            gm_table[['units']][[GM_COLUMN_NAME]]
        )
    )

    # Fill in the values
    for (i in seq_len(nrow(gm_table[['main_data']]))) {
        genotype <- gm_table[i,GENOTYPE_COLUMN_NAME]
        gm_val <- gm_table[i,GM_COLUMN_NAME]

        licor_data[,GM_COLUMN_NAME][licor_data[,GENOTYPE_COLUMN_NAME] == genotype] <- gm_val

    }

    return(licor_data)
}

# choose_input_gm_table_file: a function for interactively selecting a gm table
# file.
#
# Important note: this function is only available on MS Windows.
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# None
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a vector of file name strings representing absolute paths to Licor Excel files
#
choose_input_gm_table_file <- function()
{
    if (!interactive() | .Platform$OS.type != "windows") {
        stop(
            paste(
                "The `choose_input_gm_table_file` function is only available in",
                "interactive R sessions running in MS Windows"
            )
        )
    }

    utils::choose.files(
        default = "",
        caption = "Select a gm table input file",
        multi = FALSE,
        filters = matrix(
            c(
                "CSV files (*.csv)", "All files (*.*)",
                "*.csv", "*.*"
            ),
            ncol = 2
        ),
        index = 1
    )
}
