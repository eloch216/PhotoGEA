extract_tdl_valve <- function(
    tdl_exdf,
    valve_column_name,
    valve_number
)
{
    tdl_data <- tdl_exdf[['main_data']]
    tdl_data[tdl_data[[valve_column_name]] == valve_number,]
}
