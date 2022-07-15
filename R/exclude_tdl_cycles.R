exclude_tdl_cycles <- function(tdl_exdf, cycles_to_exclude) {
    tdl_data <- tdl_exdf[['main_data']]
    tdl_data <- tdl_data[!tdl_data[['cycle_num']] %in% cycles_to_exclude,]
    tdl_exdf[['main_data']] <- tdl_data
    return(tdl_exdf)
}
