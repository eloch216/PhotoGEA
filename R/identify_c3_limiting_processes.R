identify_c3_limiting_processes <- function(
    data_table,
    a_column_name = 'A_fit',
    ac_column_name = 'Ac',
    aj_column_name = 'Aj',
    ap_column_name = 'Ap',
    tol = 1e-3
)
{
    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]]  <- unit_dictionary('A')
    required_variables[[ac_column_name]] <- unit_dictionary('A')
    required_variables[[aj_column_name]] <- unit_dictionary('A')
    required_variables[[ap_column_name]] <- unit_dictionary('A')

    check_required_variables(data_table, required_variables, check_NA = FALSE)

    # Determine which process is limiting, allowing for the possibility of
    # co-limitation
    data_table[, 'Ac_limiting'] <- A_limiting(data_table, a_column_name, ac_column_name, tol)
    data_table[, 'Aj_limiting'] <- A_limiting(data_table, a_column_name, aj_column_name, tol)
    data_table[, 'Ap_limiting'] <- A_limiting(data_table, a_column_name, ap_column_name, tol)

    data_table[, 'limiting_process'] <- 'unknown' # Should never occur

    data_table[data_table[, 'Ac_limiting'] == TRUE  & data_table[, 'Aj_limiting'] == FALSE & data_table[, 'Ap_limiting'] == FALSE, 'limiting_process'] <- 'Ac'
    data_table[data_table[, 'Ac_limiting'] == FALSE & data_table[, 'Aj_limiting'] == TRUE  & data_table[, 'Ap_limiting'] == FALSE, 'limiting_process'] <- 'Aj'
    data_table[data_table[, 'Ac_limiting'] == FALSE & data_table[, 'Aj_limiting'] == FALSE & data_table[, 'Ap_limiting'] == TRUE,  'limiting_process'] <- 'Ap'
    data_table[data_table[, 'Ac_limiting'] == TRUE  & data_table[, 'Aj_limiting'] == TRUE  & data_table[, 'Ap_limiting'] == FALSE, 'limiting_process'] <- 'co-limited (Ac and Aj)'
    data_table[data_table[, 'Ac_limiting'] == TRUE  & data_table[, 'Aj_limiting'] == FALSE & data_table[, 'Ap_limiting'] == TRUE,  'limiting_process'] <- 'co-limited (Ac and Ap)'
    data_table[data_table[, 'Ac_limiting'] == FALSE & data_table[, 'Aj_limiting'] == TRUE  & data_table[, 'Ap_limiting'] == TRUE,  'limiting_process'] <- 'co-limited (Aj and Ap)'
    data_table[data_table[, 'Ac_limiting'] == TRUE  & data_table[, 'Aj_limiting'] == TRUE  & data_table[, 'Ap_limiting'] == TRUE,  'limiting_process'] <- 'co-limited (Ac, Aj, and Ap)'

    # Document and return
    document_variables(
        data_table,
        c('identify_c3_limits', 'Ac_limiting',      ''),
        c('identify_c3_limits', 'Aj_limiting',      ''),
        c('identify_c3_limits', 'Ap_limiting',      ''),
        c('identify_c3_limits', 'limiting_process', '')
    )
}
