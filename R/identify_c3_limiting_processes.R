identify_c3_limiting_processes <- function(
    exdf_obj,
    a_column_name = 'A_fit',
    ac_column_name = 'Ac',
    aj_column_name = 'Aj',
    ap_column_name = 'Ap',
    tol = 1e-3
)
{
    if (!is.exdf(exdf_obj)) {
        stop('identify_c3_limits requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]]  <- unit_dictionary[['A']]
    required_variables[[ac_column_name]] <- unit_dictionary[['A']]
    required_variables[[aj_column_name]] <- unit_dictionary[['A']]
    required_variables[[ap_column_name]] <- unit_dictionary[['A']]

    check_required_variables(exdf_obj, required_variables, check_NA = FALSE)

    # Determine which process is limiting, allowing for the possibility of
    # co-limitation
    exdf_obj[, 'Ac_limiting'] <- A_limiting(exdf_obj, a_column_name, ac_column_name, tol)
    exdf_obj[, 'Aj_limiting'] <- A_limiting(exdf_obj, a_column_name, aj_column_name, tol)
    exdf_obj[, 'Ap_limiting'] <- A_limiting(exdf_obj, a_column_name, ap_column_name, tol)

    exdf_obj[, 'limiting_process'] <- sapply(seq_len(nrow(exdf_obj)), function(i) {
        if (exdf_obj[i, 'Ac_limiting'] && !exdf_obj[i, 'Aj_limiting'] && !exdf_obj[i, 'Ap_limiting']) {
            'Ac'
        } else if (!exdf_obj[i, 'Ac_limiting'] && exdf_obj[i, 'Aj_limiting'] && !exdf_obj[i, 'Ap_limiting']) {
            'Aj'
        } else if (!exdf_obj[i, 'Ac_limiting'] && !exdf_obj[i, 'Aj_limiting'] && exdf_obj[i, 'Ap_limiting']) {
            'Ap'
        } else {
            'co-limited'
        }
    })

    # Document and return
    document_variables(
        exdf_obj,
        c('identify_c3_limits', 'Ac_limiting',      ''),
        c('identify_c3_limits', 'Aj_limiting',      ''),
        c('identify_c3_limits', 'Ap_limiting',      ''),
        c('identify_c3_limits', 'limiting_process', '')
    )
}
