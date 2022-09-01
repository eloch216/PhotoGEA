calculate_cc <- function(
    licor_exdf,
    a_column_name,
    ca_column_name,
    ci_column_name,
    gm_column_name,
    pa_column_name,
    deltapcham_column_name
)
{
    if (!is.exdf(licor_exdf)) {
        stop("calculate_cc requires an exdf object")
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]] <- "micromol m^(-2) s^(-1)"
    required_variables[[ca_column_name]] <- "micromol mol^(-1)"
    required_variables[[ci_column_name]] <- "micromol mol^(-1)"
    required_variables[[gm_column_name]] <- "mol m^(-2) s^(-1) bar^(-1)"
    required_variables[[pa_column_name]] <- "kPa"
    required_variables[[deltapcham_column_name]] <- "kPa"

    check_required_variables(licor_exdf, required_variables)

    # Define new column names
    cc_column_name <- "Cc"
    drawdown_m_column_name <- "drawdown_m"
    drawdown_s_column_name <- "drawdown_s"

    # Make calculations
    licor_exdf[,cc_column_name] <-
        licor_exdf[,ci_column_name] - licor_exdf[,a_column_name] /
            (licor_exdf[,gm_column_name] * (licor_exdf[,pa_column_name] + licor_exdf[,deltapcham_column_name]) / 100)

    licor_exdf[,drawdown_m_column_name] <-
        licor_exdf[,ci_column_name] - licor_exdf[,cc_column_name]

    licor_exdf[,drawdown_s_column_name] <-
        licor_exdf[,ca_column_name] - licor_exdf[,ci_column_name]

    # Document the columns that were added
    licor_exdf <- document_variables(
        licor_exdf,
        c("calculate_cc", cc_column_name,         "micromol mol^(-1)"),
        c("calculate_cc", drawdown_m_column_name, "micromol mol^(-1)"),
        c("calculate_cc", drawdown_s_column_name, "micromol mol^(-1)")
    )

    return(licor_exdf)
}
