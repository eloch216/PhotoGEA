residual_stats <- function(fit_residuals, units, nparam) {
    # Only use points that are not NA
    fit_residuals <- fit_residuals[!is.na(fit_residuals)]

    npts <- length(fit_residuals)
    dof <- npts - nparam

    RSS <- sum(fit_residuals^2)
    MSE <- RSS / npts
    RMSE <- sqrt(MSE)
    RSE <- if (dof > 0) {sqrt(RSS / dof)} else {NA}

    result <- exdf(data.frame(
        npts = npts,
        nparam = nparam,
        dof = dof,
        RSS = RSS,
        MSE = MSE,
        RMSE = RMSE,
        RSE = RSE,
        stringsAsFactors = FALSE
    ))

    document_variables(
        result,
        c('residual_stats', 'npts',   'NA'),
        c('residual_stats', 'nparam', 'NA'),
        c('residual_stats', 'dof',    'NA'),
        c('residual_stats', 'RSS',    paste0('(', units, ')^2')),
        c('residual_stats', 'MSE',    paste0('(', units, ')^2')),
        c('residual_stats', 'RMSE',   units),
        c('residual_stats', 'RSE',    units)
    )
}
