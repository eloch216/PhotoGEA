\name{fit_c4_aci_hyperbola}

\alias{fit_c4_aci_hyperbola}

\title{Fits a hyperbolic C4 assimilation model to an experimental curve}

\description{
  Fits an empirical hyperbola model to an experimentally measured C4 A-Ci curve.

  It is possible to fit the following parameters: \code{c4_curvature},
  \code{c4_slope}, \code{rL}, and \code{Vmax}.

  By default, all of these parameters are fit.

  Best-fit parameters are found using maximum likelihood fitting, where the
  optimizer (\code{optim_fun}) is used to minimize the error function (defined
  by \code{\link{error_function_c4_aci_hyperbola}}).

  Once best-fit parameters are found, confidence intervals are calculated
  using \code{\link{confidence_intervals_c4_aci_hyperbola}}.

  See below for more details.
}

\usage{
  fit_c4_aci_hyperbola(
    replicate_exdf,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    sd_A = 'RMSE',
    optim_fun = optimizer_nmkb(1e-7),
    lower = list(),
    upper = list(),
    fit_options = list(),
    relative_likelihood_threshold = 0.147,
    hard_constraints = 0,
    calculate_confidence_intervals = TRUE
  )
}

\arguments{
  \item{replicate_exdf}{
    An \code{exdf} object representing one CO2 response curve.
  }

  \item{a_column_name}{
    The name of the column in \code{replicate_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ci_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    intercellular CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{sd_A}{
    A value of the standard deviation of measured \code{A} values, or the name
    of a method for determining the deviation; currently, the only supported
    option is \code{'RMSE'}.
  }

  \item{optim_fun}{
    An optimization function that accepts the following input arguments: an
    initial guess, an error function, lower bounds, and upper bounds. It should
    return a list with the following elements: \code{par}, \code{convergence},
    \code{feval}, and \code{convergence_msg}. See \code{\link{optimizers}}
    for a list of available options.
  }

  \item{lower}{
    A list of named numeric elements representing lower bounds to use when
    fitting. Values supplied here override the default values (see details
    below). For example, \code{lower = list(Vmax = 10)} sets the lower limit for
    \code{Vmax} to 10 micromol / m^2 / s.
  }

  \item{upper}{
    A list of named numeric elements representing upper bounds to use when
    fitting. Values supplied here override the default values (see details
    below). For example, \code{upper = list(Vmax = 200)} sets the upper limit
    for \code{Vmax} to 200 micromol / m^2 / s.
  }

  \item{fit_options}{
    A list of named elements representing fit options to use for each parameter.
    Values supplied here override the default values (see details below). Each
    element must be \code{'fit'}, \code{'column'}, or a numeric value. A value
    of \code{'fit'} means that the parameter will be fit; a value of
    \code{'column'} means that the value of the parameter will be taken from a
    column in \code{exdf_obj} of the same name; and a numeric value means that
    the parameter will be set to that value. For example,
    \code{fit_options = list(rL = 0, Vmax = 'fit', c4_curvature = 'column')}
    means that \code{rL} will be set to 0, \code{Vmax} will be fit, and
    \code{c4_curvature} will be set to the values in the \code{c4_curvature}
    column of \code{replicate_exdf}.
  }

  \item{relative_likelihood_threshold}{
    To be passed to \code{\link{confidence_intervals_c4_aci_hyperbola}} when
    \code{calculate_confidence_intervals} is \code{TRUE}.
  }

  \item{hard_constraints}{
    To be passed to \code{\link{calculate_c4_assimilation_hyperbola}}; see that
    function for more details.
  }

  \item{calculate_confidence_intervals}{
    A logical value indicating whether or not to estimate confidence intervals
    for the fitting parameters using
    \code{\link{confidence_intervals_c4_aci_hyperbola}}.
  }
}

\details{
  This function calls \code{\link{calculate_c4_assimilation_hyperbola}} to
  calculate values of net assimilation. The user-supplied optimization function
  is used to vary the values of \code{c4_curvature}, \code{c4_slope}, \code{rL},
  and \code{Vmax} to find ones that best reproduce the experimentally measured
  values of net assimilation. By default, the following options are used for the
  fits:
  \itemize{
    \item \code{c4_curvature}: lower = -10, upper = 10, fit_option = \code{'fit'}
    \item \code{c4_slope}: lower = -50, upper = 1000, fit_option = \code{'fit'}
    \item \code{rL}: lower = -10, upper = 100, fit_option = \code{'fit'}
    \item \code{Vmax}: lower = -50, upper = 1000, fit_option = \code{'fit'}
  }
  With these settings, all of the parameters are fit during the process (see
  \code{fit_options} above). The bounds are chosen liberally to avoid any bias.

  An initial guess for the parameters is generated by calling
  \code{\link{initial_guess_c4_aci_hyperbola}}. Note that any fixed values
  specified in the fit options will override the values returned by the guessing
  function.

  The fit is made by creating an error function using
  \code{\link{error_function_c4_aci_hyperbola}} and minimizing its value using
  \code{optim_fun}, starting from the initial guess described above.  The
  \code{\link{optimizer_nmkb}} optimizer is used by default since it has been
  found to reliably return great fits. However, it is a fast optimizer that can
  get stuck in local minima. If it seems to be returning bad fits, consider
  using the \code{\link{optimizer_deoptim}} optimizer instead, but be aware that
  the fits will take more time to complete.

  Unlike the model represented by \code{\link{calculate_c4_assimilation}}, the
  model in \code{\link{calculate_c4_assimilation_hyperbola}} is smooth in the
  sense that small changes in the input parameters cause small changes in its
  outputs. Because of this, it is a fairly easy model to fit.

  This function assumes that \code{replicate_exdf} represents a single
  C4 A-Ci curve. To fit multiple curves at once, this function is often used
  along with \code{\link{by.exdf}} and \code{\link{consolidate}}.
}

\value{
  A list with two elements:
  \itemize{
    \item \code{fits}: An \code{exdf} object including the original contents of
          \code{replicate_exdf} along with several new columns:
          \itemize{
            \item The fitted values of net assimilation will be stored in a
                  column whose name is determined by appending \code{'_fit'} to
                  the end of \code{a_column_name}; typically, this will be
                  \code{'A_fit'}.
            \item Residuals (measured - fitted) will be stored in a column whose
                  name is determined by appending \code{'_residuals'} to the end
                  of \code{a_column_name}; typically, this will be
                  \code{'A_residuals'}.
            \item Values of fitting parameters will be stored in the
                  \code{c4_curvature}, \code{c4_slope}, \code{rL}, and
                  \code{Vmax} columns.
            \item The other outputs from
                  \code{\link{calculate_c4_assimilation_hyperbola}} will be
                  stored in columns with the usual names: \code{Ag},
                  \code{Ainitial}, \code{Amax}, \code{An}, \code{c4_curvature},
                  \code{c4_slope}, \code{rL}, \code{Vinitial}, \code{Vmax}, and
                  \code{c4_assimilation_hyperbola_msg}.
          }
    \item \code{fits_interpolated}: An \code{exdf} object including the
          calculated assimilation rates at a fine spacing of \code{Ci} values
          (step size of 1 \code{micromol mol^(-1)}).
    \item \code{parameters}: An \code{exdf} object including the identifiers,
          fitting parameters, and convergence information for the A-Ci curve:
          \itemize{
            \item The best-fit values are stored in the \code{c4_curvature},
                  \code{c4_slope}, \code{rL}, and \code{Vmax}. If
                  \code{calculate_confidence_intervals} is \code{TRUE}, upper
                  and lower limits for each of these parameters will also be
                  included.
            \item The \code{convergence} column indicates whether the fit was
                  successful (\code{==0}) or if the optimizer encountered a
                  problem (\code{!=0}).
            \item The \code{feval} column indicates how many cost function
                  evaluations were required while finding the optimal parameter
                  values.
            \item The residual stats as returned by \code{\link{residual_stats}}
                  are included as columns with the default names: \code{dof},
                  \code{RSS}, \code{RMSE}, etc.
            \item The Akaike information criterion is included in the \code{AIC}
                  column.
          }
  }
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('c4_aci_1.xlsx')
)

# Define a new column that uniquely identifies each curve
licor_file[, 'species_plot'] <-
  paste(licor_file[, 'species'], '-', licor_file[, 'plot'] )

# Organize the data
licor_file <- organize_response_curve_data(
    licor_file,
    'species_plot',
    c(9, 10, 16),
    'CO2_r_sp'
)

# Fit just one curve from the data set (it is rare to do this).
one_result <- fit_c4_aci_hyperbola(
  licor_file[licor_file[, 'species_plot'] == 'maize - 5', , TRUE]
)

# Fit all curves in the data set (it is more common to do this)
aci_results <- consolidate(by(
  licor_file,
  licor_file[, 'species_plot'],
  fit_c4_aci_hyperbola
))

# View the fitting parameters for each species / plot
col_to_keep <- c(
  'species', 'plot',                                       # identifiers
  'c4_curvature', 'c4_slope', 'rL', 'Vmax',                # best estimates for parameter values
  'dof', 'RSS', 'MSE', 'RMSE', 'RSE',                      # residual stats
  'convergence', 'convergence_msg', 'feval', 'optimum_val' # convergence info
)

aci_results$parameters[ , col_to_keep, TRUE]

# View the fits for each species / plot
plot_c4_aci_hyperbola_fit(aci_results, 'species_plot', ylim = c(0, 100))

# View the residuals for each species / plot
lattice::xyplot(
  A_residuals ~ Ci | species_plot,
  data = aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlab = paste('Intercellular CO2 concentration [', aci_results$fits$units$Ci, ']'),
  ylab = paste('Assimilation rate residuals [', aci_results$fits$units$A_residuals, ']')
)
}

\concept{exdf}
