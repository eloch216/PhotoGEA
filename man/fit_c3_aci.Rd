\name{fit_c3_aci}

\alias{fit_c3_aci}

\title{Fits a C3 assimilation model to an experimental curve}

\description{
  Fits a model to an experimentally measured C3 CO2 response curve using the
  data in the \code{exdf} object along with a few other user-supplied
  parameters. This function can accomodate alternative column names for the
  variables taken from the Licor file in case they change at some point in the
  future. This function also checks the units of each required column and will
  produce an error if any units are incorrect.
}

\usage{
  fit_c3_aci(
    replicate_exdf,
    Ca_atmospheric,
    a_column_name = 'A',
    ca_column_name = 'Ca',
    cc_column_name = 'Cc',
    ci_column_name = 'Ci',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    rd_norm_column_name = 'Rd_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    sd_A = 'RMSE',
    POc = 210000,
    atp_use = 4.0,
    nadph_use = 8.0,
    curvature_cj = 1.0,
    curvature_cjp = 1.0,
    OPTIM_FUN = optimizer_nmkb(),
    lower = list(),
    upper = list(),
    fit_options = list(),
    cj_crossover_min = NA,
    cj_crossover_max = NA,
    error_threshold_factor = 1.5,
    calculate_confidence_intervals = FALSE,
    remove_unreliable_param = FALSE
  )
}

\arguments{
  \item{replicate_exdf}{
    An \code{exdf} object representing one CO2 response curve.
  }

  \item{Ca_atmospheric}{
    The atmospheric CO2 concentration (with units of \code{micromol mol^(-1)});
    this will be used by \code{\link{estimate_operating_point}} to estimate the
    operating point.
  }

  \item{a_column_name}{
    The name of the column in \code{replicate_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ca_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    ambient CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{cc_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    chloroplastic CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{ci_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    intercellular CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{j_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{J} values (with units of \code{normalized to J at 25 degrees C}).
  }

  \item{kc_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    Michaelis-Menten constant for rubisco carboxylation in
    \code{micromol mol^(-1)}.
  }

  \item{ko_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    Michaelis-Menten constant for rubisco oxygenation in
    \code{mmol mol^(-1)}.
  }

  \item{rd_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{Rd} values (with units of \code{normalized to Rd at 25 degrees C}).
  }

  \item{total_pressure_column_name}{
    The name of the column in \code{replicate_exdf} that contains the total
    pressure in \code{bar}.
  }

  \item{vcmax_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{Vcmax} values (with units of
    \code{normalized to Vcmax at 25 degrees C}).
  }

  \item{sd_A}{
    A value of the standard deviation of measured \code{A} values, or the name
    of a method for determining the deviation; currently, the only supported
    option is \code{'RMSE'}.
  }

  \item{POc}{
    The partial pressure of O2 in the chloroplast, expressed in \code{microbar}.
    This is often assumed to be the ambient value. For air measurements, this
    would be 21 percent O2, which is about 210000 microbar at standard
    atmospheric pressure. For low oxygen measurements, this would be 2 percent
    O2, which is about 20000 microbar.
  }

  \item{atp_use}{
    The number of ATP molecules used per C3 cycle.
  }

  \item{nadph_use}{
    The number of NADPH molecules used per C3 cycle.
  }

  \item{curvature_cj}{
    A dimensionless quadratic curvature parameter greater than or equal to 0 and
    less than or equal to 1 that sets the degree of co-limitation between
    \code{Wc} and \code{Wj}. A value of 1 indicates no co-limitation.
  }

  \item{curvature_cjp}{
    A dimensionless quadratic curvature parameter greater than or equal to 0 and
    less than or equal to 1 that sets the degree of co-limitation between
    \code{Wcj} and \code{Wp}. A value of 1 indicates no co-limitation.
  }

  \item{OPTIM_FUN}{
    An optimization function that accepts the following input arguments: an
    initial guess, an error function, lower bounds, and upper bounds. It should
    return a list with the following elements: \code{par}, \code{convergence},
    \code{value}, and (optionally) \code{message}. See \code{\link{optimizers}}
    for a list of available options.
  }

  \item{lower}{
    A list of named numeric elements representing lower bounds to use when
    fitting. Values supplied here override the default values (see details
    below). For example, \code{lower = list(Vcmax_at_25 = 10)} sets the lower
    limit for \code{Vcmax_at_25} to 10 micromol / m^2 / s.
  }

  \item{upper}{
    A list of named numeric elements representing upper bounds to use when
    fitting. Values supplied here override the default values (see details
    below). For example, \code{upper = list(Vcmax_at_25 = 200)} sets the upper
    limit for \code{Vcmax_at_25} to 200 micromol / m^2 / s.
  }

  \item{fit_options}{
    A list of named elements representing fit options to use for each parameter.
    Values supplied here override the default values (see details below). Each
    element must be \code{'fit'}, \code{'column'}, or a numeric value. A value
    of \code{'fit'} means that the parameter will be fit; a value of
    \code{'column'} means that the value of the parameter will be taken from a
    column in \code{exdf_obj} of the same name; and a numeric value means that
    the parameter will be set to that value. For example, \code{fit_options =
    list(alpha_g = 0, Vcmax_at_25 = 'fit', Tp = 'column')} means that
    \code{alpha_g} will be set to 0, \code{Vcmax_at_25} will be fit, and
    \code{Tp} will be set to the values in the \code{Tp} column of
    \code{exdf_obj}.
  }

  \item{cj_crossover_min}{
    The minimum value of \code{Cc} (in ppm) where \code{Aj} is allowed to become
    the overall rate-limiting factor. If \code{cj_crossover_min} is set to
    \code{NA}, this restriction will not be applied.
  }

  \item{cj_crossover_max}{
    The maximim value of \code{Cc} (in ppm) where \code{Wj} is allowed to be
    smaller than \code{Wc}. If \code{cj_crossover_max} is set to \code{NA}, this
    restriction will not be applied.
  }

  \item{error_threshold_factor}{
    To be passed to \code{\link{confidence_intervals_c3_aci}} when
    \code{calculate_confidence_intervals} is \code{TRUE}.
  }

  \item{calculate_confidence_intervals}{
    A logical value indicating whether or not to estimate confidence intervals
    for the fitting parameters using \code{\link{confidence_intervals_c3_aci}}.
  }

  \item{remove_unreliable_param}{
    A logical value indicating whether or not to remove unreliable parameter
    estimates; see below for details.
  }
}

\details{
  This function calls \code{\link{calculate_c3_assimilation}} to calculate
  values of net assimilation. The user-supplied optimization function is used to
  vary the values of \code{alpha_g}, \code{J_at_25}, \code{Rd_at_25}, \code{Tp},
  and \code{Vcmax_at_25} to find ones that best reproduce the experimentally
  measured values of net assimilation. By default, the following options are
  used for the fits:
  \itemize{
    \item \code{alpha_g}: lower = 0, upper = 1, fit_option = 0
    \item \code{Gamma_star}: lower = 0, upper = 200, fit_option = \code{'column'}
    \item \code{J_at_25}: lower = 0, upper = 1000, fit_option = \code{'fit'}
    \item \code{Rd_at_25}: lower = 0, upper = 100, fit_option = \code{'fit'}
    \item \code{Tp}: lower = 0, upper = 40, fit_option = \code{'fit'}
    \item \code{Vcmax_at_25}: lower = 0, upper = 1000, fit_option = \code{'fit'}
  }
  With these settings, \code{alpha_g} is set to 0, values of \code{Gamma_star} are
  taken from the \code{Gamma_star} column of \code{replicate_exdf}, and the
  other parameters are fit during the process (see \code{fit_options} above).

  An initial guess for the parameters is generated by calling
  \code{\link{initial_guess_c3_aci}} as follows:
  \itemize{
    \item \code{cc_threshold_rd} is set to 100 micromol / mol.
    \item If \code{alpha_g} is being fit, the \code{alpha_g} argument of
          \code{initial_guess_c3_aci} is set to 0.5; otherwise, the argument is
          set to the value specified by the fit options.
    \item If \code{Gamma_star} is being fit, the \code{Gamma_star} argument of
          \code{initial_guess_c3_aci} is set to 40; otherwise, the argument is
          set to the value specified by the fit options.
  }
  Note that any fixed values specified in the fit options will override the
  values returned by the guessing function.

  The fit is made by creating an error function using
  \code{\link{error_function_c3_aci}} and minimizing its value using
  \code{OPTIM_FUN}, starting from the initial guess described above.

  The photosynthesis model represented by \code{calculate_c3_assimilation} is
  not smooth in the sense that small changes in the input parameters do not
  necessarily cause changes in its outputs. This is related to the final step in
  the calculations, where the overall assimilation rate is taken to be the
  minimum of three enzyme-limited rates. For example, if the assimilation rate
  is never TPU-limited, modifying \code{Tp} will not change the model's
  outputs. For this reason, derivative-based optimizers tend to struggle when
  fitting C3 A-Ci curves. Best results are obtained using derivative-free
  methods. It has been found that \code{dfoptim::nmkb} is often able to find a
  good fit.

  Sometimes the optimizer may choose a set of parameter values where one or more
  of the potential limiting carboxylation rates (\code{Wc}, \code{Wj}, or
  \code{Wp}) is never the smallest rate. In this case, the corresponding
  parameter estimates (\code{Vcmax}, \code{J}, or \code{alpha_g} & \code{Tp})
  will be unreliable. It is also possible that the upper limit of the confidence
  interval for a parameter is infinity; this also indicates an unreliable
  parameter estimate. These criteria are used to determine the reliability of
  each parameter estimate, which is indicated in the \code{Vcmax_trust},
  \code{J_trust}, \code{alpha_g_trust}, and \code{Tp_trust} columns of the
  output from \code{fit_c3_aci}, where a value of \code{0} indicates an
  unreliable estimate and \code{1} indicates a reliable estimate. If
  \code{remove_unreliable_param} is \code{TRUE}, the values of any parameters
  deemed unreliable will be set to \code{NA} in the outputs from
  \code{fit_c3_aci}.

  Once the best-fit parameters have been determined, this function also
  estimates the operating value of `\code{Cc} from the atmospheric CO2
  concentration \code{atmospheric_ca} using
  \code{\link{estimate_operating_point}}, and then uses that value to estimate
  the modeled \code{An} at the operating point via
  \code{\link{calculate_c3_assimilation}}.

  This function assumes that \code{replicate_exdf} represents a single
  C3 A-Ci curve. To fit multiple curves at once, this function is often used
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
            \item Values of fitting parameters at 25 degrees C will be stored in
                  the \code{J_at_25}, \code{Rd_at_25}, and \code{Vcmax_at_25}
                  columns.
            \item The other outputs from \code{\link{calculate_c3_assimilation}}
                  will be stored in columns with the usual names: \code{alpha_g},
                  \code{Tp}, \code{Vcmax_tl}, \code{Rd_tl}, \code{J_tl},
                  \code{Wc}, \code{Wj}, \code{Wp}, \code{Vc}, \code{Ac},
                  \code{Aj}, and \code{Ap}.
          }
    \item \code{parameters}: An \code{exdf} object including the identifiers,
          fitting parameters, and convergence information for the A-Ci curve:
          \itemize{
            \item The number of points where \code{Wc}, \code{Wj}, and \code{Wp}
                  are each the smallest potential carboxylation are stored in
                  the \code{n_Wc_smallest}, \code{n_Wj_smallest}, and
                  \code{n_Wp_smallest} columns.
            \item The best-fit values are stored in the \code{alpha_g},
                  \code{Tp}, \code{J_at_25}, \code{Rd_at_25}, and
                  \code{Vcmax_at_25} columns.
            \item For parameters that depend on leaf temperature, the average
                  leaf-temperature-dependent values are stored in
                  \code{J_tl_avg}, \code{Rd_tl_avg}, and \code{Vcmax_tl_avg}.
            \item Information about the operating point is stored in
                  \code{operating_Cc}, \code{operating_Ci}, \code{operating_An},
                  and \code{operating_An_model}.
            \item The \code{convergence} column indicates whether the fit was
                  successful (\code{==0}) or if the optimizer encountered a
                  problem (\code{!=0}).
            \item The \code{feval} column indicates how many cost function
                  evaluations were required while finding the optimal parameter
                  values.
            \item The residual stats as returned by \code{\link{residual_stats}}
                  are included as columns with the default names: \code{dof},
                  \code{RSS}, \code{RMSE}, etc.
          }
  }
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  system.file('extdata', 'c3_aci_1.xlsx', package = 'PhotoGEA', mustWork = TRUE)
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

# Specify an infinite mesophyll conductance (so `Cc` = `Ci`)
licor_file <- set_variable(
  licor_file,
  'gmc', 'mol m^(-2) s^(-1) bar^(-1)', value = Inf
)

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Calculate Cc
licor_file <- apply_gm(licor_file)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_file <- calculate_arrhenius(licor_file, c3_arrhenius_bernacchi)

# We can fit just one curve from the data set, although it is rare to do this
one_result <- fit_c3_aci(
  licor_file[licor_file[, 'species_plot'] == 'tobacco - 1', , TRUE],
  Ca_atmospheric = 420
)

# We can fit the same curve, but also allow alpha_g and Gamma_star to vary, and
# calculate confidence intervals
one_result_v2 <- fit_c3_aci(
  licor_file[licor_file[, 'species_plot'] == 'tobacco - 1', , TRUE],
  Ca_atmospheric = 420,
  fit_options = list(Gamma_star = 'fit', alpha_g = 'fit'),
  calculate_confidence_intervals = TRUE
)

# Fit all curves in the data set (it is more common to do this)
aci_results <- consolidate(by(
  licor_file,
  licor_file[, 'species_plot'],
  fit_c3_aci,
  Ca_atmospheric = 420
))

# View the fitting parameters for each species / plot
col_to_keep <- c(
  'species', 'plot',                                       # identifiers
  'n_Wc_smallest', 'n_Wj_smallest', 'n_Wp_smallest',       # number of points where
                                                           #   each process is limiting
  'Tp', 'J_at_25', 'Rd_at_25', 'Vcmax_at_25',              # parameters scaled to 25 degrees C
  'J_tl_avg', 'Rd_tl_avg', 'Vcmax_tl_avg',                 # average temperature-dependent values
  'operating_Ci', 'operating_An', 'operating_An_model',    # operating point info
  'dof', 'RSS', 'MSE', 'RMSE', 'RSE',                      # residual stats
  'convergence', 'convergence_msg', 'feval', 'optimum_val' # convergence info
)

aci_results$parameters[ , col_to_keep, TRUE]

# View the fits for each species / plot
lattice::xyplot(
  Ac + Aj + Ap + A_fit + A ~ Ci | species_plot,
  data = aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlab = paste0('Intercellular CO2 concentration (', aci_results$fits$units$Ci, ')'),
  ylab = paste0('Assimilation rate (', aci_results$fits$units$A, ')')
)

# View the residuals for each species / plot
lattice::xyplot(
  A_residuals ~ Ci | species_plot,
  data = aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlab = paste0('Intercellular CO2 concentration (', aci_results$fits$units$Ci, ')'),
  ylab = paste0('Assimilation rate residuals (', aci_results$fits$units$A_residuals, ')')
)

# In the curves above, there are no points where carboxylation is TPU limited.
# We can automatically remove the unreliable Tp estimates.
aci_results <- consolidate(by(
  licor_file,
  licor_file[, 'species_plot'],
  fit_c3_aci,
  Ca_atmospheric = 420,
  remove_unreliable_param = TRUE
))

aci_results$parameters[ , col_to_keep, TRUE]
}

\concept{exdf}
