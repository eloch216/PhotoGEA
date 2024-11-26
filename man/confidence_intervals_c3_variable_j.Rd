\name{confidence_intervals_c3_variable_j}

\alias{confidence_intervals_c3_variable_j}

\title{Calculate confidence intervals for C3 Variable J fitting parameters}

\description{
  Calculates confidence intervals for parameters estimated by a C3 A-Ci curve
  fit. It is rare for users to call this function directly, because it can be
  automatically applied to each curve when calling
  \code{\link{fit_c3_variable_j}}.
}

\usage{
  confidence_intervals_c3_variable_j(
    replicate_exdf,
    best_fit_parameters,
    lower = list(),
    upper = list(),
    fit_options = list(),
    sd_A = 1,
    error_threshold_factor = 0.147,
    atp_use = 4.0,
    nadph_use = 8.0,
    curvature_cj = 1.0,
    curvature_cjp = 1.0,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    oxygen_column_name = 'oxygen',
    phips2_column_name = 'PhiPS2',
    qin_column_name = 'Qin',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    tp_norm_column_name = 'Tp_norm',
    vcmax_norm_column_name = 'Vcmax_norm',
    cj_crossover_min = NA,
    cj_crossover_max = NA,
    hard_constraints = 0,
    require_positive_gmc = 'positive_a',
    gmc_max = Inf,
    check_j = TRUE,
    \dots
  )
}

\arguments{
  \item{replicate_exdf}{
    An \code{exdf} object representing one CO2 response curve.
  }

  \item{best_fit_parameters}{
    An \code{exdf} object representing best-fit parameters for the CO2 response
    curve in \code{replicate_exdf}, as calculated by
    \code{\link{fit_c3_variable_j}}.
  }

  \item{lower}{
    The same value that was passed to \code{\link{fit_c3_variable_j}} when
    generating \code{best_fit_parameters}.
  }

  \item{upper}{
    The same value that was passed to \code{\link{fit_c3_variable_j}} when
    generating \code{best_fit_parameters}.
  }

  \item{fit_options}{
    The same value that was passed to \code{\link{fit_c3_variable_j}} when
    generating \code{best_fit_parameters}.
  }

  \item{sd_A}{
    The same value that was passed to \code{\link{fit_c3_variable_j}} when
    generating \code{best_fit_parameters}.
  }

  \item{error_threshold_factor}{
    A multiplicative factor that sets the threshold value of the error function
    used to define the edges of the confidence intervals; see details below.
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

  \item{a_column_name}{
    The name of the column in \code{replicate_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
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

  \item{oxygen_column_name}{
    The name of the column in \code{exdf_obj} that contains the concentration
    of O2 in the ambient air, expressed as a percentage (commonly 21\% or 2\%);
    the units must be \code{percent}.
  }

  \item{phips2_column_name}{
    The name of the column in \code{replicate_exdf} that contains values of the
    operating efficiency of photosystem II (dimensionless).
  }

  \item{qin_column_name}{
    The name of the column in \code{replicate_exdf} that contains values of the
    incident photosynthetically active flux density in
    \code{micromol m^(-2) s^(-1)}.
  }

  \item{rl_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{RL} values (with units of \code{normalized to RL at 25 degrees C}).
  }

  \item{total_pressure_column_name}{
    The name of the column in \code{replicate_exdf} that contains the total
    pressure in \code{bar}.
  }

  \item{tp_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{Tp} values (with units of \code{normalized to Tp at 25 degrees C}).
  }

  \item{vcmax_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{Vcmax} values (with units of
    \code{normalized to Vcmax at 25 degrees C}).
  }

  \item{cj_crossover_min}{
    To be passed to \code{\link{error_function_c3_variable_j}}.
  }

  \item{cj_crossover_max}{
    To be passed to \code{\link{error_function_c3_variable_j}}.
  }

  \item{hard_constraints}{
    To be passed to \code{\link{calculate_c3_assimilation}} and
    \code{\link{calculate_c3_variable_j}}; see those functions for more details.
  }

  \item{require_positive_gmc}{
    To be passed to \code{\link{error_function_c3_variable_j}}.
  }

  \item{gmc_max}{
    To be passed to \code{\link{error_function_c3_variable_j}}.
  }

  \item{check_j}{
    To be passed to \code{\link{error_function_c3_variable_j}}.
  }

  \item{\dots}{
    Additional arguments to be passed to \code{\link{calculate_c3_assimilation}}.
  }
}

\details{
  In maximum likelihood fitting, each set of parameter values has an associated
  likelihood value. If the maximum likelihood is known, then it is also possible
  to define a relative likelihood \code{p} according to \code{p = L / L_max}.
  The set of all parameter values where \code{p} exceeds a threshold value
  \code{p_0} defines a region in parameter space called like a "relative
  likelihood region." When taking one-dimensional cuts through parameter space,
  the boundaries of the relative likelihood region define a relative likelihood
  interval.

  Here we calculate the upper and lower limits of the relative likelihood
  intervals for each parameter. This is done by fixing the other parameters to
  their best-fit values, and varying a single parameter to find the interval
  where the relative likelihood is above the threshold value. If the threshold
  \code{p_0} is set to 0.147, then these intervals are equivalent to 95\%
  confidence intervals in most situations. See the Wikipedia page about
  \href{https://en.wikipedia.org/wiki/Relative_likelihood}{relative likelihood}
  for more information.

  If the upper limit of a confidence interval is found to exceed ten times the
  upper limit specified when fitting that parameter, then the upper limit of the
  condfidence interval is taken to be infinity.
}

\value{
  An \code{exdf} object based on \code{best_fit_parameters} that contains lower
  and upper bounds for each parameter; for example, if \code{Vcmax_at_25} was
  fit, \code{best_fit_parameters} will contain new columns called
  \code{Vcmax_at_25_lower} and \code{Vcmax_at_25_upper}.
}

\examples{
\donttest{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('c3_aci_1.xlsx')
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

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_file <- calculate_temperature_response(licor_file, c3_temperature_param_bernacchi)

# Fit just one curve from the data set
one_result <- fit_c3_variable_j(
  licor_file[licor_file[, 'species_plot'] == 'tobacco - 1', , TRUE],
  Ca_atmospheric = 420
)

# Calculate confidence limits for the fit parameters
parameters_with_limits <- confidence_intervals_c3_variable_j(
    licor_file[licor_file[, 'species_plot'] == 'tobacco - 1', , TRUE],
    one_result$parameters
)

# View confidence limits and best estimate for Vcmax_at_25
parameters_with_limits[, c('Vcmax_at_25_lower', 'Vcmax_at_25', 'Vcmax_at_25_upper')]
}
}

\concept{exdf}
