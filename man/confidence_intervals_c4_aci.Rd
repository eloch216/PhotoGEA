\name{confidence_intervals_c4_aci}

\alias{confidence_intervals_c4_aci}

\title{Calculate confidence intervals for C4 A-Ci fitting parameters}

\description{
  Calculates confidence intervals for parameters estimated by a C4 A-Ci curve
  fit. It is rare for users to call this function directly, because it can be
  automatically applied to each curve when calling \code{\link{fit_c4_aci}}.
}

\usage{
  confidence_intervals_c4_aci(
    replicate_exdf,
    best_fit_parameters,
    lower = list(),
    upper = list(),
    fit_options = list(),
    sd_A = 1,
    relative_likelihood_threshold = 0.147,
    x_etr = 0.4,
    a_column_name = 'A',
    ao_column_name = 'ao',
    ci_column_name = 'Ci',
    gamma_star_column_name = 'gamma_star',
    gmc_norm_column_name = 'gmc_norm',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    kp_column_name = 'Kp',
    oxygen_column_name = 'oxygen',
    qin_column_name = 'Qin',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    vpmax_norm_column_name = 'Vpmax_norm',
    hard_constraints = 0
  )
}

\arguments{
  \item{replicate_exdf}{
    An \code{exdf} object representing one CO2 response curve.
  }

  \item{best_fit_parameters}{
    An \code{exdf} object representing best-fit parameters for the CO2 response
    curve in \code{replicate_exdf}, as calculated by \code{\link{fit_c4_aci}}.
  }

  \item{lower}{
    The same value that was passed to \code{\link{fit_c4_aci}} when generating
    \code{best_fit_parameters}.
  }

  \item{upper}{
    The same value that was passed to \code{\link{fit_c4_aci}} when generating
    \code{best_fit_parameters}.
  }

  \item{fit_options}{
    The same value that was passed to \code{\link{fit_c4_aci}} when generating
    \code{best_fit_parameters}.
  }

  \item{sd_A}{
    The same value that was passed to \code{\link{fit_c4_aci}} when generating
    \code{best_fit_parameters}.
  }

  \item{relative_likelihood_threshold}{
    The threshold value of relative likelihood used to define the boundaries of
    the confidence intervals; see details below.
  }

  \item{x_etr}{
    The fraction of whole-chain electron transport occurring in the mesophyll
    (dimensionless). See Equation 29 from S. von Caemmerer (2021).
  }

  \item{a_column_name}{
    The name of the column in \code{replicate_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ao_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    dimensionless ratio of solubility and diffusivity of O2 to CO2.
  }

  \item{ci_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    intercellular CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{gamma_star_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    dimensionless \code{gamma_star} values.
  }

  \item{gmc_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    mesophyll conductance values (with units of
    \code{normalized to gmc at 25 degrees C}).
  }

  \item{j_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Jmax} values (with units of
    \code{normalized to Jmax at 25 degrees C}).
  }

  \item{kc_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    Michaelis-Menten constant for rubisco carboxylation in \code{microbar}.
  }

  \item{ko_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    Michaelis-Menten constant for rubisco oxygenation in \code{mbar}.
  }

  \item{kp_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    Michaelis-Menten constant for PEP carboxylase carboxylation in
    \code{microbar}.
  }

  \item{oxygen_column_name}{
    The name of the column in \code{exdf_obj} that contains the concentration
    of O2 in the ambient air, expressed as a percentage (commonly 21\% or 2\%);
    the units must be \code{percent}.
  }

  \item{qin_column_name}{
    The name of the column in \code{exdf_obj} that contains values of the
    incident photosynthetically active flux density in
    \code{micromol m^(-2) s^(-1)}.
  }

  \item{rl_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{RL} values (with units of \code{normalized to RL at 25 degrees C}).
  }

  \item{total_pressure_column_name}{
    The name of the column in \code{exdf_obj} that contains the total pressure
    in \code{bar}.
  }

  \item{vcmax_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{Vcmax} values (with units of
    \code{normalized to Vcmax at 25 degrees C}).
  }

  \item{vpmax_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{Vpmax} values (with units of
    \code{normalized to Vpmax at 25 degrees C}).
  }

  \item{hard_constraints}{
    To be passed to \code{\link{calculate_c4_assimilation}}; see that function
    for more details.
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
  where the relative likelihood is above the threshold value (set by the
  \code{relative_likelihood_threshold} input argument). If the threshold
  is set to 0.147, then these intervals are equivalent to 95\% confidence
  intervals in most situations. See the Wikipedia page about
  \href{https://en.wikipedia.org/wiki/Relative_likelihood}{relative likelihood}
  for more information.

  Internally, this function uses \code{\link{error_function_c4_aci}} to
  calculate the negative logarithm of the likelihood (\code{-ln(L)}). It varies
  each fitting parameter independendently to find values where
  \code{ln(L) - ln(p_0) - ln(L_max) = 0}.

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

# Calculate temperature-dependent values of C4 photosynthetic parameters
licor_file <- calculate_temperature_response(licor_file, c4_temperature_param_vc)

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Fit just one curve from the data set
one_result <- fit_c4_aci(
  licor_file[licor_file[, 'species_plot'] == 'maize - 5', , TRUE],
  Ca_atmospheric = 420
)

# Calculate confidence limits for the fit parameters
parameters_with_limits <- confidence_intervals_c4_aci(
    licor_file[licor_file[, 'species_plot'] == 'maize - 5', , TRUE],
    one_result$parameters
)

# View confidence limits and best estimate for Vcmax_at_25
parameters_with_limits[, c('Vcmax_at_25_lower', 'Vcmax_at_25', 'Vcmax_at_25_upper')]
}

\concept{exdf}
