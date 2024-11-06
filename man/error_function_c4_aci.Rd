\name{error_function_c4_aci}

\alias{error_function_c4_aci}

\title{Generate an error function for C4 A-Ci curve fitting}

\description{
  Creates a function that returns an error value (the negative of the natural
  logarithm of the likelihood) representing the amount of agreement between
  modeled and measured \code{An} values.
}

\usage{
  error_function_c4_aci(
    replicate_exdf,
    fit_options = list(),
    sd_A = 1,
    absorptance = 0.85,
    f_spectral = 0.15,
    rho = 0.5,
    theta = 0.7,
    x_etr = 0.4,
    ao_column_name = 'ao',
    a_column_name = 'A',
    gamma_star_column_name = 'gamma_star',
    jmax_norm_column_name = 'Jmax_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    kp_column_name = 'Kp',
    oxygen_column_name = 'oxygen',
    pcm_column_name = 'PCm',
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

  \item{fit_options}{
    A list of named elements representing fit options to use for each parameter.
    Values supplied here override the default values (see details below). Each
    element must be \code{'fit'}, \code{'column'}, or a numeric value. A value
    of \code{'fit'} means that the parameter will be fit; a value of
    \code{'column'} means that the value of the parameter will be taken from a
    column in \code{replicate_exdf} of the same name; and a numeric value means
    that the parameter will be set to that value. For example,
    \code{fit_options = list(RL_at_25 = 0, Vcmax_at_25 = 'fit',
    Vpmax_at_25 = 'column')} means that \code{RL_at_25} will be set to 0,
    \code{Vcmax_at_25} will be fit, and \code{Vpmax_at_25} will be set to the
    values in the \code{Vpmax_at_25} column of \code{replicate_exdf}.
  }

  \item{sd_A}{
    The standard deviation of the measured values of the net CO2 assimilation
    rate, expressed in units of \code{micromol m^(-2) s^(-1)}. If
    \code{sd_A} is not a number, then there must be a column in \code{exdf_obj}
    called \code{sd_A} with appropriate units. A numeric value supplied here
    will overwrite the values in the \code{sd_A} column of \code{exdf_obj} if it
    exists.
  }

  \item{absorptance}{
    The leaf absorptance (dimensionless). See Equation 35 from S. von Caemmerer
    (2021).
  }

  \item{f_spectral}{
    The spectral quality adjustment factor (dimensionless). See Equation 35
    from S. von Caemmerer (2021).
  }

  \item{rho}{
    The fraction of light absorbed by photosystem II rather than photosystem I
    (dimensionless). See Equation 35 from S. von Caemmerer (2021).
  }

  \item{theta}{
    An empirical curvature factor (dimensionless). See Equation 34 from S. von
    Caemmerer (2021).
  }

  \item{x_etr}{
    The fraction of whole-chain electron transport occurring in the mesophyll
    (dimensionless). See Equation 29 from S. von Caemmerer (2021).
  }

  \item{ao_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    dimensionless ratio of solubility and diffusivity of O2 to CO2.
  }

  \item{a_column_name}{
    The name of the column in \code{replicate_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{gamma_star_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    dimensionless \code{gamma_star} values.
  }

  \item{jmax_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Jmax} values (with units of
    \code{normalized to Jmax at its optimal temperature}).
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

  \item{pcm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the partial
    pressure of CO2 in the mesophyll, expressed in \code{microbar}.
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
  When fitting A-Ci curves, it is necessary to define a function that calculates
  the likelihood of a given set of \code{alpha_psii}, \code{gbs},
  \code{Jmax_at_opt}, \code{RL_at_25}, \code{Rm_frac}, \code{Vcmax_at_25},
  \code{Vpmax_at_25}, and \code{Vpr} values by comparing a model prediction to a
  measured curve. This function will be passed to an optimization algorithm
  which will determine the values that produce the smallest error.

  The \code{error_function_c4_aci} returns such a function, which is based on a
  particular A-Ci curve and a set of fitting options. It is possible to just fit
  a subset of the available fitting parameters; by default, the fitting
  parameters are \code{RL_at_25}, \code{Vcmax_at_25}, and \code{Vpmax_at_25}.
  This behavior can be changed via the \code{fit_options} argument.

  For practical reasons, the function actually returns values of \code{-ln(L)},
  where \code{L} is the likelihood. The logarithm of \code{L} is simpler to
  calculate than \code{L} itself, and the minus sign converts the problem from
  a maximization to a minimization, which is important because most optimizers
  are designed to minimize a value.

  A penalty is added to the error value for any parameter combination where
  \code{An} is not a number, or where \code{\link{calculate_c4_assimilation}}
  produces an error.
}

\value{
  A function with one input argument \code{guess}, which should be a numeric
  vector representing values of the parameters to be fitted (which are specified
  by the \code{fit_options} input argument.) Each element of \code{guess} is the
  value of one parameter (arranged in alphabetical order.) For example, with the
  default settings, \code{guess} should contain values of \code{RL_at_25},
  \code{Vcmax_at_25}, and \code{Vpmax_at_25} (in that order).
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
licor_file <- calculate_temperature_response_arrhenius(licor_file, c4_arrhenius_von_caemmerer)

licor_file <- calculate_temperature_response_gaussian(licor_file, c4_peaked_gaussian_von_caemmerer)

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Calculate PCm
licor_file <- apply_gm(licor_file, 'C4')

# Define an error function for one curve from the set
error_fcn <- error_function_c4_aci(
  licor_file[licor_file[, 'species_plot'] == 'maize - 5', , TRUE]
)

# Evaluate the error for RL_at_25 = 0, Vcmax_at_25 = 35, Vpmax_at_25 = 180
error_fcn(c(0, 35, 180))

# Make a plot of error vs. Vcmax_at_25 when the other parameters are fixed to
# the values above.
vcmax_error_fcn <- function(Vcmax_at_25) {error_fcn(c(0, Vcmax_at_25, 180))}
vcmax_seq <- seq(20, 50)

lattice::xyplot(
  sapply(vcmax_seq, vcmax_error_fcn) ~ vcmax_seq,
  type = 'b',
  xlab = 'Vcmax at 25 degrees C (micromol / m^2 / s)',
  ylab = 'Negative log likelihood (dimensionless)'
)
}

\concept{exdf}
