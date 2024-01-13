\name{error_function_c3_aci}

\alias{error_function_c3_aci}

\title{Generate an error function for C3 A-Ci curve fitting}

\description{
  Creates a function that returns an error value (the sum of squared residuals)
  representing the amount of agreement between modeled and measured \code{An}
  values.
}

\usage{
  error_function_c3_aci(
    replicate_exdf,
    fit_options = list(),
    POc = 210000,
    atp_use = 4.0,
    nadph_use = 8.0,
    curvature_cj = 1.0,
    curvature_cjp = 1.0,
    a_column_name = 'A',
    cc_column_name = 'Cc',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    rd_norm_column_name = 'Rd_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    cj_crossover_min = NA,
    cj_crossover_max = NA
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
    column in \code{exdf_obj} of the same name; and a numeric value means that
    the parameter will be set to that value. For example, \code{fit_options =
    list(alpha = 0, Vcmax_at_25 = 'fit', TPU = 'column')} means that
    \code{alpha} will be set to 0, \code{Vcmax_at_25} will be fit, and
    \code{TPU} will be set to the values in the \code{TPU} column of
    \code{exdf_obj}.
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

  \item{a_column_name}{
    The name of the column in \code{replicate_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{cc_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    chloroplastic CO2 concentration in \code{micromol mol^(-1)}.
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
}

\details{
  When fitting A-Ci curves, it is necessary to define a function that calculates
  the sum of squared residuals for a given set of \code{alpha},
  \code{Gamma_star}, \code{J_at_25}, \code{Rd_at_25}, \code{TPU}, and
  \code{Vcmax} values by comparing a model prediction to a measured curve. This
  function will be passed to an optimization algorithm which will determine the
  values that produce the smallest error.

  The \code{error_function_c3_aci} returns such a function, which is based on a
  particular A-Ci curve and a set of fitting options. It is possible to just fit
  a subset of the available fitting parameters; by default, the fitting
  parameters are \code{J_at_25}, \code{Rd_at_25}, \code{TPU}, and \code{Vcmax}.
  This behavior can be changed via the \code{fit_options} argument.

  Sometimes an optimizer will choose biologically unreasonable parameter values
  that nevertheless produce good fits to the supplied assimilation values. A
  common problem is that the fit result may not indicate Ac-limited assimilation
  at low CO2 values, which should be the case for any A-Ci curves measured at
  saturating light. In this case, the optional \code{cj_crossover_min} and
  \code{cj_crossover_max} can be used to constrain the range of \code{Cc} values
  (in ppm) where \code{Aj} is allowed to be the overall rate limiting factor.
  If the crossover from Rubisco-limited to RuBP-regeneration limited
  assimilation occurs outside these bounds (when they are supplied), a heavy
  penalty will be added to the error function, preventing the optimizer from
  choosing those parameter values. See the _Analyzing C3 A-Ci Curves_ vignette
  for an example of how these arguments can be used to improve the quality of a
  fit.

  A penalty is also added for any parameter combination where \code{An} is not a
  number, or where \code{\link{calculate_c3_assimilation}} produces an error.
}

\value{
  A function with one input argument \code{guess}, which should be a numeric
  vector representing values of the parameters to be fitted (which are specified
  by the \code{fit_options} input argument.) Each element of \code{guess} is the
  value of one parameter (arranged in alphabetical order.) For example, with the
  default settings, \code{guess} should contain values of \code{J_at_25},
  \code{Rd_at_25}, \code{TPU}, and \code{Vcmax} (in that order).
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

# Define an error function for one curve from the set
error_fcn <- error_function_c3_aci(
  licor_file[licor_file[, 'species_plot'] == 'tobacco - 1', , TRUE]
)

# Evaluate the error for J_at_25 = 200, Rd_at_25 = 1, TPU = 15, Vcmax_at_25 = 100
error_fcn(c(200, 1, 15, 100))

# Make a plot of error vs. TPU when the other parameters are fixed to the values
# above. As TPU increases, it eventually stops limiting the assimilation rate
# and its value stops influencing the error.
tpu_error_fcn <- function(TPU) {error_fcn(c(200, 1, TPU, 100))}
tpu_seq <- seq(15, 25)

lattice::xyplot(
  sapply(tpu_seq, tpu_error_fcn) ~ tpu_seq,
  type = 'b',
  xlab = 'TPU (micromol / m^2 / s)',
  ylab = 'Sum of squared residuals'
)
}

\concept{exdf}
