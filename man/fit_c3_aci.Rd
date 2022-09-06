\name{fit_c3_aci}

\alias{fit_c3_aci}

\title{fit_c3_aci}

\description{
  Calculates a linear fit of stomatal conductance vs. the Ball-Berry index using
  the data in the \code{exdf} object. This function can accomodate alternative
  column names for the variables taken from the Licor file in case they change
  at some point in the future. This function also checks the units of each
  required column and will produce an error if any units are incorrect.

  Fits a model to an experimentally measured CO2 response curve using the data
  in the \code{exdf} object along with a few other user-supplied parameters.
  This function can accomodate alternative column names for the variables taken
  from the Licor file in case they change at some point in the future. This
  function also checks the units of each required column and will produce an
  error if any units are incorrect.
}

\usage{
  fit_c3_aci(
    replicate_exdf,
    a_column_name,
    cc_column_name,
    pa_column_name,
    deltapcham_column_name,
    kc_column_name,
    ko_column_name,
    gamma_star_column_name,
    vcmax_norm_column_name,
    rd_norm_column_name,
    j_norm_column_name,
    POc,
    OPTIM_FUN = function(guess, fun, lower, upper) {
      dfoptim::nmkb(guess, fun, lower, upper, control = list(
        tol = 1e-7,
        maxfeval = 2000,
        restarts.max = 10
      ))
    },
    initial_guess = c(10, 100,  0.5, 90),
    lower =         c(0,  0,    0,   0),
    upper =         c(40, 1000, 100, 1000)
  )
}

\arguments{
  \item{replicate_exdf}{
    An \code{exdf} object representing one Ball-Berry curve.
  }

  \item{a_column_name}{
    The name of the column in \code{replicate_exdf} that contains the net
    assimilation (typically \code{'A'}).
  }

  \item{cc_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    chloroplastic CO2 concentration (typically \code{'Cc'}).
  }

  \item{pa_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    atmospheric pressure (typically \code{'Pa'}).
  }

  \item{deltapcham_column_name}{
    The name of the column in \code{replicate_exdf} that contains the chamber
    overpressure (typically \code{'DeltaPcham'}).
  }

  \item{kc_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    Michaelis-Menten constant for rubisco carboxylation (typically \code{'Kc'}).
  }

  \item{ko_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    Michaelis-Menten constant for rubisco oxygenation (typically \code{'Ko'}).
  }

  \item{gamma_star_column_name}{
    The name of the column in \code{replicate_exdf} that contains the Gamma_star
    values (typically \code{'Gamma_star'}).
  }

  \item{vcmax_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    Vcmax values (typically \code{'Vcmax_norm'}).
  }

  \item{rd_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    Rd values (typically \code{'Rd_norm'}).
  }

  \item{j_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    J values (typically \code{'J_norm'}).
  }

  \item{POc}{
    The partial pressure of O2 in the chloroplast, expressed in \code{microbar}.
    This is often assumed to be the ambient value. For air measurements, this
    would be 21 percent O2, which is about 210000 microbar at standard
    atmospheric pressure. For low oxygen measurements, this would be 2 percent
    O2, which is about 20000 microbar.
  }

  \item{OPTIM_FUN}{
    An optimization function that accepts the following input arguments: an
    initial guess, an error function, lower bounds, and upper bounds. It should
    return a list with the following elements: \code{par}, \code{convergence},
    \code{value}, and (optionally) \code{message}. The bounded optimizers from
    the \code{dfoptim} package meet these requirements. The base function
    \code{optim} with the \code{L-BFGS-B} method can also be used.
  }

  \item{initial_guess}{
    A numeric vector containing initial guesses for the values of \code{TPU},
    \code{J}, \code{Rd}, and \code{Vcmax}, in that order.
  }

  \item{lower}{
    A numeric vector containing lower bounds for the values of \code{TPU},
    \code{J}, \code{Rd}, and \code{Vcmax}, in that order.
  }

  \item{upper}{
    A numeric vector containing upper bounds for the values of \code{TPU},
    \code{J}, \code{Rd}, and \code{Vcmax}, in that order.
  }
}

\details{
  This function calls \code{\link{calculate_c3_assimilation}} to calculate
  values of net assimilation. The user-supplied optimization function is used to
  vary the values of \code{TPU}, \code{J}, \code{Rd}, and \code{Vcmax} to find
  ones that best reproduce the experimentally measured values of net
  assimilation.

  The photosynthesis model represented by \code{calculate_c3_assimilation} is
  not smooth in the sense that small changes in the input parameters do not
  necessarily cause changes in its outputs. This is related to the final step in
  the calculations, where the overall assimilation rate is taken to be the
  minimum of three enzyme-limited rates. For example, if the assimilation rate
  is never phosphate-limited, modifying \code{TPU} will not change the model's
  outputs. For this reason, derivative-based optimizers tend to struggle when
  fitting C3 A-Ci curves. Best results are obtained using derivative-free
  methods. It has been found that \code{dfoptim::nmkb} is often able to find a
  good fit.

  This function assumes that \code{replicate_exdf} represents a single
  C3 A-Ci curve. To fit multiple curves at once, this function is often used
  along with \code{\link{by.exdf}} and \code{\link{consolidate}}.
}

\value{
  A list with two elements:
  \itemize{
    \item \code{fits}: An \code{exdf} object including the measured values and
          the fitted values of net assimilation. The fitted values will be
          stored in a column whose name is determined by appending
          \code{'_fits'} to the end of \code{a_column_name}; typically, this
          will be \code{'A_fits'}.
    \item \code{parameters}: An \code{exdf} object including the fitting
          parameters and convergence information for each A-Ci curve. The
          best-fit values are stored in the \code{TPU}, \code{J_at_25},
          \code{Rd_at_25}, and \code{Vcmax_at_25} columns. The
          \code{convergence} column indicates whether the fit was successful
          (\code{==0}) or if the optimizer encountered a problem (\code{!=0}).
  }
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_licor_file(
  system.file('extdata', 'c3_aci_1.xlsx', package = 'PhotoGEA'),
  c(3, 5, 7, 9, 11, 13), 14, 15, 16, 17, 'time'
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

# Calculate Cc
licor_file <- calculate_cc(
  licor_file,
  'A', 'Ca', 'Ci', 'gmc', 'Pa', 'DeltaPcham'
)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_file <- calculate_arrhenius(licor_file, c3_arrhenius_bernacchi)

# Fit just one curve from the data set (it is rare to do this)
one_result <- fit_c3_aci(
  licor_file[licor_file[, 'species_plot'] == 'tobacco - 1', , TRUE],
  'A',
  'Cc',
  'Pa',
  'DeltaPcham',
  'Kc',
  'Ko',
  'Gamma_star',
  'Vcmax_norm',
  'Rd_norm',
  'J_norm',
  210000
)

# Fit all curves in the data set (it is more common to do this)
aci_results <- consolidate(by(
  licor_file,
  licor_file[, 'species_plot'],
  fit_c3_aci,
  'A',
  'Cc',
  'Pa',
  'DeltaPcham',
  'Kc',
  'Ko',
  'Gamma_star',
  'Vcmax_norm',
  'Rd_norm',
  'J_norm',
  210000
))

# View the fitting parameters for each species / plot
col_to_keep <- c(
  'species', 'plot', 'TPU', 'J_at_25', 'Rd_at_25', 'Vcmax_at_25', 'convergence',
  'convergence_msg', 'feval', 'optimum_val'
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
}

\concept{exdf}
