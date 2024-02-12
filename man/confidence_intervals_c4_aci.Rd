\name{confidence_intervals_c4_aci}

\alias{confidence_intervals_c4_aci}

\title{Calculate confidence intervals for C3 A-Ci fitting parameters}

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
    error_threshold_factor = 1.5,
    ao_column_name = 'ao',
    a_column_name = 'A',
    gamma_star_column_name = 'gamma_star',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    kp_column_name = 'Kp',
    pcm_column_name = 'PCm',
    rd_norm_column_name = 'Rd_norm',
    vcmax_norm_column_name = 'Vcmax_norm',
    vpmax_norm_column_name = 'Vpmax_norm',
    POm = 210000,
    gbs = 0.003,
    Rm_frac = 0.5,
    alpha_psii = 0
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

  \item{error_threshold_factor}{
    A multiplicative factor that sets the threshold value of the error function
    used to define the edges of the confidence intervals; see details below.
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

  \item{pcm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the partial
    pressure of CO2 in the mesophyll, expressed in \code{microbar}.
  }

  \item{rd_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{Rd} values (with units of \code{normalized to Rd at 25 degrees C}).
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

  \item{POm}{
    The partial pressure of O2 in the mesophyll, expressed in \code{microbar}.
    This is often assumed to be the ambient value. For air measurements, this
    would be 21 percent O2, which is about 210000 microbar at standard
    atmospheric pressure. For low oxygen measurements, this would be 2 percent
    O2, which is about 20000 microbar.
  }

  \item{gbs}{
    The bundle sheath conductance to CO2 in \code{mol m^(-2) s^(-1) bar^(-1)}.
  }

  \item{Rm_frac}{
    The fraction of the total mitochondrial respiration that occurs in the
    mesophyll.
  }

  \item{alpha_psii}{
    The fraction of photosystem II activity in the bundle sheath.
  }
}

\details{
  TODO
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
  system.file('extdata', 'c4_aci_1.xlsx', package = 'PhotoGEA', mustWork = TRUE)
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
licor_file <- calculate_arrhenius(licor_file, c4_arrhenius_von_caemmerer)

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Calculate PCm
licor_file <- apply_gm(licor_file, 'C4')

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
