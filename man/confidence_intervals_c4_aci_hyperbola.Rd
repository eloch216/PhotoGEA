\name{confidence_intervals_c4_aci_hyperbola}

\alias{confidence_intervals_c4_aci_hyperbola}

\title{Calculate confidence intervals for C4 A-Ci hyperbola fitting parameters}

\description{
  Calculates confidence intervals for parameters estimated by a C4 A-Ci curve
  fit. It is rare for users to call this function directly, because it can be
  automatically applied to each curve when calling
  \code{\link{fit_c4_aci_hyperbola}}.
}

\usage{
  confidence_intervals_c4_aci_hyperbola(
    replicate_exdf,
    best_fit_parameters,
    lower = list(),
    upper = list(),
    fit_options = list(),
    sd_A = 1,
    error_threshold_factor = 0.147,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    hard_constraints = 2
  )
}

\arguments{
  \item{replicate_exdf}{
    An \code{exdf} object representing one CO2 response curve.
  }

  \item{best_fit_parameters}{
    An \code{exdf} object representing best-fit parameters for the CO2 response
    curve in \code{replicate_exdf}, as calculated by
    \code{\link{fit_c4_aci_hyperbola}}.
  }

  \item{lower}{
    The same value that was passed to \code{\link{fit_c4_aci_hyperbola}} when
    generating \code{best_fit_parameters}.
  }

  \item{upper}{
    The same value that was passed to \code{\link{fit_c4_aci_hyperbola}} when
    generating \code{best_fit_parameters}.
  }

  \item{fit_options}{
    The same value that was passed to \code{\link{fit_c4_aci_hyperbola}} when
    generating \code{best_fit_parameters}.
  }

  \item{sd_A}{
    The same value that was passed to \code{\link{fit_c4_aci_hyperbola}} when
    generating \code{best_fit_parameters}.
  }

  \item{error_threshold_factor}{
    A multiplicative factor that sets the threshold value of the error function
    used to define the edges of the confidence intervals; see details below.
  }

  \item{a_column_name}{
    The name of the column in \code{replicate_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ci_column_name}{
    The name of the column in \code{exdf_obj} that contains the intercellular
    CO2 concentration, expressed in \code{micromol mol^(-1)}.
  }

  \item{hard_constraints}{
    To be passed to \code{\link{calculate_c4_assimilation_hyperbola}}; see that
    function for more details.
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
  and upper bounds for each parameter; for example, if \code{Vmax} was fit,
  \code{best_fit_parameters} will contain new columns called \code{Vmax_lower}
  and \code{Vmax_upper}.
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

# Fit just one curve from the data set
one_result <- fit_c4_aci_hyperbola(
  licor_file[licor_file[, 'species_plot'] == 'maize - 5', , TRUE]
)

# Calculate confidence limits for the fit parameters
parameters_with_limits <- confidence_intervals_c4_aci_hyperbola(
    licor_file[licor_file[, 'species_plot'] == 'maize - 5', , TRUE],
    one_result$parameters
)

# View confidence limits and best estimate for Vmax
parameters_with_limits[, c('Vmax_lower', 'Vmax', 'Vmax_upper')]
}

\concept{exdf}
