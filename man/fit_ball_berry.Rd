\name{fit_ball_berry}

\alias{fit_ball_berry}

\title{Fits the Ball-Berry model to an experimental curve}

\description{
  Calculates a linear fit of stomatal conductance vs. the Ball-Berry index using
  the data in the \code{exdf} object. This function can accomodate alternative
  column names for the variables taken from the Licor file in case they change
  at some point in the future. This function also checks the units of each
  required column and will produce an error if any units are incorrect.
}

\usage{
  fit_ball_berry(
    replicate_exdf,
    bb_index_column_name = 'bb_index',
    gsw_column_name = 'gsw'
  )
}

\arguments{
  \item{replicate_exdf}{
    An \code{exdf} object representing one Ball-Berry curve.
  }

  \item{bb_index_column_name}{
    The name of the column in \code{replicate_exdf} that contains the Ball-Berry
    index in \code{mol m^(-2) s^(-1)}.
  }

  \item{gsw_column_name}{
    The name of the column in \code{replicate_exdf} that contains the stomatal
    conductance to water vapor in \code{mol m^(-2) s^(-1)}.
  }
}

\details{
  The Ball-Berry model is a simple way to describe the response of a leaf's
  stomata to its assimilation rate and local environmental conditions.
  Specifically, it predicts stomatal conductance to water vapor using the
  following equation:

  \code{gsw = bb_0 + bb_1 * A * h_s / C_s}

  where \code{gsw} is the stomatal conductance, \code{A} is the net assimilation
  rate, \code{h_s} is the relative humidity at the leaf surface, and \code{C_s}
  is the CO2 concentration at the leaf surface. The term \code{A * h_s / C_s} is
  commonly referred to as the Ball-Berry index, while the intercept
  (\code{bb_0}) and slope (\code{bb_1}) of the linear relationship are the
  Ball-Berry parameters which describe the stomatal response.

  Although this model is certainly an oversimplification, it does encode some
  important stomatal responses. For example, when humidity is low, the stomata
  close, reducing stomatal conductance. Likewise, if the CO2 concentration
  around the leaf is depleted, the stomata open to allow more CO2 to diffuse
  into the leaf's interior, increasing somatal conductance. For more information
  about this model and some possible alternatives, see the following papers:
  \itemize{
    \item Ball, J. T., Woodrow, I. E. and Berry, J. A. "A Model Predicting
          Stomatal Conductance and its Contribution to the Control of
          Photosynthesis under Different Environmental Conditions." in "Progress
          in Photosynthesis Research: Volume 4" (1986)
          [\doi{10.1007/978-94-017-0519-6_48}].
    \item Tardieu, F. and Davies, W. J. "Integration of hydraulic and chemical
          signalling in the control of stomatal conductance and water status of
          droughted plants." Plant, Cell & Environment 16, 341–349 (1993).
          [\doi{10.1111/j.1365-3040.1993.tb00880.x}].
    \item Leuning, R. "A critical appraisal of a combined
          stomatal-photosynthesis model for C3 plants." Plant, Cell &
          Environment 18, 339–355 (1995)
          [\doi{10.1111/j.1365-3040.1995.tb00370.x}].
    \item Dewar, R. C. "The Ball–Berry–Leuning and Tardieu–Davies stomatal
          models: synthesis and extension within a spatially aggregated picture
          of guard cell function." Plant, Cell & Environment 25, 1383–1398 (2002).
          [\doi{10.1046/j.1365-3040.2002.00909.x}].
  }

  Ball-Berry parameters are typically determined by measuring a Ball-Berry
  curve, where one or more of the factors that influence the Ball-Berry index
  is systematically varied across a range of values. At each value, care is
  taken that net assimilation and stomatal conductance have reached their
  steady-state values, and then those values are recorded. Then, a linear fit of
  the experimentally observed stomatal conductances as a function of the
  Ball-Berry index is performed to extract estimates for the Ball-Berry
  intercept and slope.
  
  This function uses \code{\link[stats]{lm}} to perform the fit.

  This function assumes that \code{replicate_exdf} represents a single
  Ball-Berry curve. To fit multiple curves at once, this function is often used
  along with \code{\link{by.exdf}} and \code{\link{consolidate}}.
}

\value{
  A list with two elements:
  \itemize{
    \item \code{fits}: An \code{exdf} object including the measured values and
          the fitted values of stomatal conductance. The fitted values will be
          stored in a column whose name is determined by appending
          \code{'_fits'} to the end of \code{gsw_column_name}; typically, this
          will be \code{'gsw_fits'}. Also includes residuals in the
          \code{gsw_residuals} column and values of the Ball-Berry slope and
          intercept.
    \item \code{parameters}: An \code{exdf} object including the fitting
          parameters and R-squared values. The Ball-Berry intercept is stored in
          the \code{bb_intercept} column and the Ball-Berry slope is stored in
          the \code{bb_slope} column. Their standard errors are stored in the
          \code{bb_intercept_err} and \code{bb_slope_err} columns. The R-squared
          value and p-value for the fit are stored in the \code{r_squared} and
          \code{p_value} columns. Other statistical descriptors of the fit as
          calculated by \code{\link{residual_stats}} are also included.
  }
}

\examples{
# Read an example Licor file included in the PhotoGEA package, calculate
# additional gas properties, calculate the Ball-Berry index, define a new column
# that uniquely identifies each curve, and then perform a fit to extract the
# Ball-Berry parameters from each curve.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- calculate_total_pressure(licor_file)

licor_file <- calculate_gas_properties(licor_file)

licor_file[,'species_plot'] <-
  paste(licor_file[,'species'], '-', licor_file[,'plot'])

licor_file <- calculate_ball_berry_index(licor_file)

# Fit just one curve from the data set (it is rare to do this)
one_result <- fit_ball_berry(
  licor_file[licor_file[, 'species_plot'] == 'soybean - 1a', , TRUE]
)

# Fit all curves in the data set (it is more common to do this)
bb_results <- consolidate(by(
  licor_file,
  licor_file[, 'species_plot'],
  fit_ball_berry
))

# View the fitting parameters for each species / plot
col_to_keep <- c('species', 'plot', 'species_plot', 'bb_intercept', 'bb_slope', 'r_squared')
bb_results$parameters[ , col_to_keep]

# View the fits for each species / plot
plot_ball_berry_fit(bb_results, 'species_plot')
}

\concept{exdf}
