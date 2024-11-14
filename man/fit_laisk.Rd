\name{fit_laisk}

\alias{fit_laisk}

\title{Calculate RL and Ci_star using the Laisk method}

\description{
  Uses the Laisk method to estimate \code{Ci_star} and \code{RL}. This function
  can accomodate alternative colum names for the variables taken from log files
  in case they change at some point in the future. This function also checks the
  units of each required column and will produce an error if any units are
  incorrect.
}

\usage{
  fit_laisk(
    replicate_exdf,
    ci_lower = 40,  # ppm
    ci_upper = 120, # ppm
    a_column_name = 'A',
    ci_column_name = 'Ci',
    ppfd_column_name = 'PPFD'
  )
}

\arguments{
  \item{replicate_exdf}{
    An \code{exdf} object containing multiple A-Ci curves measured at different
    levels of incident photosynthetically active photon flux density (PPFD).
  }

  \item{ci_lower}{
    Lower end of \code{Ci} range used for linear fits of \code{An} vs.
    \code{Ci}.
  }

  \item{ci_upper}{
    Upper end of \code{Ci} range used for linear fits of \code{An} vs.
    \code{Ci}.
  }

  \item{a_column_name}{
    The name of the column in \code{replicate_exdf} that contains the net CO2
    assimilation rate \code{An} in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ci_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    intercellular CO2 concentration \code{Ci} in \code{micromol mol^(-1)}.
  }

  \item{ppfd_column_name}{
    The name of the column in \code{replicate_exdf} that can be used to split it
    into individual response curves.
  }
}

\details{
  The Laisk method is a way to estimate \code{RL} and \code{Ci_star} for a C3
  plant. Definitions of these quantities and a description of the theory
  underpinning this method is given below.

  For a C3 plant, the net CO2 assimilation rate \code{An} is given by

  \code{An = Vc - Rp - RL},

  where \code{Vc} is the rate of RuBP carboxylation, \code{Rp} is the rate of
  carbon loss due to photorespiration, and \code{RL} is the rate of carbon loss
  due to non-photorespiratory respiration (also known as the rate of day
  respiration, the rate of mitochondrial respiration, or the rate of respiration
  in the light). Because RuBP carboxylation and photorespiration both occur due
  to Rubisco activity, these rates are actually proportional to each other:

  \code{Rp = Vc * Gamma_star / Cc},

  where \code{Cc} is the CO2 concentration in the chloroplast (where Rubisco is
  located) and \code{Gamma_star} will be discussed below. Using this expression,
  the net CO2 assimilation rate can be written as

  \code{An = Vc * (1 - Gamma_star / Cc) - RL}.

  When \code{Cc} is equal to \code{Gamma_star}, the net assimilation rate is
  equal to \code{-RL}. For this reason, \code{Gamma_star} is usually referred to
  as the CO2 compensation point in the absence of mitochondrial respiration.

  In general, \code{Cc} is related to the intercellular CO2 concentration
  \code{Ci} according to

  \code{Ci = Cc + An / gmc},

  where \code{gmc} is the mesophyll conductance to CO2 diffusion. When \code{Cc}
  is equal to \code{Gamma_star}, we therefore have
  \code{Ci = Gamma_star - RL / gmc}. This special value of \code{Ci} is referred
  to as \code{Ci_star}, and can be understood as the value of \code{Ci} where
  \code{Cc = Gamma_star} and \code{An = -RL}. Note that the values of
  \code{Gamma_star} and \code{Ci_star} depend on Rubisco properties, mesophyll
  conductance, and the ambient O2 concentration, but not on the incident light
  intensity.

  These observations suggest a method for estimating \code{RL} from a leaf:
  Measure \code{An} vs. \code{Ci} curves at several light intensities, and find
  the value of \code{Ci} where the curves intersect with each other. This will
  be \code{Ci_star}, and the corresponding value of \code{An} will be equal to
  \code{-RL}.

  In practice, it is unlikely that the measured curves will all exactly
  intersect at a single point. A method for dealing with this issue was
  developed in Walker & Ort (2015) and described in more detail in Busch et al.
  (2024). Briefly, a linear fit is first made to each A-Ci curve, enabling the
  calculation of an intercept-slope curve. Then another linear fit is made to
  the intercept-slope curve. The intercept of this fit is equal to \code{-RL}
  and its slope is equal to \code{-Ci_star}.

  Note: it is possible that \code{RL} depends on incident light intensity, an
  issue which complicates the application of the Laisk method. See the
  references for more details.

  \strong{References}:

  \itemize{
    \item Yin, X., Sun, Z., Struik, P. C. & Gu, J. "Evaluating a new method to
          estimate the rate of leaf respiration in the light by analysis of
          combined gas exchange and chlorophyll fluorescence measurements."
          Journal of Experimental Botany 62, 3489–3499 (2011)
          [\doi{10.1093/jxb/err038}].
    \item Walker, B. J. & Ort, D. R. "Improved method for measuring the
          apparent CO2 photocompensation point resolves the impact of multiple
          internal conductances to CO2 to net gas exchange." Plant, Cell &
          Environment 38, 2462–2474 (2015)
          [\doi{10.1111/pce.12562}].
    \item Busch, F. A. et al. "A guide to photosynthetic gas exchang
          measurements: Fundamental principles, best practice and potential
          pitfalls." Plant, Cell & Environment 47, 3344–3364 (2024)
          [\doi{10.1111/pce.14815}].
  }


}

\value{
  This function returns a list with the following named elements:
  \itemize{
    \item \code{first_fit_parameters}: An \code{exdf} object with the slope (and
          its standard error), intercept (and its standard error), R-squared
          value, and p-value for each linear fit of A vs. Ci. These are included
          as the \code{laisk_slope}, \code{laisk_slope_err},
          \code{laisk_intercept}, \code{laisk_intercept_err}, \code{r_squared},
          and \code{p_value} columns.
    \item \code{first_fits}: An \code{exdf} object based on
          \code{replicate_exdf} that also includes the fitted values of
          \code{An} in a new column whose name is \code{a_column_name} followed
          by \code{_fit} (for example, \code{A_fit}). The fits are extrapolated
          to \code{Ci = 0} so they can be visually checked for a common
          intersection point.
    \item \code{second_fit_parameters}: An \code{exdf} object with \code{RL}
          (and its standard error), \code{Ci_Star} (and its standard error) as
          estimated from a linear fit of \code{laisk_intercept} vs.
          \code{laisk_slope}. Also includes the R-squared and p-value of the
          fit.
    \item \code{second_fit_parameters}: An \code{exdf} object based on
          \code{first_fit_parameters} that also includes fitted values of
          \code{laisk_intercept} in the \code{laisk_intercept_fit} column.
  }
  As noted above, the estimated values of \code{RL} and \code{Ci_star} are
  included in the \code{second_fit_parameters} element of the returned list.
}

\examples{
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

# Apply the Laisk method. Note: this is a bad example because these curves were
# measured at the same light intensity, but from different species. Because of
# this, the results are not meaningful.
laisk_results <- fit_laisk(
  licor_file, 20, 150,
  ppfd_column_name = 'species_plot'
)

# Get estimated values
print(laisk_results$second_fit_parameters[, 'RL'])
print(laisk_results$second_fit_parameters[, 'Ci_star'])

# Plot the linear fits of A vs. Ci
plot_laisk_fit(laisk_results, 'instrument', 'first', ppfd_column_name = 'species_plot')

# Plot the linear fits of Laisk intercept vs. Laisk slope
plot_laisk_fit(laisk_results, 'instrument', 'second', ppfd_column_name = 'species_plot')

}

\concept{exdf}
