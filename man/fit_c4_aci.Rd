\name{fit_c4_aci}

\alias{fit_c4_aci}

\title{Fits a C4 assimilation model to an experimental curve}

\description{
  Fits a model to an experimentally measured C4 CO2 response curve using the
  data in the \code{exdf} object along with a few other user-supplied
  parameters. This function can accomodate alternative column names for the
  variables taken from the Licor file in case they change at some point in the
  future. This function also checks the units of each required column and will
  produce an error if any units are incorrect.
}

\usage{
  fit_c4_aci(
    replicate_exdf,
    ca_atmospheric,
    a_column_name = 'A',
    ca_column_name = 'Ca',
    ci_column_name = 'Ci',
    pcm_column_name = 'PCm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    kp_column_name = 'Kp',
    gamma_star_column_name = 'gamma_star',
    ao_column_name = 'ao',
    vcmax_norm_column_name = 'Vcmax_norm',
    vpmax_norm_column_name = 'Vpmax_norm',
    rd_norm_column_name = 'Rd_norm',
    POm = 210000,
    gbs = 0.003,
    Rm_frac = 0.5,
    alpha = 0,
    OPTIM_FUN = optimizer_nmkb(),
    initial_guess_fun = initial_guess_c4_aci(
      gbs = gbs,
      Rm_frac = Rm_frac,
      a_column_name = a_column_name,
      pcm_column_name = pcm_column_name,
      kp_column_name = kp_column_name,
      rd_norm_column_name = rd_norm_column_name,
      vcmax_norm_column_name = vcmax_norm_column_name,
      vpmax_norm_column_name = vpmax_norm_column_name
    ),
    lower = c(0,   0,   0,    0),
    upper = c(100, 100, 1000, 1000),
    fixed = c(NA,  NA,  NA,   1000)
  )
}

\arguments{
  \item{replicate_exdf}{
    An \code{exdf} object representing one CO2 response curve.
  }

  \item{ca_atmospheric}{
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

  \item{ci_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    intercellular CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{pcm_column_name}{
    The name of the column in \code{exdf_obj} that contains the partial pressure
    of CO2 in the mesophyll, expressed in \code{microbar}.
  }

  \item{kc_column_name}{
    The name of the column in \code{exdf_obj} that contains the Michaelis-Menten
    constant for rubisco carboxylation in \code{microbar}.
  }

  \item{ko_column_name}{
    The name of the column in \code{exdf_obj} that contains the Michaelis-Menten
    constant for rubisco oxygenation in \code{mbar}.
  }

  \item{kp_column_name}{
    The name of the column in \code{exdf_obj} that contains the Michaelis-Menten
    constant for PEP carboxylase carboxylation in \code{microbar}.
  }

  \item{gamma_star_column_name}{
    The name of the column in \code{exdf_obj} that contains the dimensionless
    \code{gamma_star} values.
  }

  \item{ao_column_name}{
    The name of the column in \code{exdf_obj} that contains the dimensionless
    ratio of solubility and diffusivity of O2 to CO2.
  }

  \item{vcmax_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Vcmax} values (with units of
    \code{normalized to Vcmax at 25 degrees C}).
  }

  \item{vpmax_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Vpmax} values (with units of
    \code{normalized to Vpmax at 25 degrees C}).
  }

  \item{rd_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Rd} values (with units of \code{normalized to Rd at 25 degrees C}).
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

  \item{alpha}{
    The fraction of photosystem II activity in the bundle sheath.
  }

  \item{OPTIM_FUN}{
    An optimization function that accepts the following input arguments: an
    initial guess, an error function, lower bounds, and upper bounds. It should
    return a list with the following elements: \code{par}, \code{convergence},
    \code{value}, and (optionally) \code{message}. The bounded optimizers from
    the \code{dfoptim} package meet these requirements. The base function
    \code{optim} with the \code{L-BFGS-B} method can also be used.
  }

  \item{initial_guess_fun}{
    A function that estimates values of \code{Rd}, \code{Vcmax}, \code{Vpmax},
    and \code{Vpr}; \code{initial_guess_fun} should have one input argument (an
    \code{exdf} object representing a single response curve and should return a
    vector of the estimated values (at 25 degrees C) in the order above.
  }

  \item{lower}{
    A numeric vector containing lower bounds for the values of \code{Rd},
    \code{Vcmax}, \code{Vpmax}, and \code{Vpr}, in that order.
  }

  \item{upper}{
    A numeric vector containing upper bounds for the values of \code{Rd},
    \code{Vcmax}, \code{Vpmax}, and \code{Vpr}, in that order.
  }

  \item{fixed}{
    A numeric vector specifying fixed values of \code{Rd}, \code{Vcmax},
    \code{Vpmax}, and \code{Vpr}, in that order. Values of \code{NA} indicate
    that the corresponding parameter is not fixed and should be varied during
    the optimization procedure. The default value of this argument fixes
    \code{Vpr} to a high value (so net assimilation is rarely limited by PEP
    carboxylase regeneration) and allows the other parameters to vary. To vary
    all of the parameters, set \code{fixed} to \code{c(NA, NA, NA, NA)}.
  }
}

\details{
  This function calls \code{\link{calculate_c4_assimilation}} to calculate
  values of net assimilation. The user-supplied optimization function is used to
  vary the values of \code{Rd}, \code{Vcmax}, \code{Vpmax}, and \code{Vpr} to
  find ones that best reproduce the experimentally measured values of net
  assimilation.

  The photosynthesis model represented by \code{calculate_c4_assimilation} is
  not smooth in the sense that small changes in the input parameters do not
  necessarily cause changes in its outputs. This is related to the calculation
  of the PEP carboxylase activity \code{Vp}, which is taken to be the minimum of
  \code{Vpr} and \code{Vpc}. For example, if \code{Vpr} is high and \code{Vp =
  Vpc} at all points along the curve, modifying \code{Vpr} by a small amount
  will not change the model's outputs. For this reason, derivative-based
  optimizers tend to struggle when fitting C4 A-Ci curves. Best results are
  obtained using derivative-free methods. It has been found that
  \code{dfoptim::nmkb} is often able to find a good fit.

  By default, this function will use \code{initial_guess_c4_aci} to make rough
  estimates of the parameters to use as an initial guess. However, it is
  possible to also use a fixed initial guess by using something like the
  following for the \code{initial_guess_fun}:
  \code{function(x){c(1, 20, 60, 90)}}, which would set the initial guess to
  \code{Rd = 1}, \code{Vcmax = 20}, \code{Vpmax = 60}, and \code{Vpr = 90}.

  By default, this function fixes \code{Vpr} to a high value (\code{800 micromol
  m^(-2) s^(-1)}) and allows the other parameters to be fitted. This behavior
  is determined by the value of the \code{fixed} input argument, and the value
  of \code{Vpr} specified in this way will override any value returned by
  \code{initial_guess_fun}.

  Once the best-fit parameters have been determined, this function also
  estimates the operating value of `\code{PCm} from the atmospheric CO2
  concentration \code{atmospheric_ca} using
  \code{\link{estimate_operating_point}}, and then uses that value to estimate
  the modeled \code{An} at the operating point via
  \code{\link{calculate_c4_assimilation}}.

  This function assumes that \code{replicate_exdf} represents a single
  C4 A-Ci curve. To fit multiple curves at once, this function is often used
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
                  the \code{Rd_at_25}, \code{Vcmax_at_25}, \code{Vpmax_at_25},
                  and \code{Vpr} columns.
            \item The other outputs from \code{\link{calculate_c4_assimilation}}
                  will be stored in columns with the usual names:
                  \code{Vcmax_tl}, \code{Vpmax_tl}, \code{Rd_tl}, \code{Rm_tl},
                  \code{Vp}, \code{Apc}, \code{Apr}, \code{Ap}, and \code{Ar}.
          }
    \item \code{parameters}: An \code{exdf} object including the identifiers,
          fitting parameters, and convergence information for the A-Ci curve:
          \itemize{
            \item The best-fit values are stored in the \code{Rd_at_25},
                  \code{Vcmax_at_25}, \code{Vpmax_at_25}, and \code{Vpr}
                  columns.
            \item For parameters that depend on leaf temperature, the average
                  leaf-temperature-dependent values are stored in
                  \code{X_tl_avg} columns: \code{Rd_tl_avg},
                  \code{Vcmax_tl_avg}, and \code{Vpmax_tl_avg}.
            \item Information about the operating point is stored in
                  \code{operating_PCm}, \code{operating_Ci},
                  \code{operating_An}, and \code{operating_An_model}.
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

# Fit just one curve from the data set (it is rare to do this)
one_result <- fit_c4_aci(
  licor_file[licor_file[, 'species_plot'] == 'maize - 5', , TRUE]
)

# Fit all curves in the data set (it is more common to do this)
aci_results <- consolidate(by(
  licor_file,
  licor_file[, 'species_plot'],
  fit_c4_aci
))

# View the fitting parameters for each species / plot
col_to_keep <- c(
  'species', 'plot',                                       # identifiers
  'Rd_at_25', 'Vcmax_at_25', 'Vpmax_at_25', 'Vpr',         # parameters scaled to 25 degrees C
  'Rd_tl_avg', 'Vcmax_tl_avg', 'Vpmax_tl_avg',             # average temperature-dependent values
  'operating_Ci', 'operating_An', 'operating_An_model',    # operating point info
  'dof', 'RSS', 'MSE', 'RMSE', 'RSE',                      # residual stats
  'convergence', 'convergence_msg', 'feval', 'optimum_val' # convergence info
)

aci_results$parameters[ , col_to_keep, TRUE]

# View the fits for each species / plot
lattice::xyplot(
  Ar + Apc + A_fit + A ~ Ci | species_plot,
  data = aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlab = paste('Intercellular CO2 concentration [', aci_results$fits$units$Ci, ']'),
  ylab = paste('Assimilation rate [', aci_results$fits$units$A, ']')
)

# View the residuals for each species / plot
lattice::xyplot(
  A_residuals ~ Ci | species_plot,
  data = aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlab = paste('Intercellular CO2 concentration [', aci_results$fits$units$Ci, ']'),
  ylab = paste('Assimilation rate residuals [', aci_results$fits$units$A_residuals, ']')
)
}

\concept{exdf}
