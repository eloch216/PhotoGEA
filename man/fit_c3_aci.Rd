\name{fit_c3_aci}

\alias{fit_c3_aci}

\title{Fits a C3 assimilation model to an experimental curve}

\description{
  Fits a model to an experimentally measured C3 CO2 response curve using the
  data in the \code{exdf} object along with a few other user-supplied
  parameters. This function can accomodate alternative column names for the
  variables taken from the Licor file in case they change at some point in the
  future. This function also checks the units of each required column and will
  produce an error if any units are incorrect.
}

\usage{
  fit_c3_aci(
    replicate_exdf,
    a_column_name = 'A',
    cc_column_name = 'Cc',
    total_pressure_column_name = 'total_pressure',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    gamma_star_column_name = 'Gamma_star',
    vcmax_norm_column_name = 'Vcmax_norm',
    rd_norm_column_name = 'Rd_norm',
    j_norm_column_name = 'J_norm',
    POc = 210000,
    atp_use = 4.0,
    nadph_use = 8.0,
    OPTIM_FUN = default_optimizer(),
    initial_guess_fun = initial_guess_c3_aci(
        Oc = POc,
        atp_use = atp_use,
        nadph_use = nadph_use,
        a_column_name = a_column_name,
        cc_column_name = cc_column_name,
        kc_column_name = kc_column_name,
        ko_column_name = ko_column_name,
        gamma_star_column_name = gamma_star_column_name,
        vcmax_norm_column_name = vcmax_norm_column_name,
        rd_norm_column_name = rd_norm_column_name,
        j_norm_column_name = j_norm_column_name
    ),
    lower = c(0,  0,    0,   0),
    upper = c(40, 1000, 100, 1000),
    fixed = c(40, NA,   NA,  NA),
    min_aj_cutoff = NA,
    max_aj_cutoff = NA,
    curvature = 0.99
  )
}

\arguments{
  \item{replicate_exdf}{
    An \code{exdf} object representing one CO2 response curve.
  }

  \item{a_column_name}{
    The name of the column in \code{replicate_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{cc_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    chloroplastic CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{total_pressure_column_name}{
    The name of the column in \code{replicate_exdf} that contains the total
    pressure in \code{bar}.
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

  \item{gamma_star_column_name}{
    The name of the column in \code{replicate_exdf} that contains the
    \code{Gamma_star} values in \code{micromol mol^(-1)}.
  }

  \item{vcmax_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{Vcmax} values (with units of
    \code{normalized to Vcmax at 25 degrees C}).
  }

  \item{rd_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{Rd} values (with units of \code{normalized to Rd at 25 degrees C}).
  }

  \item{j_norm_column_name}{
    The name of the column in \code{replicate_exdf} that contains the normalized
    \code{J} values (with units of \code{normalized to J at 25 degrees C}).
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

  \item{OPTIM_FUN}{
    An optimization function that accepts the following input arguments: an
    initial guess, an error function, lower bounds, and upper bounds. It should
    return a list with the following elements: \code{par}, \code{convergence},
    \code{value}, and (optionally) \code{message}. The bounded optimizers from
    the \code{dfoptim} package meet these requirements. The base function
    \code{optim} with the \code{L-BFGS-B} method can also be used.
  }

  \item{initial_guess_fun}{
    A function that estimates values of \code{TPU}, \code{J}, \code{Rd}, and
    \code{Vcmax}; \code{initial_guess_fun} should have one input argument (an
    \code{exdf} object representing a single response curve and should return a
    vector of the estimated values (at 25 degrees C) in the order above.
  }

  \item{lower}{
    A numeric vector containing lower bounds for the values of \code{TPU},
    \code{J}, \code{Rd}, and \code{Vcmax}, in that order.
  }

  \item{upper}{
    A numeric vector containing upper bounds for the values of \code{TPU},
    \code{J}, \code{Rd}, and \code{Vcmax}, in that order.
  }

  \item{fixed}{
    A numeric vector specifying fixed values of \code{TPU}, \code{J}, \code{Rd},
    and \code{Vcmax}, in that order. Values of \code{NA} indicate that the
    corresponding parameter is not fixed and should be varied during the
    optimization procedure. The default value of this argument fixes \code{TPU}
    to its upper limit (so net assimilation is never limited by phosphate
    utilization) and allows the other parameters to vary. To vary all of the
    parameters, set \code{fixed} to \code{c(NA, NA, NA, NA)}.
  }

  \item{min_aj_cutoff}{
    The minimum value of \code{Cc} (in ppm) where \code{Aj} is allowed to become
    the overall rate-limiting factor. If \code{min_aj_cutoff} is set to
    \code{NA}, this restriction will not be applied.
  }

  \item{max_aj_cutoff}{
    The maximum value of \code{Cc} (in ppm) where \code{Aj} is allowed to be
    smaller than \code{Ac}. If \code{max_aj_cutoff} is set to \code{NA}, this
    restriction will not be applied.
  }

  \item{curvature}{
    A dimensionless quadratic curvature parameter to be passed to
    \code{\link{calculate_c3_assimilation}} as the input argument with the same
    name.
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

  Sometimes the optimizer will choose biologically unreasonable parameter values
  that nevertheless produce good fits to the supplied assimilation values. A
  common problem is that the fit result may not follow the expected progression
  of \code{Ac-limited -> Aj-limited -> Ap-limited} assimilation as \code{Cc}
  increases. In this case, the optional \code{min_aj_cutoff} and
  \code{max_aj_cutoff} can be used to constrain the range of \code{Cc} values
  (in ppm) where \code{Aj} is allowed to be the overall rate limiting factor.
  See the _Analyzing C3 A-Ci Curves_ vignette for an example of how these
  arguments can be used to improve the quality of a fit.

  By default, this function will use \code{initial_guess_c3_aci} to make rough
  estimates of the parameters to use as an initial guess. However, it is
  possible to also use a fixed initial guess by using something like the
  following for the \code{initial_guess_fun}:
  \code{function(x){c(10, 100, 0.5, 90)}}, which would set the initial guess to
  \code{TPU = 10}, \code{J = 100}, \code{Rd = 0.5}, and \code{Vcmax = 90}.

  By default, this function fixes \code{TPU} to a high value (\code{40 micromol
  m^(-2) s^(-1)}) and allows the other parameters to be fitted. This behavior
  is determined by the value of the \code{fixed} input argument, and the value
  of \code{TPU} specified in this way will override any value returned by
  \code{initial_guess_fun}.

  This function assumes that \code{replicate_exdf} represents a single
  C3 A-Ci curve. To fit multiple curves at once, this function is often used
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
            \item The other outputs from \code{\link{calculate_c3_assimilation}}
                  will be stored in columns with the usual names: \code{TPU},
                  \code{Vcmax_tl}, \code{Rd_tl}, \code{J_tl}, \code{Ac},
                  \code{Aj}, and \code{Ap}.
          }
    \item \code{parameters}: An \code{exdf} object including the identifiers,
          fitting parameters, and convergence information for the A-Ci curve:
          \itemize{
            \item The best-fit values are stored in the \code{TPU},
                  \code{J_at_25}, \code{Rd_at_25}, and \code{Vcmax_at_25}
                  columns.
            \item For parameters that depend on leaf temperature, the average
                  leaf-temperature-dependent values are stored in
                  \code{X_tl_avg} columns: \code{J_tl_avg}, \code{Rd_tl_avg},
                  and \code{Vcmax_tl_avg}.
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
  system.file('extdata', 'c3_aci_1.xlsx', package = 'PhotoGEA')
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

# Fit just one curve from the data set (it is rare to do this)
one_result <- fit_c3_aci(
  licor_file[licor_file[, 'species_plot'] == 'tobacco - 1', , TRUE]
)

# Fit all curves in the data set (it is more common to do this)
aci_results <- consolidate(by(
  licor_file,
  licor_file[, 'species_plot'],
  fit_c3_aci
))

# View the fitting parameters for each species / plot
col_to_keep <- c(
  'species', 'plot',                                       # identifiers
  'TPU', 'J_at_25', 'Rd_at_25', 'Vcmax_at_25',             # parameters scaled to 25 degrees C
  'J_tl_avg', 'Rd_tl_avg', 'Vcmax_tl_avg',                 # average temperature-dependent values
  'dof', 'RSS', 'MSE', 'RMSE', 'RSE',                      # residual stats
  'convergence', 'convergence_msg', 'feval', 'optimum_val' # convergence info
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

# View the residuals for each species / plot
lattice::xyplot(
  A_residuals ~ Ci | species_plot,
  data = aci_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE,
  grid = TRUE,
  xlab = paste0('Intercellular CO2 concentration (', aci_results$fits$units$Ci, ')'),
  ylab = paste0('Assimilation rate residuals (', aci_results$fits$units$A_residuals, ')')
)
}

\concept{exdf}
