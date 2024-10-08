\name{initial_guess_c3_aci}

\alias{initial_guess_c3_aci}

\title{Make an initial guess of FvCB model parameter values for one curve}

\description{
  Creates a function that makes an initial guess of FvCB model parameter values
  for one curve. This function is used internally by \code{\link{fit_c3_aci}}.

  Values estimated by this guessing function should be considered inaccurate,
  and should always be improved upon by an optimizer.
}

\usage{
  initial_guess_c3_aci(
    alpha_g,
    alpha_old,
    alpha_s,
    Gamma_star,
    cc_threshold_rd = 100,
    atp_use = 4.0,
    nadph_use = 8.0,
    a_column_name = 'A',
    cc_column_name = 'Cc',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    oxygen_column_name = 'oxygen',
    rl_norm_column_name = 'RL_norm',
    vcmax_norm_column_name = 'Vcmax_norm'
  )
}

\arguments{
  \item{alpha_g}{
    A dimensionless parameter where \code{0 <= alpha_g <= 1}, representing the
    proportion of glycolate carbon taken out of the photorespiratory pathway as
    glycine. \code{alpha_g} is often assumed to be 0. If \code{alpha_g} is not a
    number, then there must be a column in \code{rc_exdf} called \code{alpha_g}
    with appropriate units. A numeric value supplied here will overwrite the
    values in the \code{alpha_g} column of \code{rc_exdf} if it exists.
  }

  \item{alpha_old}{
    A dimensionless parameter where \code{0 <= alpha_old <= 1}, representing the
    fraction of remaining glycolate carbon not returned to the chloroplast
    after accounting for carbon released as CO2. \code{alpha_old} is often
    assumed to be 0. If \code{alpha_old} is not a number, then there must be a
    column in \code{rc_exdf} called \code{alpha_old} with appropriate units.
    A numeric value supplied here will overwrite the values in the
    \code{alpha_old} column of \code{rc_exdf} if it exists.
  }

  \item{alpha_s}{
    A dimensionless parameter where \code{0 <= alpha_s <= 0.75 * (1 - alpha_g)}
    representing the proportion of glycolate carbon taken out of the
    photorespiratory pathway as serine. \code{alpha_s} is often assumed to be 0.
    If \code{alpha_s} is not a number, then there must be a column in
    \code{rc_exdf} called \code{alpha_s} with appropriate units. A numeric
    value supplied here will overwrite the values in the \code{alpha_s} column
    of \code{rc_exdf} if it exists.
  }

  \item{Gamma_star}{
    The CO2 compensation point in the absence of day respiration, expressed in
    \code{micromol mol^(-1)}. If \code{Gamma_star} is not a number, then there
    must be a column in \code{rc_exdf} called \code{Gamma_star} with
    appropriate units. A numeric value supplied here will overwrite the values
    in the \code{Gamma_star} column of \code{rc_exdf} if it exists.
  }

  \item{cc_threshold_rd}{
    An upper cutoff value for the chloroplast CO2 concentration in
    \code{micromol mol^(-1)} to be used when estimating \code{RL}.
  }

  \item{atp_use}{
    The number of ATP molecules used per C3 cycle.
  }

  \item{nadph_use}{
    The number of NADPH molecules used per C3 cycle.
  }

  \item{a_column_name}{
    The name of the column in \code{rc_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{cc_column_name}{
    The name of the column in \code{rc_exdf} that contains the
    chloroplastic CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{j_norm_column_name}{
    The name of the column in \code{rc_exdf} that contains the normalized
    \code{J} values (with units of \code{normalized to J at 25 degrees C}).
  }

  \item{kc_column_name}{
    The name of the column in \code{rc_exdf} that contains the
    Michaelis-Menten constant for rubisco carboxylation in
    \code{micromol mol^(-1)}.
  }

  \item{ko_column_name}{
    The name of the column in \code{rc_exdf} that contains the
    Michaelis-Menten constant for rubisco oxygenation in
    \code{mmol mol^(-1)}.
  }

  \item{oxygen_column_name}{
    The name of the column in \code{rc_exdf} that contains the concentration
    of O2 in the ambient air, expressed as a percentage (commonly 21\% or 2\%);
    the units must be \code{percent}.
  }

  \item{rl_norm_column_name}{
    The name of the column in \code{rc_exdf} that contains the normalized
    \code{RL} values (with units of \code{normalized to RL at 25 degrees C}).
  }

  \item{vcmax_norm_column_name}{
    The name of the column in \code{rc_exdf} that contains the normalized
    \code{Vcmax} values (with units of
    \code{normalized to Vcmax at 25 degrees C}).
  }
}

\details{
  Here we estimate values of \code{J_at_25}, \code{RL_at_25}, \code{Tp}, and
  \code{Vcmax_at_25} from a measured C3 CO2 response curve. It is difficult to
  estimate values of \code{alpha_g}, \code{alpha_old}, \code{alpha_s}, and
  \code{Gamma_star} from a curve, so they must be supplied beforehand. For more
  information about these parameters, see the documentation for
  \code{\link{calculate_c3_assimilation}}.

  \itemize{
    \item \strong{Estimating RL}: Regardless of which process is limiting at low
          \code{Cc}, it is always true that \code{An = -RL} when
          \code{Cc = Gamma_star_ag}. Here we make a linear fit of the measured
          \code{An} vs. \code{Cc} values, and evaluate it at at
          \code{Cc = Gamma_star_ag} to estimate \code{RL}. If the linear fit
          predicts a negative value for \code{RL}, we use a typical value
          instead (0.5 \code{micromol m^(-2) s^(-1)}).
    \item \strong{Estimating Vc}: Once an estimate for \code{RL} has been found,
          the RuBP carboxylation rate \code{Vc} can be estimated using
          \code{Vc = (An + RL) / (1 - Gamma_star_ag / Cc)}. This is useful for
          the remaining parameter estimates.
    \item \strong{Estimating Vcmax}: An estimate for \code{Vcmax} can be
          obtained by solving the equation for \code{Wc} for \code{Vcmax}, and
          evaluating it with \code{Wc = Vc} as estimated above. In the
          rubisco-limited part of the curve, \code{Vc = Wc} and the estimated
          values of \code{Vcmax} should be reasonable. In other parts of the
          curve, \code{Wc} is not the limiting rate, so \code{Vc < Wc}.
          Consequently, the estimated values of \code{Vcmax} in these parts of
          the curve will be smaller. So, to make an overall estimate, we choose
          the the largest estimated \code{Vcmax} value.
    \item \strong{Estimating J and Tp}: Estimates for these parameters can be
          made using the equations for \code{Wj} and \code{Wp}, similar to the
          approach followed for \code{Vcmax}.
  }

  For the parameter values estimated above, the values of \code{RL_norm},
  \code{Vcmax_norm}, and \code{J_norm} are used to convert the values at leaf
  temperature to the values at 25 degrees C.
}

\value{
  A function with one input argument \code{rc_exdf}, which should be an
  \code{exdf} object representing one C3 CO2 response curve. The return value of
  this function will be a numeric vector with eight elements, representing the
  values of \code{alpha_g}, \code{alpha_old}, \code{alpha_s}, \code{Gamma_star},
  \code{J_at_25}, \code{RL_at_25}, \code{Tp}, and \code{Vcmax_at_25} (in that
  order).
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

# Create the guessing function; here we set all alpha values to 0 and use the
# temperature-dependent values of Gamma_star calculated above
guessing_func <-
  initial_guess_c3_aci(alpha_g = 0, alpha_old = 0, alpha_s = 0, Gamma_star = '')

# Apply it and see the initial guesses for each curve
str(by(licor_file, licor_file[, 'species_plot'], guessing_func))

# Calculate simulated A-Ci curves based on the guesses and compare them to the
# actual data
calculated_aci <- do.call(rbind, by(
  licor_file,
  licor_file[, 'species_plot'],
  function(x) {
    param <- guessing_func(x)
    calculate_c3_assimilation(
      x,
      param[1], param[2], param[3], param[4], param[5], param[6], param[7], param[8]
    )
  }
))

lattice::xyplot(
  Ac + Aj + Ap + An + A ~ Ci | species_plot,
  data = cbind(licor_file, calculated_aci)$main_data,
  type = 'b',
  auto = TRUE,
  grid = TRUE
)
}

\concept{exdf}
