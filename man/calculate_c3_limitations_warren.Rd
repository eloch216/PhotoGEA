\name{calculate_c3_limitations_warren}

\alias{calculate_c3_limitations_warren}

\title{Estimate the relative limiting factors to C3 photosynthesis}

\description{
  Uses the method from Warren et al. (2003) to estimate the relative limitations
  to C3 photosynthesis due to stomatal conductance and mesophyll conductance.
  This function can accomodate alternative column names for the variables taken
  from the data file in case they change at some point in the future. This
  function also checks the units of each required column and will produce an
  error if any units are incorrect.
}

\usage{
  calculate_c3_limitations_warren(
    exdf_obj,
    atp_use = 4.0,
    nadph_use = 8.0,
    curvature_cj = 1.0,
    curvature_cjp = 1.0,
    ca_column_name = 'Ca',
    cc_column_name = 'Cc',
    ci_column_name = 'Ci',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    oxygen_column_name = 'oxygen',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    hard_constraints = 0,
    \dots
  )
}

\arguments{
  \item{exdf_obj}{
    An \code{exdf} object representing gas exchange data. Typically this should
    be an \code{exdf} object returned from \code{\link{fit_c3_aci}}; it will be
    expected to have columns for \code{alpha_g}, \code{Gamma_star},
    \code{J_at_25}, \code{RL_at_25}, \code{Tp}, and \code{Vcmax_at_25}.
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

  \item{ca_column_name}{
    The name of the column in \code{exdf_obj} that contains the ambient
    CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{cc_column_name}{
    The name of the column in \code{exdf_obj} that contains the chloroplastic
    CO2 concentration in \code{micromol mol^(-1)}. Typically these are
    calculated using \code{\link{apply_gm}}.
  }

  \item{ci_column_name}{
    The name of the column in \code{exdf_obj} that contains the intercellular
    CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{j_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{J} values (with units of \code{normalized to J at 25 degrees C}).
    Typically these are the leaf-temperature dependent values calculated using
    \code{\link{calculate_temperature_response}}.
  }

  \item{kc_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    Michaelis-Menten constant for rubisco carboxylation in
    \code{micromol mol^(-1)}. Typically these are the leaf-temperature dependent
    values calculated using
    \code{\link{calculate_temperature_response}}.
  }

  \item{ko_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    Michaelis-Menten constant for rubisco oxygenation in
    \code{mmol mol^(-1)}. Typically these are the leaf-temperature dependent
    values calculated using
    \code{\link{calculate_temperature_response}}.
  }

  \item{oxygen_column_name}{
    The name of the column in \code{exdf_obj} that contains the concentration
    of O2 in the ambient air, expressed as a percentage (commonly 21\% or 2\%);
    the units must be \code{percent}.
  }

  \item{rl_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{RL} values (with units of \code{normalized to RL at 25 degrees C}).
  }

  \item{total_pressure_column_name}{
    The name of the column in \code{exdf_obj} that contains the total pressure
    in \code{bar}. Typically this is calculated using
    \code{\link{calculate_total_pressure}}.
  }

  \item{vcmax_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{Vcmax} values (with units of
    \code{normalized to Vcmax at 25 degrees C}).
  }

  \item{hard_constraints}{
    To be passed to \code{\link{calculate_c3_assimilation}}; see that function
    for more details.
  }

  \item{\dots}{
    Additional arguments to be passed to
    \code{\link{calculate_c3_assimilation}}.
  }
}

\details{
  When analyzing or interpreting C3 gas exchange data, it is often useful to
  estimate the relative limitations to assimilation that are due to stomatal
  conductance or mesophyll conductance. This can be done using a framework first
  introduced by Warren et al. (2003). In this framework, the relative limitation
  due to stomatal conductance (\code{ls}) is

  \code{ls = (An_inf_gsc - A_modeled) / An_inf_gsc}

  and the relative limitation due to mesophyll conductance (\code{lm}) is

  \code{lm = (An_inf_gmc - A_modeled) / An_inf_gmc}. These are equations 10 and
  11 in Warren et al. (2003).

  In these equations \code{A_modeled} is the net assimilation rate calculated
  using the Farquhar-von-Caemmerer-Berry (FvCB) model at the measured value of
  the chloroplast CO2 concentration (\code{Cc}). The other two assimilation
  rates (\code{An_inf_gsc} and \code{An_inf_gmc}) are also calculated using the
  FvCB model, but under different assumptions: \code{An_inf_gsc} assumes that
  stomatal conductance is infinite while mesophyll conductance is as measured,
  while \code{An_inf_gmc} assumes that mesophyll conductance is infinite while
  stomatal conductance is as measured.

  In other words, \code{ls} expresses the observed assimilation rate as a
  fractional decrease relative to a hypothetical plant with infinite stomatal
  conductance, while \code{lm} expresses the observed assimilation rate as a
  fractional decrease relative to a hypothetical plant with infinite mesophyll
  conductance.

  For example, if \code{lm = 0.4}, this means that the observed assimilation
  rate is 40 percet lower than a hypothetical plant with infinite mesophyll
  conductance. If mesophyll conductance were to increase (all else remaining
  the same), then \code{lm} would decrease. This is not the case with other
  estimations of limiting factors, such as the one used in
  \code{\link{calculate_c3_limitations_grassi}}. (See Leverett & Kromdijk for
  more details.)

  To actually calculate \code{An_inf_gsc} and \code{An_inf_gmc}, it is first
  necessary to estimate the corresponding values of \code{Cc} that would occur
  with infinite stomatal or mesophyll conductance. This can be done with a 1D
  diffusion equation expressed using drawdown values:

  \code{Cc = Ca - drawdown_cs - drawdown_cm},

  where \code{drawdown_cs = Ca - Ci} is the drawdown of CO2 across the stomata
  (assuming infinite boundary layer conductance) and
  \code{drawdown_cm = Ci - Cc} is the drawdown of CO2 across the mesophyll. If
  one conductance is infinite, the corresponding drawdown becomes zero. Thus, we
  have:

  \code{Cc_inf_gsc = Ca - 0 - (Ci - Cc) = Ca - Ci + Cc}

  and

  \code{Cc_inf_gmc = Ca - (Ca - Ci) - 0 = Ci},

  where \code{Cc_inf_gsc} is the value of \code{Cc} that would occur with
  infinite stomatal conductance and the measured mesophyll conductance, and
  \code{Cc_inf_gmc} is the value of \code{Cc} that would occur with infinite
  mesophyll conductance and the measured stomatal conductance.

  Once values of \code{Cc}, \code{Cc_inf_gsc}, and \code{Cc_inf_gmc}, the
  corresponding assimilation rates are calculated using
  \code{\link{calculate_c3_assimilation}}, and then the limitation factors are
  calculated as described above.

  References:

  Warren, C. R. et al. "Transfer conductance in second growth Douglas-fir
  (Pseudotsuga menziesii (Mirb.)Franco) canopies." Plant, Cell & Environment 26,
  1215–1227 (2003) [\doi{10.1046/j.1365-3040.2003.01044.x}].

  Leverett, A. & Kromdijk, J. "The long and tortuous path towards improving
  photosynthesis by engineering elevated mesophyll conductance."
  [\doi{10.22541/au.170016201.13513761/v1}].

}

\value{
  This function returns an \code{exdf} object based on \code{exdf_obj} but with
  several new columns representing the quantities discussed above:
  \code{Cc_inf_gsc}, \code{Cc_inf_gmc}, \code{An_inf_gsc}, \code{An_inf_gmc},
  \code{ls_warren}, and \code{lm_warren}.
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

# Specify mesophyll conductance
licor_file <- set_variable(
  licor_file,
  'gmc', 'mol m^(-2) s^(-1) bar^(-1)', value = 0.5
)

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Calculate additional gas properties
licor_file <- calculate_gas_properties(licor_file)

# Calculate Cc
licor_file <- apply_gm(licor_file)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_file <- calculate_temperature_response(licor_file, c3_temperature_param_bernacchi)

# Fit all curves in the data set. Here we use a faster optimizer than the
# default one to ensure the example runs quickly.
aci_results <- consolidate(by(
  licor_file,
  licor_file[, 'species_plot'],
  fit_c3_aci,
  Ca_atmospheric = 420,
  OPTIM_FUN = optimizer_nmkb(1e-7)
))

# Get a subset of fitting results corresponding to the first measured point
# in each curve (where CO2_r_sp = 400 ppm)
aci_fit_subset <- aci_results$fits[aci_results$fits[, 'CO2_r_sp'] == 400, , TRUE]

# Calculate limiting factors
aci_fit_subset <- calculate_c3_limitations_warren(aci_fit_subset)

# View the limiting factors for each species / plot
col_to_keep <- c(
  'species', 'plot',       # identifiers
  'ls_warren', 'lm_warren' # limitation info
)

aci_fit_subset[ , col_to_keep, TRUE]
}

\concept{exdf}
