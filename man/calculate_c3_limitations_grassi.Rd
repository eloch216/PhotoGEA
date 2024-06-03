\name{calculate_c3_limitations_grassi}

\alias{calculate_c3_limitations_grassi}

\title{Estimate the relative limiting factors to C3 photosynthesis}

\description{
  Uses the method from Grassi & Magnani (2005) to estimate the relative
  limitations to C3 photosynthesis due to stomatal conductance, mesophyll
  conductance, and biochemistry. This function can accomodate alternative
  column names for the variables taken from the data file in case they change
  at some point in the future. This function also checks the units of each
  required column and will produce an error if any units are incorrect.
}

\usage{
  calculate_c3_limitations_grassi(
    exdf_obj,
    atp_use = 4.0,
    nadph_use = 8.0,
    cc_column_name = 'Cc',
    gamma_star_column_name = 'Gamma_star',
    gmc_column_name = 'gmc',
    gsc_column_name = 'gsc',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    oxygen_column_name = 'oxygen',
    total_pressure_column_name = 'total_pressure',
    vcmax_column_name = 'Vcmax_tl',
    j_column_name = NULL
  )
}

\arguments{
  \item{exdf_obj}{
    An \code{exdf} object representing gas exchange data.
  }

  \item{atp_use}{
    The number of ATP molecules used per C3 cycle.
  }

  \item{nadph_use}{
    The number of NADPH molecules used per C3 cycle.
  }

  \item{cc_column_name}{
    The name of the column in \code{exdf_obj} that contains the chloroplastic
    CO2 concentration in \code{micromol mol^(-1)}. Typically these are
    calculated using \code{\link{apply_gm}}.
  }

  \item{gamma_star_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    \code{Gamma_star} values in \code{micromol mol^(-1)}. Typically these are
    the leaf-temperature dependent values calculated using
    \code{\link{calculate_arrhenius}}.
  }

  \item{gmc_column_name}{
    The name of the column in \code{exdf_obj} that contains the mesophyll
    conductance to CO2 in \code{mol m^(-2) s^(-1) bar^(-1)}.
  }

  \item{gsc_column_name}{
    The name of the column in \code{exdf_obj} that contains the stomatal
    conductance to CO2 in \code{mol m^(-2) s^(-1)}. Typically this column is
    calculated using \code{\link{calculate_gas_properties}}.
  }

  \item{kc_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    Michaelis-Menten constant for rubisco carboxylation in
    \code{micromol mol^(-1)}. Typically these are the leaf-temperature dependent
    values calculated using \code{\link{calculate_arrhenius}}.
  }

  \item{ko_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    Michaelis-Menten constant for rubisco oxygenation in
    \code{mmol mol^(-1)}. Typically these are the leaf-temperature dependent
    values calculated using \code{\link{calculate_arrhenius}}.
  }

  \item{oxygen_column_name}{
    The name of the column in \code{exdf_obj} that contains the concentration
    of O2 in the ambient air, expressed as a percentage (commonly 21\% or 2\%);
    the units must be \code{percent}.
  }

  \item{total_pressure_column_name}{
    The name of the column in \code{exdf_obj} that contains the total pressure
    in \code{bar}. Typically this is calculated using
    \code{\link{calculate_total_pressure}}.
  }

  \item{vcmax_column_name}{
    The name of the column in \code{exdf_obj} that contains values of the
    maximum Rubisco carboxylation rate (\code{Vcmax}) in
    \code{micromol m^(-2) s^(-1)}. Typically these are the leaf-temperature
    adjusted values that are automatically calculated by
    \code{\link{fit_c3_aci}}.
  }

  \item{j_column_name}{
    The name of the column in \code{exdf_obj} that contains values of the
    RuBP regeneration rate (\code{J}) in \code{micromol m^(-2) s^(-1)}.
    Typically these are the leaf-temperature adjusted values that are
    automatically calculated by \code{\link{fit_c3_aci}}.
  }
}

\details{
  When analyzing or interpreting C3 gas exchange data, it is often useful to
  estimate the relative limitations to assimilation that are due to stomatal
  conductance, mesophyll conductance, and biochemistry. This can be done using
  a framework first introduced by Grassi & Magnani (2005). In this framework,
  the relative limitation due to stomatal conductance (\code{ls}) is

  \code{ls = [(g_t / g_sc) * (dAdC)] / [g_t + dAdC]},

  the relative limitation due to mesophyll conductance (\code{lm}) is

  \code{lm = [(g_t / g_mc) * (dAdC)] / [g_t + dAdC]},

  and the relative limitation due to biochemistry (\code{lb}) is

  \code{ln = [g_t] / [g_t + dAdC]},

  where \code{g_sc} is the stomatal conductance to CO2, \code{g_mc} is the
  mesophyll conductance to CO2, \code{gt = 1 / (1 / g_mc + 1 / g_sc)} is the
  total conductance to CO2, and \code{dAdC} is the partial derivative of the
  net CO2 assimilation rate (\code{An}) with respect to the chloroplast CO2
  concentration (\code{Cc}). These can be found in Equation 7 from Grassi &
  Magnani (2005).

  These equations were derived by assuming that CO2 assimilation is limited by
  Rubisco activity; in other words, that the net CO2 assimilation rate is given
  by

  \code{Ac = Vcmax * (Cc - Gamma_star) / (Cc + Km) - Rd},

  where \code{Vcmax} is the maximum Rubisco carboxylation rate,
  \code{Gamma_star} is the CO2 compensation point in the absence of day
  respiration, \code{Rd} is the day respiration rate, \code{Km} is the effective
  Michaelis-Menten constant for Rubisco carboxylation. In turn, \code{Km} is
  given by \code{Km = Kc * (1 + O / Ko)}, where \code{Kc} is the
  Michaelis-Menten constant for carboxylation, \code{Ko} is the Michaelis-Menten
  constant for oxygenation, and \code{O} is the oxygen concentration in the
  chloroplast.

  Under this assumption, it is possible to analytically determine the partial
  derivative \code{dAdC}:

  \code{dAdC_rubisco = Vcmax * (Gamma_star + Km) / (Cc + Km)^2}

  In this case, the limitation due to "biochemistry" actually refers to
  limitation due to the value of \code{Vcmax}. Note that sometimes this
  derivative is estimated from the initial slope of a measured A-Ci curve rather
  than calculated analytically. (See, for example, Pathare et al. (2023).)
  However, we do not take that approach here. Also note that the value of
  \code{Vcmax} can be estimated using different approaches. For example, Xiong
  (2023) uses single-point gas exchange measurements. When possible, it would be
  better to use an estimate from fitting an entire A-Ci curve, as shown in the
  example below.

  To understand the meaning of these limiting factors, note that simultaneously
  making small fractional increases to \code{g_sc}, \code{g_mc}, and
  \code{Vcmax} will generally cause an associated small fractional increase in
  \code{An}. The limiting factors describe the fraction of the increase in
  \code{An} that can be attributed to each of \code{g_sc}, \code{g_mc}, and
  \code{Vcmax}. For example, \code{ls = 0.2, lm = 0.3, lb = 0.5} would mean that
  20 percent of the increase in \code{An} would be due to an increase in
  stomatal conductance, 30 percent due to an increase in mesophyll conductance,
  and 50 percent due to an increase in \code{Vcmax}. Note that \code{ls},
  \code{lm}, and \code{lb} always add up to 1.

  Thus, when one of the factors is large, changes in the related parameter
  produce relatively larger changes in the assimilation rate. In that case, it
  can be said that that parameter is setting a large limit on the assimilation
  rate. On the other hand, if a factor is small, small changes in the related
  parameter produce relatively small changes in \code{An}, and therefore that
  parameter is not setting a large limit on the assimilation rate.

  It is also possible to calculate \code{dAdC} when assimilation is limited by
  RuBP regeneration. In this case, we have

  \code{Aj = J * (Cc - Gamma_star) / (4 * Cc + 8 * Gamma_star) - Rd},

  where \code{J} is the RuBP regeneration rate, and the limitation due to
  "biochemistry" actually refers to limitation due to the value of \code{J}
  (rather than \code{Vcmax}. The same equations as before can be used to
  calculate the limiting factors (\code{ls}, \code{lm}, \code{lb}), but the
  partial derivative is now given by

  \code{dAdC_j = J * Gamma_star * 12 / (4 * Cc + 8 * Gamma_star)^2}.

  Most users will want the limitations assuming Rubisco-limited assimilation.
  However, if \code{j_column_name} is not \code{NULL}, values of \code{J} will
  be used to calculate the limiting factors assuming RuBP-regeneration-limited
  assimilation. For an example of how these additional factors can be used, see
  Sakoda et al. (2021).

  References:

  Grassi, G. & Magnani, F. "Stomatal, mesophyll conductance and biochemical
  limitations to photosynthesis as affected by drought and leaf ontogeny in ash
  and oak trees." Plant, Cell & Environment 28, 834–849 (2005)
  [\doi{10.1111/j.1365-3040.2005.01333.x}].

  Pathare, V. S. et al. "Altered cell wall hydroxycinnamate composition impacts
  leaf- and canopy-level CO2 uptake and water use in rice." Plant Physiology
  kiad428 (2023) [\doi{10.1093/plphys/kiad428}].

  Xiong, D. "Leaf anatomy does not explain the large variability of mesophyll
  conductance across C3 crop species." The Plant Journal 113, 1035–1048 (2023)
  [\doi{10.1111/tpj.16098}].

  Sakoda, K., Yamori, W., Groszmann, M. & Evans, J. R. "Stomatal, mesophyll
  conductance, and biochemical limitations to photosynthesis during induction."
  Plant Physiology 185, 146–160 (2021)
  [\doi{10.1093/plphys/kiaa011}].

}

\value{
  This function returns an \code{exdf} object based on \code{exdf_obj} but with
  several new columns representing the partial derivatives and limiting factors
  discussed above: \code{dAdC_rubisco}, \code{ls_rubisco_grassi}, \code{lm_rubisco_grassi},
  and \code{lb_rubisco_grassi}. If \code{j_column_name} is not \code{NULL}, the output
  will also include \code{dAdC_j}, \code{ls_j_grassi}, \code{lm_j_grassi}, and \code{lb_j_grassi}.
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
licor_file <- calculate_arrhenius(licor_file, c3_arrhenius_bernacchi)

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
aci_fit_subset <- calculate_c3_limitations_grassi(aci_fit_subset)

# View the limiting factors for each species / plot
col_to_keep <- c(
  'species', 'plot',                                            # identifiers
  'ls_rubisco_grassi', 'lm_rubisco_grassi', 'lb_rubisco_grassi' # limitation info
)

aci_fit_subset[ , col_to_keep, TRUE]

# Display as a bar chart
lattice::barchart(
  ls_rubisco_grassi + lm_rubisco_grassi + lb_rubisco_grassi ~ species_plot,
  data = aci_fit_subset$main,
  stack = TRUE,
  auto = TRUE,
  ylab = 'Factors limiting assimilation'
)

}

\concept{exdf}
