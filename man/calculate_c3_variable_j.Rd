\name{calculate_c3_variable_j}

\alias{calculate_c3_variable_j}

\title{Calculate C3 variable J}

\description{
  Calculates values of mesophyll conductance and chloroplast CO2 concentration
  using the "variable J" equation, as originally described in Harley et al.
  (1992) and modified in Moualeu-Ngangue, Chen, & Stutzel (2016). This function
  can accomodate alternative colum names for the variables taken from Licor
  files in case they change at some point in the future. This function also
  checks the units of each required column and will produce an error if any
  units are incorrect.
}

\usage{
  calculate_c3_variable_j(
    exdf_obj,
    alpha_g,
    alpha_s,
    alpha_t,
    Gamma_star,
    RL_at_25,
    tau,
    atp_use = 4.0,
    nadph_use = 8.0,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    phips2_column_name = 'PhiPS2',
    qin_column_name = 'Qin',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    hard_constraints = 0,
    perform_checks = TRUE,
    return_exdf = TRUE
  )
}

\arguments{
  \item{exdf_obj}{An \code{exdf} object.}

  \item{alpha_g}{
    A dimensionless parameter where \code{0 <= alpha_g <= 1}, representing the
    proportion of glycolate carbon taken out of the photorespiratory pathway as
    glycine. \code{alpha_g} is often assumed to be 0. If \code{alpha_g} is not a
    number, then there must be a column in \code{exdf_obj} called \code{alpha_g}
    with appropriate units. A numeric value supplied here will overwrite the
    values in the \code{alpha_g} column of \code{exdf_obj} if it exists.
  }

  \item{alpha_s}{
    A dimensionless parameter where \code{0 <= alpha_s <= 0.75 * (1 - alpha_g)}
    representing the proportion of glycolate carbon taken out of the
    photorespiratory pathway as serine. \code{alpha_s} is often assumed to be 0.
    If \code{alpha_s} is not a number, then there must be a column in
    \code{exdf_obj} called \code{alpha_s} with appropriate units. A numeric
    value supplied here will overwrite the values in the \code{alpha_s} column
    of \code{exdf_obj} if it exists.
  }

  \item{alpha_t}{
    A dimensionless parameter where \code{0 <= alpha_t <= 1} representing the
    proportion of glycolate carbon taken out of the photorespiratory pathway as
    CH2-THF. \code{alpha_t} is often assumed to be 0. If \code{alpha_t} is not a
    number, then there must be a column in \code{exdf_obj} called \code{alpha_t}
    with appropriate units. A numeric value supplied here will overwrite the
    values in the \code{alpha_t} column of \code{exdf_obj} if it exists.
  }

  \item{Gamma_star}{
    The CO2 compensation point in the absence of day respiration, expressed in
    \code{micromol mol^(-1)}. If \code{Gamma_star} is not a number, then there
    must be a column in \code{exdf_obj} called \code{Gamma_star} with
    appropriate units. A numeric value supplied here will overwrite the values
    in the \code{Gamma_star} column of \code{exdf_obj} if it exists.
  }

  \item{RL_at_25}{
    The respiration rate at 25 degrees C, expressed in
    \code{micromol m^(-2) s^(-1)}. If \code{RL_at_25} is not a number, then
    there must be a column in \code{exdf_obj} called \code{RL_at_25} with
    appropriate units. A numeric value supplied here will overwrite the values
    in the \code{RL_at_25} column of \code{exdf_obj} if it exists.
  }

  \item{tau}{
    The proportionality factor used to calculate the RuBP regeneration rate from
    chlorophyll fluorescence measurements (dimensionless). If \code{tau} is not
    a number, then there must be a column in \code{exdf_obj} called \code{tau}
    with appropriate units. A numeric value supplied here will overwrite the
    values in the \code{tau} column of \code{exdf_obj} if it exists.
  }

  \item{atp_use}{
    The number of ATP molecules used per C3 cycle.
  }

  \item{nadph_use}{
    The number of NADPH molecules used per C3 cycle.
  }

  \item{a_column_name}{
    The name of the column in \code{exdf_obj} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ci_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    intercellular CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{phips2_column_name}{
    The name of the column in \code{exdf_obj} that contains values of the
    operating efficiency of photosystem II (dimensionless).
  }

  \item{qin_column_name}{
    The name of the column in \code{exdf_obj} that contains values of the
    incident photosynthetically active flux density in
    \code{micromol m^(-2) s^(-1)}.
  }

  \item{rl_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    \code{RL} values (with units of \code{normalized to RL at 25 degrees C}).
  }

  \item{total_pressure_column_name}{
    The name of the column in \code{exdf_obj} that contains the total pressure
    in \code{bar}.
  }

  \item{hard_constraints}{
    An integer numerical value indicating which types of hard constraints to
    place on the values of input parameters; see below for more details.
  }

  \item{perform_checks}{
    A logical value indicating whether to check units for the required columns.
    This should almost always be \code{TRUE}. The option to disable these checks
    is only intended to be used when \code{\link{fit_c3_variable_j}} calls this
    function, since performing these checks many times repeatedly slows down the
    fitting procedure.
  }

  \item{return_exdf}{
    A logical value indicating whether to return an \code{exdf} object. This
    should almost always be \code{TRUE}. The option to return a vector is mainly
    intended to be used when  \code{\link{fit_c3_variable_j}} calls this
    function, since creating an \code{exdf} object to return will slow down the
    fitting procedure.
  }
}

\details{
  The "Variable J" method is a way to estimate the chloroplast CO2 concentration
  \code{Cc} and the mesophyll conductance to CO2 \code{gmc} from combined gas
  exchange and chlorophyll fluorescence measurements, and was originally
  described in Harley et al. (1992). The main idea is that along with \code{Cc},
  the net CO2 assimilation rate (\code{An}), day respiration rate (\code{RL}),
  and CO2 compensation point in the absence of day respiration
  (\code{Gamma_star}) determine the actual RuBP regeneration rate
  (\code{J_actual}) required to support the Calvin-Benson cycle:

  \code{J_actual = (A + RL) * (4 * Cc + 8 * Gamma_star) / (Cc - Gamma_star)}

  This is Equation 6 in Harley et al. (1992). (Note: this equation can be
  derived by solving the equation for \code{Aj} from the FvCB model for
  \code{J}. However, this relationship holds true even when CO2 assimilation is
  not limited by RuBP regeneration. Hence, we distinguish between the actual
  regeneration rate \code{J_actual} and the maximum regeneration rate for a
  given incident light level \code{J}.)

  This equation can be rewritten by using a 1D diffusion equation to replace
  \code{Cc} with \code{Cc = Ci - An / gmc} and then solving for the mesophyll
  conductance. The result is Equation 7 in Harley et al. (1992), which we do not
  reproduce here. The importance of Equation 7 is that it calculates \code{gmc}
  from several quantities that can be measured using gas exchange (\code{Ci},
  \code{An}, and \code{RL}), a quantity whose values can be known beforehand
  (\code{Gamma_star}), and \code{J_actual} (which can be estimated from
  chlorophyll fluorescence measurements). Here we update Equation 7 to include
  \code{alpha_g} and \code{alpha_s} following Busch et al. (2018) (also see
  \code{\link{calculate_c3_assimilation}}.)

  The actual RuBP regeneration rate is related to the incident
  photosynthetically active flux density \code{Qin} and the operating efficiency
  of photosystem II \code{PhiPSII} according to:

  \code{J_actual = alpha_g * beta * Qin * PhiPSII},

  where \code{alpha_g} is the leaf absorptance and \code{beta} is the fraction of
  absorbed light energy directed to photosystem II. \code{Qin} is set by the
  measurement conditions, while \code{PhiPSII} can be estimated from chlorophyll
  fluorescence. However, the values of \code{alpha_g} and \code{beta} are
  generally unknown; \code{beta} in particular is difficult or impossible to
  measure and is often assumed to be 0.5. Thus, while Equation 7 from Harley et
  al. (1992) can be used to estimate \code{gmc}, there is a practical
  uncertainty associated with determining a value of \code{J_actual} to be used
  in Equation 7.

  Moualeu-Ngangue, Chen, & Stutzel (2016) developed a way to address this issue.
  The method from that paper replaces the product of \code{alpha_g} and
  \code{beta} by a single new parameter \code{tau}, and uses it to estimate the
  actual RuBP regeneration from fluoresence (\code{J_F}):

  \code{J_F = tau * Qin * PhiPSII}.

  This new parameter \code{tau} is assumed to be constant across an A-Ci curve,
  and is treated as an unknown whose value will be determined during a fitting
  procedure.

  In this function, the supplied values of \code{Qin}, \code{PhiPSII}, and
  \code{tau} are used to calculate values of \code{J_F}. Then, the values of
  \code{J_F} are used along with Equation 7 from Harley et al. (1992) to
  calculate \code{gmc}. Finally, a 1D diffusion equation is used to calculate
  \code{Cc}.

  \strong{Hard constraints:}

  Most input parameters to the Variable J equations have hard constraints on
  their values which are set by their biochemical or physical interpretation;
  for example, \code{RL} cannot be negative and \code{tau} must lie between 0
  and 1. Yet, because of measurement noise, sometimes it is necessary to use
  values outside these ranges when fitting an A-Ci curve with
  \code{\link{fit_c3_variable_j}}. To accomodate different potential use cases,
  it is possible to selectively apply these hard constraints by specifying
  different values of the \code{hard_constraints} input argument:

  \itemize{
    \item \code{hard_constraints = 0}: Constraints are only placed on inputs
          that are user-supplied and cannot be fit, such as \code{Qin}.

    \item \code{hard_constraints = 1}: Includes the same constraints as when
          \code{hard_constraints} is 0, with the additional constraint that
          all \code{Ci} values must be non-negative.

    \item \code{hard_constraints = 2}: Includes the same constraints as when
          \code{hard_constraints} is 1, which additional constraints on the
          parameters that can be fitted. For example, \code{RL_at_25} must be
          non-negative and \code{tau} must lie between 0 and 1.
  }

  If any input values violate any of the specified constraints, an error message
  will be thrown.

  References:
  \itemize{
    \item Harley, P. C., Loreto, F., Di Marco, G. & Sharkey, T. D. "Theoretical
          Considerations when Estimating the Mesophyll Conductance to CO2 Flux
          by Analysis of the Response of Photosynthesis to CO2" Plant Physiology
          98, 1429–1436 (1992) [\doi{10.1104/pp.98.4.1429}].
    \item Moualeu-Ngangue, D. P., Chen, T.-W. & Stutzel, H. "A new method to
          estimate photosynthetic parameters through net assimilation
          rate-intercellular space CO2 concentration (A-Ci) curve and
          chlorophyll fluorescence measurements" New Phytologist 213, 1543–1554
          (2017) [\doi{10.1111/nph.14260}].
    \item Busch, Sage, & Farquhar, G. D. "Plants increase CO2 uptake by
          assimilating nitrogen via the photorespiratory pathway." Nature Plants
          4, 46–54 (2018) [\doi{10.1038/s41477-017-0065-x}].
  }
}

\value{
  The return value depends on the value of \code{return_exdf}:
  \itemize{
    \item If \code{return_exdf} is \code{TRUE}, the return value is an
          \code{exdf} object with the following columns, calculated as described
          above: \code{J_F}, \code{gmc}, \code{Cc}, \code{tau}, and
          \code{RL_tl}. The category for each of these new columns is
          \code{calculate_c3_variable_j} to indicate that they were created
          using this function.
    \item If \code{return_exdf} is \code{FALSE}, the return value is a list with
          the following named elements: \code{gmc}, \code{Cc}, and \code{J_F}.
          Each element is a numeric vector.
  }
}

\examples{
# Read an example Licor file included in the PhotoGEA package. This file
# includes gas exchange and chlorophyll fluorescence data.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('c3_aci_1.xlsx')
)

# Define a new column that uniquely identifies each curve
licor_file[, 'species_plot'] <-
  paste(licor_file[, 'species'], '-', licor_file[, 'plot'])

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_file <- calculate_temperature_response(licor_file, c3_temperature_param_bernacchi)

# Calculate values of J_F, gmc, and Cc assuming alpha_g = alpha_s = alpha_t = 0,
# Gamma_star_at_25 = 40, tau = 0.2, and RL_at_25 = 0.5.
vj_res <- calculate_c3_variable_j(licor_file, 0, 0, 0, 40, 1.5, 0.55)

# Plot mesophyll conductance against Cc. Note: this information is not very
# meaningful since the values of Gamma_star, tau and RL used above are
# arbitrary.
lattice::xyplot(
  gmc ~ Cc | licor_file[, 'species_plot'],
  data = vj_res$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE,
  xlab = paste0('Chloroplast CO2 concentration (', vj_res$units$Cc, ')'),
  ylab = paste0('Mesophyll conductance to CO2 (', vj_res$units$gmc, ')')
)
}

\concept{exdf}
