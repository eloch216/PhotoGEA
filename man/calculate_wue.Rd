\name{calculate_wue}

\alias{calculate_wue}

\title{Calculate intrinsic water use efficiency}

\description{
  Calculates the intrinsic water use efficiency (\code{iWUE}). This function can
  accomodate alternative column names for the variables taken from the data
  file in case they change at some point in the future. This function also
  checks the units of each required column and will produce an error if any
  units are incorrect.
}

\usage{
  calculate_wue(
    exdf_obj,
    calculate_c3 = FALSE,
    a_column_name = 'A',
    ca_column_name = 'Ca',
    cc_column_name = 'Cc',
    ci_column_name = 'Ci',
    e_column_name = 'E',
    gmc_column_name = 'gmc_tl',
    gsw_column_name = 'gsw',
    h2o_a_column_name = 'H2O_s',
    h2o_i_column_name = 'H2O_i',
    total_pressure_column_name = 'total_pressure'
  )
}

\arguments{
  \item{exdf_obj}{An \code{exdf} object.}

  \item{calculate_c3}{
    A logical variable indicating whether to calculate additional variables that
    can be useful for C3 plants (\code{g_ratio} and \code{drawdown_ct}). Note
    that these quantities require values of mesophyll conductance and \code{Cc},
    so it is not always possible to calculate them.
  }

  \item{a_column_name}{
    The name of the column in \code{exdf_obj} that contains the net CO2
    assimilation rate in \code{micromol m^(-2) s^(-1)}.
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

  \item{e_column_name}{
    The name of the column in \code{licor_exdf} that contains the transpiration
    rate in \code{mol m^(-2) s^(-1)}.
  }

  \item{gmc_column_name}{
    The name of the column in \code{licor_exdf} that contains the mesophyll
    conductance to CO2 at leaf temperature in \code{mol m^(-2) s^(-1) bar^(-1)}.
  }

  \item{gsw_column_name}{
    The name of the column in \code{licor_exdf} that contains the stomatal
    conductance to water vapor in \code{mol m^(-2) s^(-1)}.
  }

  \item{h2o_a_column_name}{
    The name of the column in \code{exdf_obj} that contains the water vapor
    concentration in the air surrounding the leaf (i.e., the ambient water
    vapor concentration) in \code{mmol mol^(-1)}.
  }

  \item{h2o_i_column_name}{
    The name of the column in \code{exdf_obj} that contains the water vapor
    concentration in the leaf's intercellular air spaces in
    \code{mmol mol^(-1)}. Typically this value is calculated using
    \code{\link{calculate_gas_properties}}.
  }

  \item{total_pressure_column_name}{
    The name of the column in \code{exdf_obj} that contains the total pressure
    in \code{bar}. Typically this value is calculated using
    \code{\link{calculate_total_pressure}}.
  }
}

\details{
  Leaf-level water use efficiency (\code{lWUE}) is defined as the ratio of net
  CO2 assimilation (\code{An}) to transpiration (\code{E}):

  \code{lWUE = An / E}.

  This quantity can also be expressed in terms of water and CO2 concentrations:

  \code{lWUE = 0.6 * Ca * (1 - Ci / Ca) / (H2Oi - H2Oa)}.

  Here, \code{Ca} and \code{Ci} are the atmospheric and intercellular CO2
  concentrations, and \code{H2Oa} and \code{H2Oi} are the atmospheric and
  intercellular water vapor concentrations. If differences in \code{lWUE} are
  measured between different groups of plants, it can be helpful to separately
  investigate \code{Ci / Ca} and \code{H2Oi - H2Oa} to see which factor is
  driving the differences.

  The intrinsic water use efficiency \code{iWUE} is a measure of leaf-level
  water use efficiency, and it is defined to be the ratio \code{An} and the
  stomatal conductance to H2O diffusion (\code{gsw}):

  \code{iWUE = An / gsw}.

  For C3 plants, \code{iWUE} can be reexpressed as

  \code{iWUE = (gmc / gsw) / (1 + (gmc / gsw)) * (Ca - Cc)},

  where \code{gmc} is the mesophyll conductance to CO2 diffusion and \code{Cc}
  is the chloroplast CO2 concentration. If differences in \code{iWUE} are
  measured between different groups of plants, it can be helpful to separately
  investigate \code{gmc / gsw} and \code{Ca - Cc} to see which factor is driving
  the differences.

  Note: both measures of water use efficiency depend directly or indirectly on
  stomatal conductance. Stomata are notoriously slow to reach steady-state, but
  water use efficiency is only reliable at steady-state. For this reason, it is
  recommended to only analyze water use efficiency for gas exchange measurements
  where stomatal conductance has stabilized. For an A-Ci or A-Q curve, only the
  first measured point has typically reached steady-state stomatal conductance.
  On the other hand, for a Ball-Berry curve, all measured points should have
  reached steady-state stomatal conductance.

  For more details about these quantities, see Leakey et al. "Water Use
  Efficiency as a Constraint and Target for Improving the Resilience and
  Productivity of C3 and C4 Crops." Annual Review of Plant Biology 70 (1):
  781–808 (2019) [\doi{10.1146/annurev-arplant-042817-040305}].

  In this function, the following variables are calculated:
  \itemize{
    \item \code{lWUE}, given by \code{iWUE = An / E}
    \item \code{Cia_ratio}, given by \code{Cia_ratio = Ci / Ca}
    \item \code{drawdown_sw}, given by \code{drawdown_sw = H2Oi - H2Oa} (this
          is the drawdown of water vapor across the stomata)
    \item \code{iWUE}, given by \code{iWUE = An / gsw}
    \item \code{g_ratio}, given by \code{g_ratio = gmc / gsw}
    \item \code{drawdown_ct}, given by \code{drawdown_ct = Ca - Cc} (this is the
          total drawdown of CO2 from the ambient air to the chloroplast)
  }

  Note: \code{g_ratio} and \code{drawdown_ct} are only calculated if
  \code{calculate_c3} is \code{TRUE}.
}

\value{
  An \code{exdf} object based on \code{exdf_obj} that includes the quantities
  listed above, along with their units. The category of each of these new
  columns is \code{calculate_wue} to indicate that it was created using this
  function.
}

\examples{
# Read an example Licor file included in the PhotoGEA package and calculate the
# water use efficiency.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- calculate_total_pressure(licor_file)

licor_file <- calculate_temperature_response(licor_file, c3_temperature_param_sharkey)

licor_file <- calculate_gas_properties(licor_file)

licor_file <- apply_gm(licor_file, gmc_at_25 = 0.5)

licor_file <- calculate_wue(licor_file, calculate_c3 = TRUE)

licor_file$units$iWUE      # View the units of the new `iWUE` column
licor_file$categories$iWUE # View the category of the new `iWUE` column
licor_file[, 'iWUE']       # View the values of the new `iWUE` column
}

\concept{exdf}
