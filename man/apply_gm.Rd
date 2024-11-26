\name{apply_gm}

\alias{apply_gm}

\title{Calculate CO2 concentration in the chloroplast or mesophyll}

\description{
  Calculates CO2 concentration in the chloroplast or mesophyll, the CO2 drawdown
  across the stomata, and the CO2 drawdown across the mesophyll. This function
  can accomodate alternative column names for the variables taken from the Licor
  file in case they change at some point in the future. This function also
  checks the units of each required column and will produce an error if any
  units are incorrect.
}

\usage{
  apply_gm(
    exdf_obj,
    gmc_at_25 = '',
    photosynthesis_type = 'C3',
    calculate_drawdown = TRUE,
    a_column_name = 'A',
    ca_column_name = 'Ca',
    ci_column_name = 'Ci',
    gmc_norm_column_name = 'gmc_norm',
    total_pressure_column_name = 'total_pressure',
    perform_checks = TRUE,
    return_exdf = TRUE
  )
}

\arguments{
  \item{exdf_obj}{
    An \code{exdf} object, typically representing data from a Licor gas exchange
    measurement system.
  }

  \item{gmc_at_25}{
    The mesophyll conductance to CO2 diffusion at 25 degrees C, expressed in
    \code{mol m^(-2) s^(-1) bar^(-1)}. In the absence of other reliable
    information, \code{gmc_at_25} is often assumed to be infinitely large. If
    \code{gmc_at_25} is not a number, then there must be a column in
    \code{exdf_obj} called \code{gmc_at_25} with appropriate units. A numeric
    value supplied here will overwrite the values in the \code{gmc_at_25}
    column of \code{exdf_obj} if it exists.
  }

  \item{photosynthesis_type}{
    A string indicating the type of photosynthesis being considered (either
    \code{'C3'} or \code{'C4'}).
  }

  \item{calculate_drawdown}{
    A logical value indicating whether to calculate drawdown values.
  }

  \item{a_column_name}{
    The name of the column in \code{exdf_obj} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ca_column_name}{
    The name of the column in \code{exdf_obj} that contains the ambient CO2
    concentration in the chamber in \code{micromol mol^(-1)}.
  }

  \item{ci_column_name}{
    The name of the column in \code{exdf_obj} that contains the intercellular
    CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{gmc_norm_column_name}{
    The name of the column in \code{exdf_obj} that contains the normalized
    mesophyll conductance values (with units of
    \code{normalized to gmc at 25 degrees C}).
  }

  \item{total_pressure_column_name}{
    The name of the column in \code{exdf_obj} that contains the total pressure
    in \code{bar}.
  }

  \item{perform_checks}{
    A logical value indicating whether to check units for the required columns.
    This should almost always be \code{TRUE}. The option to disable these checks
    is only intended to be used when \code{\link{fit_c3_aci}} calls this
    function, since performing these checks many times repeatedly slows down the
    fitting procedure.
  }

  \item{return_exdf}{
    A logical value indicating whether to return an \code{exdf} object. This
    should almost always be \code{TRUE}. The option to return a vector is mainly
    intended to be used when  \code{\link{fit_c3_aci}} calls this function,
    since creating an \code{exdf} object to return will slow down the fitting
    procedure.
  }
}

\details{
  For a C3 plant, the mesophyll conductance to CO2 (\code{gmc}) is said to be
  the conductance satisfying the following one-dimensional flux-conductance
  equation:

  (1) \code{An = gmc * (PCi - PCc)}

  where \code{An} is the net CO2 assimilation rate, \code{PCi} is the partial
  pressure of CO2 in the intercellular spaces, and \code{PCc} is the partial
  pressure of CO2 in the chloroplast.  A key underlying assumption for this
  equation is that the flow of CO2 has reached a steady state; in this case, the
  flow across the stomata is equal to the flow across the mesophyll.

  This equation can be rearranged to calculate \code{PCc}:

  (2) \code{PCc = PCi - An / gmc}

  This version of the equation can be found in many places, for example, as
  Equation 4 in Sharkey et al. "Fitting photosynthetic carbon dioxide response
  curves for C3 leaves" Plant, Cell & Environment 30, 1035–1040 (2007)
  [\doi{10.1111/j.1365-3040.2007.01710.x}].

  It is common to express the partial pressures in \code{microbar} and the
  assimilation rate in \code{micromol m^(-2) s^(-1)}; in this case, the units of
  mesophyll conductance become \code{mol m^(-2) s^(-1) bar^(-1)}.

  Licor measurement systems provide CO2 levels as relative concentrations with
  units of parts per million (\code{ppm}), or equivalently,
  \code{micromol mol^(-1)}. Concentrations and partial pressures are related by
  the total gas pressure according to:

  (3) \code{partial_pressure = total_pressure * relative_concentration}

  Thus, it is also possible to calculate the CO2 concentration in the
  choloroplast (\code{Cc}) using the following equation:

  (4) \code{Cc = Ci - An / (gmc * P)}

  where \code{Ci} is the intercellular CO2 concentration and \code{P} is the
  total pressure. In this function, Equation (4) is used to calculate \code{Cc},
  where the total pressure is given by the sum of the atmospheric pressure and
  the chamber overpressure.

  When a plant is photosynthesizing, it draws CO2 into its chloroplasts, and
  this flow is driven by a concentration gradient. In other words, as CO2 flows
  from the ambient air across the stomata to the intercellular spaces and then
  across the mesophyll into the chloroplast, there is a decrease in CO2
  concentration at each step. Sometimes it is useful to calculate these changes,
  which are usually referred to as "CO2 drawdown" values. So, in addition to
  \code{Ci}, this function (optionally) calculates the drawdown of CO2 across
  the stomata (\code{drawndown_cs = Ca - Ci}) and the drawdown of CO2 across the
  mesophyll (\code{drawdown_cm = Ci - Cc}).

  _Note_: mesophyll conductance is not specified in typical Licor files, so it
  usually must be added using \code{\link{set_variable}} before calling
  \code{apply_gm}.

  For a C4 plant, mesophyll conductance instead refers to the conductance
  associated with the flow of CO2 from the intercellular spaces into the
  mesophyll (rather than into the chloroplast). In this case, the equations
  above just require a small modification where \code{Pcc} and \code{Cc} are
  replaced by \code{PCm} and \code{Cm}, the partial pressure and concentration
  of CO2 in the mesophyll.
}

\value{
  The return value depends on the value of \code{return_exdf}:
  \itemize{
    \item If \code{return_exdf} is \code{TRUE}, the return value is an
          \code{exdf} object based on \code{exdf_obj} with the following
          columns, calculated as described above: \code{Pci} and \code{Ci}
          (for C3 plants) or \code{PCm} and \code{Cm} (for C4 plants),
          \code{drawndown_s}, and \code{drawdown_cm}. The category for each of
          these new columns is \code{apply_gm} to indicate that they were
          created using this function.
    \item If \code{return_exdf} is \code{FALSE}, the return value is a list with
          a single named element (\code{internal_c}), which contains values of
          \code{Cc} or \code{Cm} as a numeric vector.
  }
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Calculate temperature-dependent parameter values, including gmc_norm
licor_file <- calculate_temperature_response(licor_file, c3_temperature_param_sharkey)

# Calculate Cc and drawdowns assuming a mesophyll conductance of
# 1 mol / m^2 / s / bar at 25 degrees C
licor_file <- apply_gm(licor_file, 1)

licor_file$units$Cc      # View the units of the new `Cc` column
licor_file$categories$Cc # View the category of the new `Cc` column
licor_file[, 'Cc']       # View the values of the new `Cc` column
}

\concept{exdf}
