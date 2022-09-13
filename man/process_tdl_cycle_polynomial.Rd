\name{process_tdl_cycle_polynomial}

\alias{process_tdl_cycle_polynomial}

\title{Process cycles from the ERML TDL}

\description{
  Uses the 12C and 13C signal from the calibration lines of a tunable diode
  laser (TDL) to determine correction factors and apply them to the sample
  lines. Applicable for a system with two or more reference tanks whose 12C and
  13C concentrations are known beforehand.
}

\usage{
  process_tdl_cycle_polynomial(
    tdl_cycle,
    poly_order,
    reference_tanks,
    f_other = 0.00474,
    R_VPDB = 0.0111797,
    valve_column_name = 'valve_number',
    raw_12c_colname = 'Conc12C_Avg',
    raw_13c_colname = 'Conc13C_Avg'
  )
}

\arguments{
  \item{tdl_cycle}{
    An \code{exdf} object representing one cycle of TDL data.
  }

  \item{poly_order}{
    The order of the polynomial to fit, where 1 indicates a linear fit, 2
    indicates a quadratic fit, etc. This argument will be passed to
    \code{stats::poly} during the fitting procedure.
  }

  \item{reference_tanks}{
    A list where each element is a list with three named elements: \code{valve},
    \code{conc_12C}, and \code{conc_13C}. \code{valve} should indicate the valve
    number for the reference tank, and the other two elements should indicate
    the known concentrations of 12C and 13C in the tank.
  }

  \item{f_other}{
    The fraction of CO2 that is not 13C16O16O or 12C16O16O.
  }

  \item{R_VPDB}{
    The Vienna Pee Dee Belemnite reference value.
  }

  \item{valve_column_name}{
    The name of the column in \code{tdl_cycle} that contains the valve number.
  }

  \item{raw_12c_colname}{
    The name of the column in \code{tdl_cycle} that contains the 12C signal.
  }

  \item{raw_13c_colname}{
    The name of the column in \code{tdl_cycle} that contains the 13C signal.
  }
}

\details{
  This function applies a simple correction to the measured values of 12C and
  13C. This correction is based on the fact that each reference tank has both a
  true concentration (which is known beforehand) and a measured concentration
  (from the TDL) of each isotope. Using this information, it is possible to
  perform a polynomial fit of true vs. measured concentrations; in other words,
  it is possible to identify a polynomial function that determines true
  concentrations from measured ones. This function can then be applied to tanks
  whose concentration is not known beforehand; in this case, it provides an
  estimate of the true concentration, otherwise referred to as a calibrated
  value.

  This function assumes that \code{tdl_cycle} represents a single TDL
  measurement cycle. To process multiple cycles at once, this function is often
  used along with \code{\link{by.exdf}} and \code{\link{consolidate}}.
}

\value{
  A list with two elements:
  \itemize{
    \item \code{tdl_data}: An \code{exdf} object containing the original content
          of \code{tdl_cycle} and several new columns:
          \code{'calibrated_12c'}, \code{'calibrated_13c'},
          \code{'total_mixing_ratio'}, and \code{'total_isotope_ratio'}.
    \item \code{calibration_parameters}: An \code{exdf} object describing the
          fitted polynomial coefficients.
  }
}

\examples{
# TODO: add an example

}

\concept{exdf}
