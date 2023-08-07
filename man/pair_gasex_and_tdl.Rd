\name{pair_gasex_and_tdl}

\alias{pair_gasex_and_tdl}

\title{Pair gas exchange and TDL data}

\description{
  Identifies the closest TDL cycle corresponding to each entry in the gas
  exchange data and adds the TDL data to the gas exchange data.
}

\usage{
  pair_gasex_and_tdl(
    gasex_exdf,
    tdl_exdf,
    max_allowed_time_difference = 1,
    gasex_timestamp_column_name = 'time',
    tdl_timestamp_column_name = 'TIMESTAMP'
  )
}

\arguments{
  \item{gasex_exdf}{
    An \code{exdf} object representing data from a photosynthetic gas exchange
    measurement system.
  }

  \item{tdl_exdf}{
    An \code{exdf} object representing calibrated data from a tunable diode
    laser absorption spectroscopy system. Typically this is the output from
    applying \code{\link{process_tdl_cycle_erml}} or
    \code{\link{process_tdl_cycle_polynomial}} to a set of uncalibrated TDL data.
  }

  \item{max_allowed_time_difference}{
    The maximum time difference (in minutes) to allow between gas exchange and
    TDL timestamp values.
  }

  \item{gasex_timestamp_column_name}{
    The name of the column in \code{gasex_exdf} that contains the timestamp
    values.
  }

  \item{tdl_timestamp_column_name}{
    The name of the column in \code{tdl_exdf} that contains the timestamp
    values.
  }
}

\details{
  When making combined gas exchange and isotope discrimination measurements
  using a portable photosynthetic gas exchange system (such as a Licor LI-6800)
  coupled with a tunable diode laser (TDL) absorption spectroscopy system, the
  TDL's gas handling system cycles through several gas lines (or sites) by
  opening and closing valves. When analyzing such data, a key step is to
  combine TDL and gas exchange data that were measured at the same times.

  The \code{pair_gasex_and_tdl} function performs this operation by locating the
  TDL cycle whose timestamp is closest to each Licor file entry. Then, the
  12C, 13C, total CO2, and delta_13C values measured by the TDL from the Licor's
  sample and reference lines during that cycle are added to the gas exchange
  data as new columns.
}

\value{
  An \code{exdf} object based on \code{gasex_exdf} that includes TDL values
  measured at the same times as the original gas exchange logs. Several new
  columns are added: \code{'cycle_num'}, \code{'tdl_time_s'},
  \code{'calibrated_12c_s'}, \code{'calibrated_13c_s'}, \code{'total_CO2_s'},
  \code{'delta_C13_s'}, \code{'tdl_time_r'}, \code{'calibrated_12c_r'},
  \code{'calibrated_13c_r'}, \code{'total_CO2_r'}, and \code{'delta_C13_r'}.
  Variables with \code{'_s'} in the name refer to TDL measurements from the
  Licor sample line, and \code{'_r'} indicates the reference line. The category
  of each new column is \code{pair_gasex_and_tdl} to indicate that it was
  created using this function.
}

\concept{exdf}
