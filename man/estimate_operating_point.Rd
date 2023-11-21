\name{estimate_operating_point}

\alias{estimate_operating_point}

\title{Estimate the operating point from an A-Ci curve}

\description{
  Uses linear interpolation to estimate \code{Cc}, \code{Ci}, and \code{An} at
  atmospheric CO2 concentration from the data in the \code{exdf} object, which
  should represent a single A-Ci curve. This function can accomodate alternative
  column names for the variables taken from the data file in case they change
  at some point in the future. This function also checks the units of each
  required column and will produce an error if any units are incorrect.
}

\usage{
  estimate_operating_point(
    aci_exdf,
    ca_atmospheric = 420,
    type = 'c3',
    a_column_name = 'A',
    ca_column_name = 'Ca',
    cc_column_name = 'Cc',
    ci_column_name = 'Ci',
    pcm_column_name = 'PCm',
    return_list = FALSE
  )
}

\arguments{
  \item{aci_exdf}{
    An \code{exdf} object representing one CO2 response curve.
  }

  \item{ca_atmospheric}{
    The atmospheric CO2 concentration (with units of \code{micromol mol^(-1)});
    this will be used to estimate the operating point. The default value is the
    approximate global average during the early 2020s.
  }

  \item{type}{
    The type of photosynthesis: either \code{'c3'} or \code{'c4'}.
  }

  \item{a_column_name}{
    The name of the column in \code{aci_exdf} that contains the net
    assimilation in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ca_column_name}{
    The name of the column in \code{aci_exdf} that contains the
    ambient CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{cc_column_name}{
    The name of the column in \code{aci_exdf} that contains the
    chloroplastic CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{ci_column_name}{
    The name of the column in \code{aci_exdf} that contains the
    intercellular CO2 concentration in \code{micromol mol^(-1)}.
  }

  \item{pcm_column_name}{
    The name of the column in \code{aci_exdf} that contains the partial pressure
    of CO2 in the mesophyll, expressed in \code{microbar}.
  }

  \item{return_list}{
    A logical value indicating whether or not to return the results as a list.
    Most users will only need to use \code{return_list = TRUE};
    \code{return_list = FALSE} is used internally by other functions in the
    \code{PhotoGEA} package.
  }

}

\details{
  When analyzing or interpreting A-Ci curves, it is often useful to determine
  the values of \code{Ci} and \code{An} that correspond to typical growth
  conditions (where \code{Ca} is set to the atmospheric value). Together, these
  special values of \code{Ci} and \code{An} specify the "operating point" of the
  leaf.

  However, for a variety of practical reasons, most A-Ci curves do not actually
  contain a measurement point where \code{Ca} is at the atmospheric value.
  Nevertheless, it is possible to apply linear interpolation to the observed
  \code{Ci - Ca} and \code{An - Ca} relations to estimate the operating point.
  This function automates that procedure. It also calculates the operating
  values of \code{Cc} (for \code{c3} A-Ci curves) and \code{PCm} (for \code{c4}
  A-Ci curves).

  This function assumes that \code{aci_exdf} represents a single
  A-Ci curve. Typically, this function is not directly called by users because
  the fitting functions \code{\link{fit_c3_aci}} and \code{\link{fit_c4_aci}}
  automatically use this function to determine the operating point.
}

\value{
  The return value depends on \code{return_list} and \code{type}.

  When \code{return_list} is \code{FALSE}, this function returns an \code{exdf}
  object based on \code{aci_exdf} that includes its identifier columns as well
  as values of \code{operating_Ci}, \code{operating_An}, and \code{operating_Cc}
  (or \code{operating_PCm}) in columns with those names.

  When \code{return_list} is \code{TRUE}, this function returns a list with the
  following named elements: \code{operating_Ci}, \code{operating_An},
  \code{operating_Cc} (or \code{operating_PCm}), and \code{operating_exdf}. The
  first three are numeric values as described above, while \code{operating_exdf}
  is an \code{exdf} object with one row that can be passed to
  \code{\link{calculate_c3_assimilation}} or
  \code{\link{calculate_c4_assimilation}} in order to estimate the operating
  \code{An} from a photosynthesis model.
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  system.file('extdata', 'c3_aci_1.xlsx', package = 'PhotoGEA', mustWork = TRUE)
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

# Determine the operating point for just one curve from the data set
one_result <- estimate_operating_point(
  licor_file[licor_file[, 'species_plot'] == 'tobacco - 1', , TRUE]
)

one_result[, 'operating_Cc']
one_result[, 'operating_Ci']
one_result[, 'operating_An']
}

\concept{exdf}
