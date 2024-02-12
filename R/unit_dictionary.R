# Specify standard units for some types of quantities
dimensionless     <- 'dimensionless'
micromol_flux     <- 'micromol m^(-2) s^(-1)'
micromol_fraction <- 'micromol mol^(-1)'
millimol_fraction <- 'mmol mol^(-1)'
pressure          <- 'bar'

# Specify units for some important parameters
unit_dictionary <- list(
    alpha_g        = dimensionless,
    Ca             = micromol_fraction,
    Cc             = micromol_fraction,
    Ci             = micromol_fraction,
    Gamma_star     = micromol_fraction,
    J_at_25        = micromol_flux,
    J_norm         = 'normalized to J at 25 degrees C',
    Kc             = micromol_fraction,
    Ko             = millimol_fraction,
    Rd_at_25       = micromol_flux,
    Rd_norm        = 'normalized to Rd at 25 degrees C',
    tau            = dimensionless,
    total_pressure = pressure,
    TPU            = micromol_flux,
    Vcmax_at_25    = micromol_flux,
    Vcmax_norm     = 'normalized to Vcmax at 25 degrees C',
    Vpmax_at_25    = micromol_flux,
    Vpr            = micromol_flux
)
