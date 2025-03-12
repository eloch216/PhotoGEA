# Specify standard units for some types of quantities
conductance       <- 'mol m^(-2) s^(-1)'
conductance_bar   <- 'mol m^(-2) s^(-1) bar^(-1)'
dimensionless     <- 'dimensionless'
micromol_flux     <- 'micromol m^(-2) s^(-1)'
micromol_fraction <- 'micromol mol^(-1)'
millimol_fraction <- 'mmol mol^(-1)'
pressure          <- 'bar'
temperature       <- 'degrees C'

# Specify units for some important parameters
unit_dictionary_list <- list(
    A              = micromol_flux,
    Ainitial       = micromol_flux,
    alpha_g        = dimensionless,
    alpha_j_at_25  = dimensionless,
    alpha_j_norm   = 'normalized to alpha_j at 25 degrees C',
    alpha_old      = dimensionless,
    alpha_psii     = dimensionless,
    alpha_s        = dimensionless,
    alpha_t        = dimensionless,
    bb_index       = conductance,
    c4_curvature   = dimensionless,
    c4_slope       = conductance, # this is not a conductance, but happens to have the same units as a conductance
    Ca             = micromol_fraction,
    Cc             = micromol_fraction,
    Ci             = micromol_fraction,
    CO2_r          = micromol_fraction,
    CO2_s          = micromol_fraction,
    CorrFact       = NA,
    Csurface       = micromol_fraction,
    Flow           = 'micromol s^(-1)',
    Gamma_star     = micromol_fraction,
    gbs            = conductance_bar,
    gmc_at_25      = conductance_bar,
    gmc_norm       = 'normalized to gmc at 25 degrees C',
    gsw            = conductance,
    H2O_r          = millimol_fraction,
    H2O_s          = millimol_fraction,
    I2             = micromol_flux,
    I2_at_25       = micromol_flux,
    I2_tl          = micromol_flux,
    J              = micromol_flux,
    J_at_25        = micromol_flux,
    J_at_25_lower  = micromol_flux,
    J_at_25_upper  = micromol_flux,
    J_norm         = 'normalized to J at 25 degrees C',
    J_tl_avg       = micromol_flux,
    J_tl_avg_lower = micromol_flux,
    J_tl_avg_upper = micromol_flux,
    Jmax_at_25     = micromol_flux,
    Jmax_norm      = 'normalized to Jmax at 25 degrees C',
    Jmax_tl        = micromol_flux,
    Kc             = micromol_fraction,
    Ko             = millimol_fraction,
    oxygen         = 'percent',
    Qin            = micromol_flux,
    rL             = micromol_flux,
    RL_at_25       = micromol_flux,
    RL_norm        = 'normalized to RL at 25 degrees C',
    Rm_frac        = dimensionless,
    S              = 'cm^2',
    tau            = dimensionless,
    theta_j_at_25  = dimensionless,
    theta_j_norm   = 'normalized to theta_j at 25 degrees C',
    Tleaf_avg      = temperature,
    TleafCnd       = temperature,
    total_pressure = pressure,
    Tp_at_25       = micromol_flux,
    Tp_norm        = 'normalized to Tp at 25 degrees C',
    Vcmax_at_25    = micromol_flux,
    Vcmax_norm     = 'normalized to Vcmax at 25 degrees C',
    Vmax           = micromol_flux,
    VPDleaf        = 'kPa',
    Vpmax_at_25    = micromol_flux,
    Vpr            = micromol_flux
)

unit_dictionary <- function(quantity_name) {
    if (!quantity_name %in% names(unit_dictionary_list)) {
        msg <- paste0(
            'Units were requested for a quantity named `', quantity_name,
            '`, but it is not included in the unit dictionary'
        )

        stop(msg)
    }

    unit_dictionary_list[[quantity_name]]
}
