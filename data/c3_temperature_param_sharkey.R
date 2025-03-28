c3_temperature_param_sharkey <- list(
    Kc              = list(type = 'Arrhenius', c = 35.9774 + PhotoGEA:::c_pa_to_ppm, Ea = 80.99, units = 'micromol mol^(-1)'),
    Ko              = list(type = 'Arrhenius', c = 12.3772 + PhotoGEA:::c_pa_to_ppm, Ea = 23.72, units = 'mmol mol^(-1)'),
    Gamma_star_norm = list(type = 'Arrhenius', c = 24.46 / PhotoGEA:::f,             Ea = 24.46, units = 'normalized to Gamma_star at 25 degrees C'),
    Vcmax_norm      = list(type = 'Arrhenius', c = 26.355,                           Ea = 65.33, units = 'normalized to Vcmax at 25 degrees C'),
    J_norm          = list(type = 'Arrhenius', c = 17.71,                            Ea = 43.9,  units = 'normalized to J at 25 degrees C'),
    RL_norm         = list(type = 'Arrhenius', c = 18.7145,                          Ea = 46.39, units = 'normalized to RL at 25 degrees C'),
    gmc_norm        = list(type = 'Johnson',   c = 20.01, Ha = 49.6, Hd = 437.4, S = 1.4,        units = 'normalized to gmc at 25 degrees C'),
    Tp_norm         = list(type = 'Johnson',   c = 21.46, Ha = 53.1, Hd = 201.8, S = 0.65,       units = 'normalized to Tp at 25 degrees C')
)
