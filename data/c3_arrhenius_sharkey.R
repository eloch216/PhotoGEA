c3_arrhenius_sharkey <- list(
    Kc =         list(c = 35.9774 + PhotoGEA:::c_pa_to_ppm, Ea = 80.99, units = 'micromol mol^(-1)'),
    Ko =         list(c = 12.3772 + PhotoGEA:::c_pa_to_ppm, Ea = 23.72, units = 'mmol mol^(-1)'),
    Gamma_star = list(c = 11.187  + PhotoGEA:::c_pa_to_ppm, Ea = 24.46, units = 'micromol mol^(-1)'),
    Vcmax_norm = list(c = 26.355,                           Ea = 65.33, units = 'normalized to Vcmax at 25 degrees C'),
    J_norm =     list(c = 17.71,                            Ea = 43.9,  units = 'normalized to J at 25 degrees C'),
    RL_norm =    list(c = 18.7145,                          Ea = 46.39, units = 'normalized to RL at 25 degrees C')
)
