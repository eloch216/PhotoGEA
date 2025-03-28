c3_temperature_param_flat <- list(
    Kc              = list(type = 'Arrhenius', c = log(269.339), Ea = 0, units = 'micromol mol^(-1)'),
    Ko              = list(type = 'Arrhenius', c = log(163.715), Ea = 0, units = 'mmol mol^(-1)'),
    Gamma_star_norm = list(type = 'Arrhenius', c = 0,            Ea = 0, units = 'normalized to Gamma_star at 25 degrees C'),
    Vcmax_norm      = list(type = 'Arrhenius', c = 0,            Ea = 0, units = 'normalized to Vcmax at 25 degrees C'),
    J_norm          = list(type = 'Arrhenius', c = 0,            Ea = 0, units = 'normalized to J at 25 degrees C'),
    RL_norm         = list(type = 'Arrhenius', c = 0,            Ea = 0, units = 'normalized to RL at 25 degrees C'),
    gmc_norm        = list(type = 'Arrhenius', c = 0,            Ea = 0, units = 'normalized to gmc at 25 degrees C'),
    Tp_norm         = list(type = 'Arrhenius', c = 0,            Ea = 0, units = 'normalized to Tp at 25 degrees C')
)
