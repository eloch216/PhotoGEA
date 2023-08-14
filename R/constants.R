# Here we define some constants so they can be used consistently across
# calculations in separate functions.

ISOTOPE_CONSTANTS <- list(
    a_b             = 2.9,  # ppt; isotopic fractionation during diffusion through the laminar boundary layer
    a_s             = 4.4,  # ppt; isotopic fractionation during diffusion through air
    a_m             = 1.8,  # ppt; isotopic fractionation during diffusion through water
    b               = 29,   # ppt; isotopic fractionation during Rubisco carboxylation
    b_prime_4       = -5.7, # ppt; combined isotopic fractionation during CO2 dissolution, hydration, and PEPc activity at 25 degrees C
    delta_Ca_growth = -8,   # ppt; isotope ratio of ambient air during plant growth
    s               = 1.8   # ppt; isotopic fractionation during leakage from bundle-sheath cells
)
