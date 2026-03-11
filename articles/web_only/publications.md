# PhotoGEA Publications

## Overview

This page contains a list of all publications about or using the
PhotoGEA R package, including links to publicly-available analysis
scripts whenever possible. If we are missing any publications, please
let us know!

## Publications About PhotoGEA

The PhotoGEA R package was first described in:

1.  Lochocki, E. B., Salesse-Smith, C. E. & McGrath, J. M. “PhotoGEA: An
    R Package for Closer Fitting of Photosynthetic Gas Exchange Data
    With Non-Gaussian Confidence Interval Estimation.” *Plant, Cell &
    Environment* **48**, 5104–5119 (2025).

    > - [DOI:10.1111/pce.15501](https://doi.org/10.1111/pce.15501)
    >
    > - PhotoGEA version: `v1.1.0`
    >
    > - Analysis scripts:
    >   <https://github.com/ripeproject/PhotoGEA-paper/>

If you use PhotoGEA for your own work, please cite this publication and
specify the version of the package you used.

## Publications Using PhotoGEA

2.  Salesse-Smith, C. E. et *al*. “Greater mesophyll conductance and
    leaf photosynthesis in the field through modified cell wall porosity
    and thickness via AtCGR3 expression in tobacco.” *Plant
    Biotechnology Journal* **22**, 2504–2517 (2024).

    > - [DOI:10.1111/pbi.14364](https://doi.org/10.1111/pbi.14364)
    >
    > - PhotoGEA used for: Calculating mesophyll conductance from
    >   tunable diode laser absorption measurements
    >   (`calculate_gm_busch`), C₃*A-C_(i)* fitting (`fit_c3_aci`),
    >   C₃*A-C_(i)* + CF Variable *J* fitting (`fit_c3_variable_j`)
    >
    > - PhotoGEA version: `v0.10.0`
    >
    > - Analysis scripts:
    >   <https://github.com/ripeproject/CGR3-tobacco-2024>

3.  Pelech, E. A., Stutz, S. S., Wang, Y., Lochocki, E. B. & Long, S. P.
    “Have We Selected for Higher Mesophyll Conductance in Domesticating
    Soybean?” *Plant, Cell & Environment* **48**, 1594–1607 (2025).

    > - [DOI:10.1111/pce.15206](https://doi.org/10.1111/pce.15206)
    >
    > - PhotoGEA used for: Calculating limitations to C₃ photosynthesis
    >   (`calculate_c3_limitations_grassi` and
    >   `calculate_c3_limitations_warren`)
    >
    > - PhotoGEA version: `v1.0.0`
    >
    > - Analysis scripts: <https://doi.org/10.13012/B2IDB-7809185_V2>

4.  Salesse-Smith, C. E. et al. “Adapting C₄ photosynthesis to
    atmospheric change and increasing productivity by elevating Rubisco
    content in sorghum and sugarcane.” *Proceedings of the National
    Academy of Sciences* **122**, e2419943122 (2025).

    > - [DOI:10.1073/pnas.2419943122](https://doi.org/10.1073/pnas.2419943122)
    >
    > - PhotoGEA used for: Fitting C₄*A-C_(i)* curves using the
    >   mechanistic model (`fit_c4_aci`), reading and plotting induction
    >   curves (`read_gasex_file`)
    >
    > - PhotoGEA version: `v0.11.0`
    >
    > - Analysis scripts:
    >   <https://github.com/cabbi-bio/sorghum-sugarcane-RBCS-RAF1-2024>

5.  Lochocki, E. B. & McGrath, J. M. “Widely Used Variants of the
    Farquhar-von-Caemmerer-Berry Model Can Cause Errors in Parameter
    Estimation.” *in silico Plants* diaf014 (2025)

    > - [DOI:10.1093/insilicoplants/diaf014](https://doi.org/10.1093/insilicoplants/diaf014)
    >
    > - PhotoGEA used for: Fitting C₃*A-C_(i)* curves (`fit_c3_aci`)
    >   using several variants of the FvCB model
    >
    > - PhotoGEA version: `v1.2.0`
    >
    > - Analysis scripts: <https://github.com/ripeproject/FvCB-min-A>
    >
    > - Preprint available from *bioRxiv*:
    >   [DOI:10.1101/2025.03.11.642611](https://doi.org/10.1101/2025.03.11.642611)

6.  Tamang, B. G., Bernard, G., Bernacchi, C. J., Diers, B. W. &
    Ainsworth, E. A. “Bigger is not always better: Optimizing leaf area
    index with narrow leaf shape in soybean.” *bioRxiv* preprint (2025).

    > - [DOI:10.1101/2025.07.07.663573](https://doi.org/10.1101/2025.07.07.663573)
    >
    > - PhotoGEA used for: Fitting C₃*A-C_(i)* curves (`fit_c3_aci`)
    >
    > - PhotoGEA version: `v1.1.0`
    >
    > - Analysis scripts: Not yet available
