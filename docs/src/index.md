# MarineEcosystemsJuliaCon2021.jl

[JuliaCon 2021](https://juliacon.org/2021/) Workshop

**Title:** Modeling Marine Ecosystems At Multiple Scales Using Julia

**Speakers:** [Gael Forget](https://github.com/gaelforget), [Benoit Pasquier](https://github.com/briochemc), [Zhen Wu](https://github.com/zhenwu0728)

## 2021/07/25 Workshop Recording

streaming : https://www.youtube.com/watch?v=UCIRrXz2ZS0

webpage : https://pretalx.com/juliacon2021/talk/FEZW9Q/

[<img src="https://user-images.githubusercontent.com/20276764/132381907-1ab7d682-ea3d-4db7-b245-3cdb9dd2dcd3.png" width="75%">](https://www.youtube.com/watch?v=UCIRrXz2ZS0)

## Abstract

Life in the oceans is strongly connected to our climate. In this workshop, you will learn to use packages from the [JuliaOcean](https://github.com/JuliaOcean) and [JuliaClimate](https://github.com/JuliaClimate) organizations that provide a foundation for studying marine ecosystems across a wide range of scales. We will run agent-based models to explore individual microbes and processes that drive species interactions. On the other end of the model hierarchy, we will simulate planetary-scale transports that control ocean biogeography and climate change.

## Notebooks

Any example found in the online documentation is most easily run using [Pluto.jl](https://github.com/fonsp/Pluto.jl). Just copy the corresponding `notebook url` link below and paste into the [Pluto.jl interface](https://github.com/fonsp/Pluto.jl/wiki/🔎-Basic-Commands-in-Pluto).

- [AIBECSExample.html](AIBECSExample.html) (---> [notebook url](AIBECSExample.jl))
- [PlanktonIndividualExample.html](PlanktonIndividualExample.html) (---> [notebook url](PlanktonIndividualExample.jl))
- [MITgcm_tutorial_global_oce_biogeo.html](MITgcm_tutorial_global_oce_biogeo.html) (---> [notebook url](MITgcm_tutorial_global_oce_biogeo.jl))
- [IndividualDisplacementsExample.html](IndividualDisplacementsExample.html) (---> [notebook url](IndividualDisplacementsExample.jl))

## Description

Packages covered in this workshop will include:

- [AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl): global steady-state biogeochemistry and gridded transport models that run fast for long time scales (centuries or even millennia).
- [PlanktonIndividuals.jl](https://github.com/JuliaOcean/PlanktonIndividuals.jl): local to global agent-based model, particularly suited to study microbial communities, plankton physiology, and nutrient cycles.
- [MITgcmTools.jl](https://github.com/gaelforget/MITgcmTools.jl): interface to full-featured, Fortran-based, general circulation model and its output (transports, chemistry, ecology, ocean, sea-ice, atmosphere, and more).
- [IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl): local to global particle tracking, for simulating dispersion, connectivity, transports in the ocean or atmosphere, etc.

