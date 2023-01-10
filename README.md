# [JuliaCon 2021](https://juliacon.org/2021/) Workshop

**Title:** Modeling Marine Ecosystems At Multiple Scales Using Julia

**Speakers:** [Gael Forget](https://github.com/gaelforget), [Benoit Pasquier](https://github.com/briochemc), [Zhen Wu](https://github.com/zhenwu0728)

## 2021/07/25 Workshop Recording

streaming : https://www.youtube.com/watch?v=UCIRrXz2ZS0

webpage : https://pretalx.com/juliacon2021/talk/FEZW9Q/

documentation : https://juliaocean.github.io/MarineEcosystemsJuliaCon2021.jl/dev/

[<img src="https://user-images.githubusercontent.com/20276764/132381907-1ab7d682-ea3d-4db7-b245-3cdb9dd2dcd3.png" width="75%">](https://www.youtube.com/watch?v=UCIRrXz2ZS0)

## Abstract

Life in the oceans is strongly connected to our climate. In this workshop, you will learn to use packages from the [JuliaOcean](https://github.com/JuliaOcean) and [JuliaClimate](https://github.com/JuliaClimate) organizations that provide a foundation for studying marine ecosystems across a wide range of scales. We will run agent-based models to explore individual microbes and processes that drive species interactions. On the other end of the model hierarchy, we will simulate planetary-scale transports that control ocean biogeography and climate change.

## Notebooks

Any example found in the online documentation is most easily run using [Pluto.jl](https://github.com/fonsp/Pluto.jl) . 

Just copy the corresponding `notebook url` link below and paste into the [Pluto.jl interface](https://github.com/fonsp/Pluto.jl/wiki/🔎-Basic-Commands-in-Pluto) (v0.15 or later).

- [AIBECSExample.html](https://juliaocean.github.io/MarineEcosystemsJuliaCon2021.jl/dev/AIBECSExample.html) (---> [notebook url](https://juliaocean.github.io/MarineEcosystemsJuliaCon2021.jl/dev/AIBECSExample.jl))
- [PlanktonIndividualExample.html](https://juliaocean.github.io/MarineEcosystemsJuliaCon2021.jl/dev/PlanktonIndividualExample.html) (---> [notebook url](https://juliaocean.github.io/MarineEcosystemsJuliaCon2021.jl/dev/PlanktonIndividualExample.jl))
- [MITgcm_tutorial_global_oce_biogeo.html](https://juliaocean.github.io/MarineEcosystemsJuliaCon2021.jl/dev/MITgcm_tutorial_global_oce_biogeo.html) (---> [notebook url](https://juliaocean.github.io/MarineEcosystemsJuliaCon2021.jl/dev/MITgcm_tutorial_global_oce_biogeo.jl))
- [IndividualDisplacementsExample.html](https://juliaocean.github.io/MarineEcosystemsJuliaCon2021.jl/dev/IndividualDisplacementsExample.html) (---> [notebook url](https://juliaocean.github.io/MarineEcosystemsJuliaCon2021.jl/dev/IndividualDisplacementsExample.jl))

## Description

Packages covered in this workshop include:

- [AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl): global steady-state biogeochemistry and gridded transport models that run fast for long time scales (centuries or even millennia).
- [PlanktonIndividuals.jl](https://github.com/JuliaOcean/PlanktonIndividuals.jl): local to global agent-based model, particularly suited to study microbial communities, plankton physiology, and nutrient cycles.
- [MITgcmTools.jl](https://github.com/gaelforget/MITgcmTools.jl): interface to full-featured, Fortran-based, general circulation model and its output (transports, chemistry, ecology, ocean, sea-ice, atmosphere, and more).
- [IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl): local to global particle tracking, for simulating dispersion, connectivity, transports in the ocean or atmosphere, etc.

The workshop was organized around tutorials and self-contained Pluto notebooks for the different packages.

## Schedule

- Introduction of the topics covered, presenters, installation, and workshop roadmap (15 minutes).

- [AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl): concept, implementation, tutorial workthough (30 minutes + 10' for questions; [this notebook URL](https://raw.githubusercontent.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl/main/src/AIBECSExample.jl))

- [PlanktonIndividuals.jl](https://github.com/JuliaOcean/PlanktonIndividuals.jl): concept, implementation, tutorial workthough (30 minutes + 10' for questions; [this notebook URL](https://raw.githubusercontent.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl/main/src/PlanktonIndividualExample.jl))

- [ClimateModels.jl](https://github.com/gaelforget/ClimateModels.jl) and [MITgcmTools.jl](https://github.com/gaelforget/MITgcmTools.jl): concept, implementation, tutorial workthough (20 minutes + 10' for questions; [this notebook URL](https://raw.githubusercontent.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl/main/src/MITgcm_tutorial_global_oce_biogeo.jl))

- [ClimateModels.jl](https://github.com/gaelforget/ClimateModels.jl) and [IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl): concept, implementation, tutorial workthough (20 minutes + 10' for questions; [this notebook URL](https://raw.githubusercontent.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl/main/src/IndividualDisplacementsExample.jl))

- Q&A, tutorials, etc wrap-up

Workshop materials are available ahead of time @ https://github.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl



## Setup instructions


To run the notebooks of this workshop on your machine, you need to:

1. **Install Julia** from <https://julialang.org/> (latest version is v1.6.2)

1. **Start Julia**

1. **Add Pluto.jl** (v0.15.0 or later)**

    This is simply done by typing, in the julia REPL,

    ```julia
    import Pkg
    Pkg.add("Pluto")
    ```

    > *Note*: Please make sure you get version 0.15.0 or later.
    > If you get an older version then you can add Pluto in a clean, temporary, environment as follows:
    > ```julia
    > import Pkg
    > Pkg.activate(mktempdir())
    > Pkg.add("Pluto")
    > ```

1. **Use Pluto to run the notebooks.**
    This is as simple as copy-pasting one of the following lines, depending on which notebook you want to run:

    ```julia
    using Pluto
    Pluto.run(notebook="https://raw.githubusercontent.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl/main/src/AIBECSExample.jl")
    ```
**Alternatively**, instead of your own computer, you can just launch a Pluto instance in the cloud using [JuliaHub.com](https://juliahub.com/ui/Home), paste a notebook URL in the [Pluto start page](https://github.com/fonsp/Pluto.jl), and click open.
