# [JuliaCon 2021](https://juliacon.org/2021/) Workshop 

**Title:** Modeling Marine Ecosystems At Multiple Scales Using Julia

**Speakers:** [Gael Forget](https://github.com/gaelforget), [Benoit Pasquier](https://github.com/briochemc), [Zhen Wu](https://github.com/zhenwu0728)






## Abstract

Life in the oceans is strongly connected to our climate. In this workshop, you will learn to use packages from the [JuliaOcean](https://github.com/JuliaOcean) and [JuliaClimate](https://github.com/JuliaClimate) organizations that provide a foundation for studying marine ecosystems across a wide range of scales. We will first run agent-based models to explore individual microbes and processes that drive species interactions. On the other end of the model hierarchy, we will simulate planetary-scale transports that control ocean biogeography and climate change.






## Description

Packages covered in this workshop will include: 

- [AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl): global steady-state biogeochemistry and gridded transport models that run fast for long time scales (centuries or even millennia).
- [PlanktonIndividuals.jl](https://github.com/JuliaOcean/PlanktonIndividuals.jl): local to global agent-based model, particularly suited to study microbial communities, plankton physiology, and nutrient cycles.
- [IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl): local to global particle tracking, for simulating dispersion, connectivity, transports in the ocean or atmosphere, etc.
- [MITgcmTools.jl](https://github.com/gaelforget/MITgcmTools.jl): interface to full-featured, Fortran-based, general circulation model and its output (transports, chemistry, ecology, ocean, sea-ice, atmosphere, and more).

The workshop's first two hours will be organized around tutorials and self-contained Pluto notebooks for the different packages.

The third hour will provide the opportunity for attendees to further explore the models in breakout rooms and via exercises.





## Schedule

- Introduction of the topics covered, presenters, installation, and workshop roadmap (15 minutes).

- [AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl): concept, implementation, tutorial workthough (20 minutes + 10' for questions)

- [PlanktonIndividuals.jl](https://github.com/JuliaOcean/PlanktonIndividuals.jl): concept, implementation, tutorial workthough (20 minutes + 10' for questions)

- [IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl): concept, implementation, tutorial workthough (10 minutes + 10' for questions)

- [MITgcmTools.jl](https://github.com/gaelforget/MITgcmTools.jl): concept, implementation, tutorial workthough (10 minutes + 10' for questions)

- 5-minute break

- Breakout rooms for deeper dive in tutorials, exercises, or trying out your own idea with guidance from the presenters (1 hour)

Workshop materials will be made available ahead of time @ https://github.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl





## Setup instructions

### Run the notebooks on your local machine

To run the notebooks of this workshop on your machine, you need to:

1. **Install Julia v1.6.1** from https://julialang.org/.

1. **Start Julia.**

1. **Add the Pluto package (v0.14.8 or later).**
    Making sure you get the latest version can easily be done by making sure Pluto is added in a clean temporary environment. If you are unsure, just copy-paste these lines in your Julia REPL:

    ```julia
    import Pkg
    Pkg.activate(mktempdir())
    Pkg.add("Pluto")
    ```

1. **Use Pluto to run the notebooks.**
    This is as simple as copy-pasting one of the following lines, depending on which notebook you want to run:

    ```julia
    using Pluto
    Pluto.run(notebook="https://raw.githubusercontent.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl/main/src/AIBECSExample.jl")
    ```

### Run the notebooks remotely

TBC.
