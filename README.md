# JuliaCon 2021 Workshop 

**Title:** Modeling Marine Ecosystems At Multiple Scales Using Julia

**Speakers:** Gael Forget, Benoit Pasquier, Zhen Wu

## Abstract

Life in the oceans is strongly connected to our climate. In this workshop, you will learn to use packages from the JuliaOcean and JuliaClimate organizations that provide a foundation for studying marine ecosystems across a wide range of scales. We will first run agent-based models to explore individual microbes and processes that drive species interactions. On the other end of the model hierarchy, we will simulate planetary-scale transports that control ocean biogeography and climate change.

## Description

Packages covered in this workshop will include: 

- `AIBECS.jl` : global steady-state biogeochemistry and gridded transport models that run fast for long time scales (centuries or even millenia).
- `PlanktonIndividuals.jl` : local to global agent based model, particluarly suited to study microbial communities, plankton physiology, and nutrient cycles.
- `IndividualDisplacements.jl` : local to global particle tracking, for simulating dispersion, connectivity, transports in the ocean or atmosphere, etc.
- `MITgcmTools.jl` : interface to full-featured, fortran-based, general circulation model and its output (transports, chemistry, ecology, ocean, seaice, atmosphere, and more).

The workshop's first two hours will be organized around tutorials and self-contained Pluto notebooks for the different packages.

The third hour will provide the opportunity for attendees to further explore the models in breakout rooms and via exercises.

### Workshop schedule

- Introduction of the topics covered, presenters, installation, and workshop roadmap (15 minutes).

- `AIBECS.jl` : concept, implementation, tutorial workthough (20 minutes + 10' for questions)

- `PlanktonIndividuals.jl` : concept, implementation, tutorial workthough (20 minutes + 10' for questions)

- `IndividualDisplacements.jl` : concept, implementation, tutorial workthough (10 minutes + 10' for questions)

- `MITgcmTools.jl` : concept, implementation, tutorial workthough (10 minutes + 10' for questions)

- 5-minute break

- Breakout rooms for deeper dive in tutorials, exercises, or trying out your own idea with guidance from the presenters (1 hour)

Workshop materials will be made available ahead of time @ https://github.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl