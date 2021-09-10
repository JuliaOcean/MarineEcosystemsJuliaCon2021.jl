# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Julia 1.6.0
#     language: julia
#     name: julia-1.6
# ---

# + [markdown] slideshow={"slide_type": "slide"}
# # [JuliaCon 2021](https://juliacon.org/2021/) Workshop
#
# **Title:** Modeling Marine Ecosystems At Multiple Scales Using Julia
#
# **Speakers:** [Gael Forget](https://github.com/gaelforget), [Benoit Pasquier](https://github.com/briochemc), [Zhen Wu](https://github.com/zhenwu0728)
#
# - material : <https://github.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl>
# - streaming : <https://www.youtube.com/watch?v=UCIRrXz2ZS0>
# - Q & A : <https://pigeonhole.at/MARINE>

# + [markdown] slideshow={"slide_type": "subslide"}
# Life in the oceans is strongly connected to our climate. In this workshop, you will learn to use packages from the [JuliaOcean](https://github.com/JuliaOcean) and [JuliaClimate](https://github.com/JuliaClimate) organizations that provide a foundation for studying marine ecosystems across a wide range of scales. We will run agent-based models to explore individual microbes and processes that drive species interactions. On the other end of the model hierarchy, we will simulate planetary-scale transports that control ocean biogeography and climate change.

# + [markdown] slideshow={"slide_type": "slide"}
# ## Schedule
#
# - Introduction of the topics covered, presenters, installation, and workshop roadmap (15 minutes).
#
# - [AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl): concept, implementation, tutorial workthough (30 minutes + 10' for questions; [this notebook URL](https://raw.githubusercontent.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl/main/src/AIBECSExample.jl))
#
# - [PlanktonIndividuals.jl](https://github.com/JuliaOcean/PlanktonIndividuals.jl): concept, implementation, tutorial workthough (30 minutes + 10' for questions; [this notebook URL](https://raw.githubusercontent.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl/main/src/PlanktonIndividualExample.jl))
#
# - [ClimateModels.jl](https://github.com/gaelforget/ClimateModels.jl) and [MITgcmTools.jl](https://github.com/gaelforget/MITgcmTools.jl): concept, implementation, tutorial workthough (20 minutes + 10' for questions; [this notebook URL](https://raw.githubusercontent.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl/main/src/MITgcm_tutorial_global_oce_biogeo.jl))
#
# - [ClimateModels.jl](https://github.com/gaelforget/ClimateModels.jl) and [IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl): concept, implementation, tutorial workthough (20 minutes + 10' for questions; [this notebook URL](https://raw.githubusercontent.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl/main/src/IndividualDisplacementsExample.jl))
#
# - Q&A, tutorials, etc wrap-up
#

# + [markdown] slideshow={"slide_type": "slide"}
# ## Context
#
# - marine ecosystems and climate: 
#     - from microbes O($10^{-6}m$) to whole Earth O($10^6$m); from seconds to millenia
#     - life in a moving fluid; interactions between physics, chemistry, and biology
#     - many fascinating and crucial questions w.r.t. climate change (carbon, food, biodiversity, ...)
# - hierarchy of models :
#     - single ODE to massively parallel HPC simulations
#     - various languages; julia as the future and to unify

# + [markdown] slideshow={"slide_type": "-"} cell_style="split"
# <img src="https://github.com/JuliaOcean/PlanktonIndividuals.jl/raw/master/docs/src/PI_Quota.jpeg" width="90%">
#
#

# + [markdown] slideshow={"slide_type": "-"} cell_style="split"
# <img src="https://mitgcm.readthedocs.io/en/latest/_images/co2flux.png" width="90%">
#
#

# + [markdown] slideshow={"slide_type": "subslide"}
# - Github Organizations : 
#     - JuliaOcean
#     - JuliaClimate
#     - ...
# - Funding / Support : 
#     - NASA (PO, MAP, IDS)
#     - Simons Foundation (CBIOMES, Scope, Gradients)
#     - ...

# + [markdown] slideshow={"slide_type": "slide"}
# ## How To 
#
# To run the notebooks of this workshop on your machine, you need to:
#
# 1. **Install julia** (v1.6.0 or later, <https://julialang.org/>)
#
# 1. **Start julia**
#
# 1. **Add Pluto.jl** (v0.15.0 or later, https://github.com/fonsp/Pluto.jl)
#
#     ```julia
#     import Pkg
#     Pkg.add("Pluto")
#     ```
#
# 1. **Use Pluto to run a notebook**
#
#     ```julia
#     using Pluto
#     Pluto.run(notebook="https://raw.githubusercontent.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl/main/src/AIBECSExample.jl")
#     ```

# + [markdown] slideshow={"slide_type": "subslide"}
# **Alternatively**, instead of your own computer, you can just launch a Pluto instance in the cloud using [JuliaHub.com](https://juliahub.com/ui/Home), paste a notebook URL in the [Pluto start page](https://github.com/fonsp/Pluto.jl), and click open.

# + [markdown] slideshow={"slide_type": "slide"}
# ## Schedule
#
# - Introduction of the topics covered, presenters, installation, and workshop roadmap (15 minutes).
#
# - [AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl): concept, implementation, tutorial workthough (30 minutes + 10' for questions; [this notebook URL](https://raw.githubusercontent.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl/main/src/AIBECSExample.jl))
#
# - [PlanktonIndividuals.jl](https://github.com/JuliaOcean/PlanktonIndividuals.jl): concept, implementation, tutorial workthough (30 minutes + 10' for questions; [this notebook URL](https://raw.githubusercontent.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl/main/src/PlanktonIndividualExample.jl))
#
# - [ClimateModels.jl](https://github.com/gaelforget/ClimateModels.jl) and [MITgcmTools.jl](https://github.com/gaelforget/MITgcmTools.jl): concept, implementation, tutorial workthough (20 minutes + 10' for questions; [this notebook URL](https://raw.githubusercontent.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl/main/src/MITgcm_tutorial_global_oce_biogeo.jl))
#
# - [ClimateModels.jl](https://github.com/gaelforget/ClimateModels.jl) and [IndividualDisplacements.jl](https://github.com/JuliaClimate/IndividualDisplacements.jl): concept, implementation, tutorial workthough (20 minutes + 10' for questions; [this notebook URL](https://raw.githubusercontent.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl/main/src/IndividualDisplacementsExample.jl))
#
# - Q&A, tutorials, etc wrap-up
#