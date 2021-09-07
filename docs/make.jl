using Documenter, MarineEcosystemsJuliaCon2021
import PlutoSliderServer
using Plots

makedocs(;
    modules=[MarineEcosystemsJuliaCon2021],
    format=Documenter.HTML(),
    repo="https://github.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl/blob/{commit}{path}#L{line}",
    sitename="MarineEcosystemsJuliaCon2021.jl",
    authors="gaelforget <gforget@mit.edu>",
)

pth = joinpath(@__DIR__, "build")
lst=("AIBECSExample.jl","PlanktonIndividualExample.jl","MITgcm_tutorial_global_oce_biogeo.jl","IndividualDisplacementsExample.jl")
for i in lst
    fil_in=joinpath(@__DIR__,"..","src",i)
    fil_out=joinpath(pth,i[1:end-2]*"html")
    PlutoSliderServer.export_notebook(fil_in)
    mv(fil_in[1:end-2]*"html",fil_out)
end

deploydocs(;
    repo="github.com/JuliaOcean/MarineEcosystemsJuliaCon2021.jl",
)
