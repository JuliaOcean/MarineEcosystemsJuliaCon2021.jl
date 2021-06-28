### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 153caed2-f735-4508-8e80-f0461a6af88c
using PlanktonIndividuals, Plots

# ╔═╡ 418c95a7-0125-40b6-acb5-e2011ea60df2
md"""
# Vertical 2-Dimensional Example
"""


# ╔═╡ 6ec492fc-d622-11eb-2bfb-9d24a91685b7
md"""
Here we simulate phytoplankton cells as Lagrangian particles in a 2D flow field, with 
one horizontal direction (x) and one vertical one (z), like in an ocean transect.

The domain is periodic in the x direction while it is bounded in the z direction.
"""

# ╔═╡ 139108e4-0829-4710-8419-dc05808f173f
md"""
## 1. Import packages
"""

# ╔═╡ 15bc275c-7d15-479c-b644-d588e67cf3f5
md"""
Full documentation is available at [https://juliaocean.github.io/PlanktonIndividuals.jl/dev/](https://juliaocean.github.io/PlanktonIndividuals.jl/dev/)
"""

# ╔═╡ 3c05c373-af76-4d28-979c-e91aaa7321de
p=dirname(pathof(PlanktonIndividuals));

# ╔═╡ b7f4cd3e-78aa-44cd-bb6c-3f41ed9f8a96
include(joinpath(p,"../examples/helper_functions.jl"));

# ╔═╡ 37f42666-276f-43bf-b4b4-a245af7edb29
md"""
## 2.Grid & Flow Fields
"""

# ╔═╡ 49c16379-d570-4849-bb74-0b7ccedb86a5
md"""
First we generate grid information (128 by 128 grid boxes, 1m thick, and 1m wide) and the computational architecture (CPU).
"""

# ╔═╡ 6e45867e-e12e-4e44-9795-0a9ee01ecccb
arch = CPU()

# ╔═╡ 3b1861b8-4653-4c68-990d-7059dc187518
grid = RegularRectilinearGrid(size=(128, 1, 128), spacing=(1, 1, 1))

# ╔═╡ bee232fe-b10b-4ce8-81b6-6053611bd2e9
md"""
Then we use a stream function (see helper_functions.jl) to generate a simple flow field (displayed below) in a 2D vertical plane.
"""

# ╔═╡ bc16d909-f629-42b9-952b-777470b00a52
(uvels, vvels, wvels, ϕcenters) = streamfunction_xz();

# ╔═╡ 5f5a05c6-c76d-4f57-848b-73e391abe7dd
md"""
## 3. Model Setup
"""

# ╔═╡ 097232e8-dbe7-4c0e-8a3e-25638c28a623
md"""
Next we setup the individual-based model by specifying the architecture, grid, and plankton community.
"""

# ╔═╡ 467d1cfd-aadc-4e16-86d4-94cabedcd15a
model = PlanktonModel(arch, grid; N_species = 1, N_individual = 2^7, max_individuals = 2^7*8)

# ╔═╡ d4dd3b10-8953-42c9-9274-50b201e08a95
md"""
Finally we setup the duration of the model simulation and the kind of output we want.
"""

# ╔═╡ 5931e9f6-c028-4933-9f63-2c8881221f25
sim = PlanktonSimulation(model, ΔT = 60, nΔT = 1, vels=(u=uvels, v=vvels, w=wvels), ΔT_vel=60*120)

# ╔═╡ 7c58dd75-c636-4c02-8647-38bae5e8ab6a
md"""
## 4. Model Run
"""

# ╔═╡ 8aa211b9-cdb8-4504-b859-e768d6f842dd
md"""
We run the model for 120 time steps (2 hours) and then plot individuals and nutrients in their final state (stored in model).
"""

# ╔═╡ c79a07b2-e448-44d8-b0a0-c42b5fc1db89
for i in 1:120
    update!(sim)
end

# ╔═╡ ce973a4d-b612-450a-8dd9-ddcfbec59c65
md"""
To plot the distribution of individuals as well as nutrient fields we use Plots.jl and create a function that can easily be re-used e.g. to create an animation.
"""

# ╔═╡ d87581c9-63ec-4519-8c84-f8442f05e0a8
function plot_model(model::PlanktonModel)
    ## Coordinate arrays for plotting
    xC, zC = collect(model.grid.xC)[3:130], collect(model.grid.zC)[3:130]

    ## contour of the flow field
    fl_plot = Plots.contourf(xC, reverse(zC), rotl90(ϕcenters), xlabel="x (m)", ylabel="z (m)", color=:balance, fmt=:png, colorbar=false)

    ## a scatter plot embeded in the flow fields
    px = Array(model.individuals.phytos.sp1.data.x) .* 1 # convert fractional indices to degree
    pz = Array(model.individuals.phytos.sp1.data.z) .* -1# convert fractional indices to degree
    Plots.scatter!(fl_plot, px, pz, ms=5, color = :red, legend=:none)

    ## DOC field
    trac1 = Plots.contourf(xC, reverse(zC), rotl90(Array(model.nutrients.DOC.data)[3:130,3,3:130]), xlabel="x (m)", ylabel="z (m)", clims=(0.5, 1.1), fmt=:png)

    ## Arrange the plots side-by-side.
    plt = Plots.plot(fl_plot, trac1, size=(800, 400),
        title=[lpad(model.t÷86400,2,"0")*"day "*lpad(model.t÷3600-24*(model.t÷86400),2,"0")*"hour" "DOC (mmolC/L)"])

    return plt
end

# ╔═╡ 581809a1-6f8b-40d0-9a6e-d943c424be07
plot_model(model)

# ╔═╡ 3073ceb2-dccb-400e-ae0d-8ff73569cd54
md"""
You can also generate a animation if you like.
"""

# ╔═╡ 2c119e21-afbf-48ae-aebe-fc161c18f344
anim = @animate for i in 1:120
	update!(sim)
   plot_model(model)
end

# ╔═╡ a685d6e7-88aa-427b-be61-f1e5da15653a
gif(anim, "anim_fps15.gif", fps = 15)

# ╔═╡ fa502a3f-ced0-4e82-afce-3baf4cacf4c7
md"""
## 5. Exercise
"""

# ╔═╡ af8be494-afdf-4e2e-a1be-b069c3a14f18
md"""
1. Add the second species
2. Run for a longer period
3. Add diagnostics
4. ...
"""

# ╔═╡ 964cb129-9d7c-47bc-945a-7639e80da21f
md"""
#### 1. Add the second species

Run the code below with `N_species = 2` to add the second species.

```julia
model = PlanktonModel(arch, grid; N_species = 2, N_individual = 2^7, max_individuals = 2^7*8)
```

In this model setup, there're **two** species, each species has `2^7` individuals.
"""

# ╔═╡ ed092d52-8afb-4c97-a1fe-4e6cbae33420
md"""
#### 2. Run for a longer period

In order to run a longer period, you need first to change `ΔT_vel` to tell the model to use the flow fields for more time steps.

For example, if we want to run the model for 1 day:

```julia
sim = PlanktonSimulation(model, ΔT = 60, nΔT = 60*24, vels=(u=uvels, v=vvels, w=wvels), ΔT_vel=60*60*24)
```

As you may noticed, there're two ways to update the model for multiple time steps.

The first one is that `nΔT` is set to `1`, and we use a external for loop to update the model. In this way, you are able to generate animations.

The second one is that `nΔT` is set to the number of time steps you want to run (see example above). In this way, you don't need external for loops.
"""

# ╔═╡ db9910c4-b5a9-44a9-9596-d5003bcb8c24
md"""
#### 3. Add diagnostics

Diagnostics is a good way to check if the model is run correctly.
Diagnoscits in the model include two parts: Eulerian tracers and individuals.
Diagnostics of individuals are aggregated into gridded fields.

A full list of available diagnostics can be found [here](https://juliaocean.github.io/PlanktonIndividuals.jl/dev/model_run/#Model-Diagnostics)

Below is a example of model diagnostics. In this example, we are interested in photosynthetically active radiation(`:PAR`), ammonia(`:NH4`), nitrate(`:NO3`), dissolved organic carbon(`:DOC`) as well as the distribution of individuals(`:num`), grazing rate(`:graz`), mortality rate(`:mort`), and division rate(`:dvid`).

```julia
diags = PlanktonDiagnostics(model; tracer=(:PAR, :NH4, :NO3, :DOC),
                    plankton=(:num, :graz, :mort, :dvid),
                    frequency = 1)
```

`frequency` indicates the frequency of diagnostics (in numbers of time steps), diagnose every time step by default.

Then we need to include this `PlanktonDiagnostics` in `PlanktonSimulation`.

```julia
sim = PlanktonSimulation(model, ΔT = 60, nΔT = 60*24,
						 diags = diags,
						 vels=(u=uvels, v=vvels, w=wvels),
						 ΔT_vel=60*120)
```
"""

# ╔═╡ 9c8ae4fe-b682-4afc-9c76-0ac3a4b7af7b


# ╔═╡ Cell order:
# ╟─418c95a7-0125-40b6-acb5-e2011ea60df2
# ╟─6ec492fc-d622-11eb-2bfb-9d24a91685b7
# ╟─139108e4-0829-4710-8419-dc05808f173f
# ╠═153caed2-f735-4508-8e80-f0461a6af88c
# ╟─15bc275c-7d15-479c-b644-d588e67cf3f5
# ╠═3c05c373-af76-4d28-979c-e91aaa7321de
# ╠═b7f4cd3e-78aa-44cd-bb6c-3f41ed9f8a96
# ╟─37f42666-276f-43bf-b4b4-a245af7edb29
# ╟─49c16379-d570-4849-bb74-0b7ccedb86a5
# ╠═6e45867e-e12e-4e44-9795-0a9ee01ecccb
# ╠═3b1861b8-4653-4c68-990d-7059dc187518
# ╟─bee232fe-b10b-4ce8-81b6-6053611bd2e9
# ╠═bc16d909-f629-42b9-952b-777470b00a52
# ╟─5f5a05c6-c76d-4f57-848b-73e391abe7dd
# ╟─097232e8-dbe7-4c0e-8a3e-25638c28a623
# ╠═467d1cfd-aadc-4e16-86d4-94cabedcd15a
# ╟─d4dd3b10-8953-42c9-9274-50b201e08a95
# ╠═5931e9f6-c028-4933-9f63-2c8881221f25
# ╟─7c58dd75-c636-4c02-8647-38bae5e8ab6a
# ╟─8aa211b9-cdb8-4504-b859-e768d6f842dd
# ╠═c79a07b2-e448-44d8-b0a0-c42b5fc1db89
# ╟─ce973a4d-b612-450a-8dd9-ddcfbec59c65
# ╟─d87581c9-63ec-4519-8c84-f8442f05e0a8
# ╠═581809a1-6f8b-40d0-9a6e-d943c424be07
# ╟─3073ceb2-dccb-400e-ae0d-8ff73569cd54
# ╠═2c119e21-afbf-48ae-aebe-fc161c18f344
# ╠═a685d6e7-88aa-427b-be61-f1e5da15653a
# ╟─fa502a3f-ced0-4e82-afce-3baf4cacf4c7
# ╟─af8be494-afdf-4e2e-a1be-b069c3a14f18
# ╟─964cb129-9d7c-47bc-945a-7639e80da21f
# ╟─ed092d52-8afb-4c97-a1fe-4e6cbae33420
# ╟─db9910c4-b5a9-44a9-9596-d5003bcb8c24
# ╠═9c8ae4fe-b682-4afc-9c76-0ac3a4b7af7b
