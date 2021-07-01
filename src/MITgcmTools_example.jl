### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 8cf4d8ca-84eb-11eb-22d2-255ce7237090
begin
	using MITgcmTools, ClimateModels, PlutoUI, Printf
	exps=verification_experiments()
	iexp=findall([exps[i].configuration=="tutorial_global_oce_biogeo" for i in 1:length(exps)])[1]
	myexp=exps[iexp]
	
	build_path=joinpath(MITgcm_path,"verification",myexp.configuration,"build")
	run_path=joinpath(myexp.folder,"run")
	
	🏁 = "🏁"
end

# ╔═╡ 6ef93b0e-859f-11eb-1b3b-d76b26d678dc
begin
	md"""# global ocean biogeochemistry tutorial

	###

	Here we experiment with [MIT general circulation model (MITgcm)](https://mitgcm.readthedocs.io/en/latest/?badge=latest)'s global ocean biogeochemistry tutorial (`tutorial_global_oce_biogeo`) via the [ClimateModels.jl](https://github.com/gaelforget/ClimateModels.jl) interface as implemented in [MITgcmTools.jl](https://gaelforget.github.io/MITgcmTools.jl/dev/).
	"""
end

# ╔═╡ 7fa8a460-89d4-11eb-19bb-bbacdd32719a
md"""## Model Configuration

The model will compile in a subfolder of `MITgcm` :

*$build_path*

And it will run in a temporary folder :

*$run_path*

The model configuration can be summarized as shown below.
	
"""

# ╔═╡ 1f3096d3-ca68-4a71-9411-fe3b201cf5a9
myexp

# ╔═╡ d90039c4-85a1-11eb-0d82-77db4decaa6e
md"""## Workflow Summary

- the `build` method will compile MITgcm (mitgcmuv)
- the `setup` method will set up the run directory
- the `launch` method will start the model run


_Note : the `clean` method can be called to remove a previous run directory before going back to `setup`, `launch`, etc._ 
"""

# ╔═╡ 76291182-86d1-11eb-1524-73dc02ca7b64
@bind do_build Button("Build MITgcm (mitgcmuv)")

# ╔═╡ 848241fe-86d1-11eb-3b30-b94aa0b4431d
let
	do_build
	#build(exps[iexp])
	🏁
end

# ╔═╡ 8569269c-859c-11eb-1ab1-2d874dfa741b
@bind do_cleanup Button("Clean up run/ directory")

# ╔═╡ f008ccaa-859c-11eb-1188-114843d333e6
let
	do_cleanup
	clean(exps[iexp])
	🏁
end

# ╔═╡ 11b024ac-86d1-11eb-1db9-47a5e41398e3
@bind do_link Button("Link files to run/")

# ╔═╡ 31829f08-86d1-11eb-3e26-dfae038b4c01
let
	do_link
	setup(exps[iexp])
	🏁
end

# ╔═╡ 5d826e4c-859d-11eb-133d-859c3abe3ebe
@bind do_run Button("Run model in run/")

# ╔═╡ 550d996a-859d-11eb-34bf-717389fbf809
let
	do_run
	launch(exps[iexp])
	🏁
end

# ╔═╡ 3f58906b-4f50-4df6-84dd-56abe4c9a4c3
md"""## Contents of `run/` directory
"""

# ╔═╡ a04c1cd6-3b9e-4e69-b986-c863b120bb0b
begin
	rundir=joinpath(exps[iexp].folder,string(exps[iexp].ID),"run")
	readdir(rundir)
end

# ╔═╡ Cell order:
# ╟─6ef93b0e-859f-11eb-1b3b-d76b26d678dc
# ╟─8cf4d8ca-84eb-11eb-22d2-255ce7237090
# ╟─7fa8a460-89d4-11eb-19bb-bbacdd32719a
# ╟─1f3096d3-ca68-4a71-9411-fe3b201cf5a9
# ╟─d90039c4-85a1-11eb-0d82-77db4decaa6e
# ╟─76291182-86d1-11eb-1524-73dc02ca7b64
# ╟─848241fe-86d1-11eb-3b30-b94aa0b4431d
# ╟─8569269c-859c-11eb-1ab1-2d874dfa741b
# ╟─f008ccaa-859c-11eb-1188-114843d333e6
# ╟─11b024ac-86d1-11eb-1db9-47a5e41398e3
# ╟─31829f08-86d1-11eb-3e26-dfae038b4c01
# ╟─5d826e4c-859d-11eb-133d-859c3abe3ebe
# ╟─550d996a-859d-11eb-34bf-717389fbf809
# ╟─3f58906b-4f50-4df6-84dd-56abe4c9a4c3
# ╟─a04c1cd6-3b9e-4e69-b986-c863b120bb0b
