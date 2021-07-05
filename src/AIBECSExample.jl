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

# ╔═╡ 0df28961-4e34-496b-86c0-a4740b368127
begin
	# create temporary environment
	import Pkg
	Pkg.activate(mktempdir())
	# Add package specific versions for reproducibility
	# (This should work in the future unless one of the dependencies broke semver)
	Pkg.add(name="AIBECS", version="0.9.3")
	Pkg.add(name="Plots", version="1.16.8")
	Pkg.add(name="PlutoUI", version="0.7.9")
	Pkg.add(name="OceanBasins", version="0.1.7")
end

# ╔═╡ 0a889312-d965-11eb-28e0-a72a1379850a
using AIBECS

# ╔═╡ 4d57848c-a2e7-4d5a-b848-f3822437da6a
using Plots

# ╔═╡ 361bbe46-bf9a-4a6c-834f-5a0660f9a0e6
using OceanBasins

# ╔═╡ 325084cb-b233-4a63-8f8c-7ba32b7e8d43
using PlutoUI

# ╔═╡ c91da4d9-a3d1-4c38-aadb-c4548b6cc0e3
md"""
# [AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl) tutorial: Radiocarbon age

In this notebook, you will simulate the radiocarbon age in the global ocean in just a few minutes!

Note that this notebook is heavily inspired from [the radiocarbon tutorial in AIBECS.jl's documentation](https://juliaocean.github.io/AIBECS.jl/stable/tutorials/2_radiocarbon/).
"""

# ╔═╡ 41e825fa-43e2-4b33-9e25-98145cb4eb63
md"""
## Introduction

Radiocarbon, ¹⁴C, is produced by cosmic rays hitting ¹⁴N atoms in the lower stratosphere and upper troposphere.
¹⁴C then quickly reacts with oxygen to produce ¹⁴CO₂, which eventually enters the ocean through air–sea gas exchange.

![](https://wserv4.esc.cam.ac.uk/pastclimate/wp-content/uploads/2014/09/Radiocarbon-cycle_2-e1410282981969.jpg)
*image taken from [a "Past Oceans and Climate" blog post by Luke Skinner](https://wserv4.esc.cam.ac.uk/pastclimate/?page_id=19)*

As it travels the oceans, radiocarbon decays (halflife ~5730 yr).
Deviations of ¹⁴C concentrations away from atmoshperic values serve as a tracer label for the time passed since a water parcel was last in contact with the atmosphere.
We call this mean time the "radiocarbon age", and we are going to calculate it!
"""

# ╔═╡ 1d5df523-879e-4398-849e-d4b21003dd12
html"""
<h2>Tracer equation</h2>
"""

# ╔═╡ fc00184b-9bcc-4833-9fc8-49ef2efbd3af
md"""
The tracer equation we solve the steady-state of is

$$\frac{\partial \boldsymbol{R}}{\partial t} + \mathbf{T} \, \boldsymbol{R} = \frac{\lambda}{h} (R_\mathsf{atm} - \boldsymbol{R}) (\boldsymbol{z} ≤ h) - \boldsymbol{R} / \tau.$$

where now $\boldsymbol{R}$ is a column vector of radiocarbon concentrations and $\mathbf{T}$ is a sparse matrix representing the ocean circulation.
"""

# ╔═╡ 2b531ea8-87eb-4a63-aabc-e7520f18a893
md"""
#### Continuous PDE

The "real" three-dimensional PDE for $R(x,y,z)$ is

$$\frac{\partial R}{\partial t} + \underbrace{\nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \cdot \nabla \right] R}_\mathsf{circulation} = \underbrace{ \frac{\lambda}{h} (R_\mathsf{atm} - R) (z ≤ h)}_\mathsf{air–sea exchange} - \underbrace{R / \tau}_\mathsf{decay},$$

where $\frac{\partial R}{\partial t}$ is the rate of change of $R$ and $\nabla \cdot \left[ \boldsymbol{u} - \mathbf{K} \cdot \nabla \right] R$ is the flux divergence of $R$ due to the ocean circulation, here given by the mean advective currents $\boldsymbol{u}$ and an eddy-diffusivity matrix $\mathbf{K}$. The first term on the right-hand-side represents the air–sea gas exchange with a piston velocity $λ$ over a depth $h$ and $R / \tau$ is the radioactive decay rate, where $\tau$ is the radioactive decay timescale.

In this notebook, the goal is to solve for the steady-state of the tracer equation.
"""

# ╔═╡ 59304232-cf0f-47ec-b993-86afcd2cc227
html"""

<h4> Discretize and linearize the PDE</h4>

<p>For any model of ocean or atmospheric processes, we must discretize the equation onto a 3D grid. Below is an example of a very coarse ocean circulation model, with just 5 boxes for water"</p>

<img src="https://user-images.githubusercontent.com/4486578/58314610-3b130b80-7e53-11e9-9fe8-9527cdcca2d0.png" width=80%><br>

<em>image from François Primeau and Louis Primeau</em>
"""

# ╔═╡ 5473f343-4dd5-4aef-88c9-4287a876d7f9
html"""

<img src="https://user-images.githubusercontent.com/4486578/61757212-fe3ba480-ae02-11e9-8d17-d72866eaafb5.gif" width=60% align="right">
<p>With AIBECS, the 3D field is rearranged into a column vector. Taking the toy model of ocean circulation above, this means that the 5 boxes are stacked together to form a 5-element column, like illustrated on the right.</p>

<p>The tracer equation given after discretization and reorganization into column vectors, is the "discrete" tracer equation. It is given by</p>
"""

# ╔═╡ ba54f82c-9e6c-409d-86bd-a85cc4609357
md"""
## Simulate radiocarbon age
"""

# ╔═╡ 38129859-df5e-4faa-a75c-8304b2bed55b
md"""
##### 1. Load AIBECS

Note I recommend using Pluto's main branch (v0.15) to get self-contained package management. Otherwise, you'll either need to manually handle the packages dependencies or add them to your default environment.
"""

# ╔═╡ 56bf3614-092d-4891-b1aa-7b2c5e1ddcbc
md"""
##### 2. Chose a circulation model
"""

# ╔═╡ 4bbffe2c-0520-4e7b-bf38-0857454336ef
md"""
Load the circulation grid and matrix. Here I chose the `OCCA` matrix, which was built using the MIT general circulation model, but you can select other circulations, like the `OCIM2` one, for example (FYI, the OCCA matrix was made by Gael Forget, the organizer of this workshop).
"""

# ╔═╡ 35641e0f-537c-476a-86bd-d690d0989736
grd, T = OCCA.load() ;

# ╔═╡ bcc5fa0f-7595-4a21-a057-1cbabf94a049
md"""
Note that you can also test the different OCIM2 matrices using the `version` keyword argument, .g., with `version="KiHIGH_He"` or `version="KvHIGH_He"` for increased isopycnal and vertical diffusivities, respectively)
"""

# ╔═╡ d76ec9fb-77f4-4a14-afb2-713cbd03ffab
md"""
The grid, `grd`, contains a bunch of properties, like latitude, longitude, depth, and so on.
"""

# ╔═╡ ba5fb3fa-b9f8-43c3-881c-0dc26c399b70
fieldnames(typeof(grd))

# ╔═╡ 5862c0fe-ff2c-45c9-aef4-ab0635442ee5
md"""
And `T` is the sparse matrix of the ocean circulation.
"""

# ╔═╡ ac71b40d-5f79-4881-9086-c71c89077499
T

# ╔═╡ fe54d279-6d03-4892-9c0f-9dccb42b77a1
md"""
##### 3. Local sources and sinks

We lump the local sources and sinks (air–sea exchange and decay) into a `RHS` function.
"""

# ╔═╡ 3c300825-ef08-4e5d-ab29-9264bb8e1024
md"""
Above you'll note that I "unpack" a bunch of scalar parameters. We'll deal with those in the next section. But first, let us define `z`, the column vector of all the depths of all the boxes, for which there is the `depthvec` function.
"""

# ╔═╡ 5b99fcb5-c3e1-4a01-84e6-334b26562301
z = depthvec(grd)

# ╔═╡ daef9e51-08b5-4ab3-946a-703052ef1a34
function RHS(R,p)
    @unpack λ, h, Ratm, τ = p
    return @. λ / h * (Ratm - R) * (z ≤ h) - R / τ
end

# ╔═╡ c9f07270-f1ce-4b1b-b8b0-a4ab53440366
md"""
##### 4. Parameters

We define our variable parameters (those that I unpacked earlier) with a special `struct`. Here we do that with a little add-on package for units. But to do that we must first import some unit functions
"""

# ╔═╡ 021a82d2-4a7b-4768-9e49-bdcb2a2b5da4
import AIBECS: @units, units

# ╔═╡ 7679b0ba-5a04-4af3-9f58-80644da08f9e
md"""
Now we define our parameter structure called `Params`.
"""

# ╔═╡ 202247cd-7a29-4f24-a5f9-831badce334f
@units struct Params{U} <: AbstractParameters{U}
    λ::U    | u"m/yr"
    h::U    | u"m"
    τ::U    | u"yr"
    Ratm::U | u"M"
end

# ╔═╡ 55292329-200f-4cf3-bb1d-e4cff6def15b
md"""
For the air–sea gas exchange, we use a constant piston velocity $\lambda$ of 50m / 10years.
And for the radioactive decay we use a timescale $\tau$ of 5730/log(2) years.
"""

# ╔═╡ e533c789-73e5-4843-b58e-9b00c950aa0b
p = Params(λ = 50u"m"/10u"yr",
           h = grd.δdepth[1],
           τ = 5730u"yr"/log(2),
           Ratm = 42.0u"nM")

# ╔═╡ fa2265a6-3bf5-4b29-9fcb-271d1215ea62
md"""
##### 5. The "state" function and its Jacobian

To simulate radiocarbon, AIBECS first collects all the terms of the rate of change into a big "state" function that depends on the state (here just the radiocarbon concentration vector $\boldsymbol{R}$ and the parameters $\boldsymbol{p}$),

$$\boldsymbol{F}(\boldsymbol{R}, \boldsymbol{p}) = \frac{\partial \boldsymbol{R}}{\partial t}$$

and then AIBECS solves for the steady-state by looking for roots of $\boldsymbol{F}$, i.e., the $\boldsymbol{R}$ values such that 

$$\boldsymbol{F}(\boldsymbol{R}, \boldsymbol{p}) = 0$$

To find that root, AIBECS uses a quasi-Newton solver, which requires the Jacobian of $\boldsymbol{F}$ w.r.t. $\boldsymbol{R}$, which we denote by $\nabla_\boldsymbol{x}\boldsymbol{F}$, which is a sparse matrix just like $\mathbf{T}$. But don't worry, AIBECS can compute both of them for you using the `F_and_∇ₓF` function.
"""

# ╔═╡ 9a40e63d-bcb0-47d9-9daa-f98e6b2d90d7
F, ∇ₓF = F_and_∇ₓF(T, RHS)

# ╔═╡ 4e82b6b7-abb8-4248-ad92-65e082c2eedc
md"""
##### 6. Run the simulation!

We chose an initial guess for the values of $\boldsymbol{R}$ (here just zeros).
"""

# ╔═╡ 8875b4e6-95ca-4c31-a6c5-70fd7193f5a5
R₀ = zeros(length(z)) # an initial guess

# ╔═╡ 77d30cef-dc91-4db1-886f-15b4787f1ec7
md"""
We generate the steady-state problem (using SciML lingo)"""

# ╔═╡ fa03c28c-5c2a-4764-9260-6acd9b5ec5a3
prob = SteadyStateProblem(F, ∇ₓF, R₀, p) ;

# ╔═╡ da45685f-db9b-4c8a-a62b-0bc9cfc04311
md"""And we solve it."""

# ╔═╡ d3074335-0625-402a-b4d3-6ae65c8cd168
R = solve(prob, CTKAlg()).u # in SciML lingo, the solution array is in `.u`

# ╔═╡ 5ef6f676-ead9-486d-afdd-b145d48e5a24
md"""
Note that in order to get back to the radiocarbon "age", we must compare $R$ to atmospheric values. We do this using the parameters $\tau$ and $R_\mathsf{atm}$, and, while we're at it, we convert it to units of years. (The default unit for time is seconds in AIBECS.)
"""

# ╔═╡ 59dd2d1e-78e2-431d-ad92-a7368b8e54ce
C14age = let
	@unpack τ, Ratm = p
	@. log(Ratm / R) * τ * u"s" |> u"yr"
end

# ╔═╡ fdc77b5c-5a56-43bc-b304-89f1f064c9a5
md"""
That's it, now let's plot this radiocarbon age!
"""

# ╔═╡ 5d8b2518-66a5-4ab2-992b-b80b98d1a498
md"""
### Plot it from all angles
"""

# ╔═╡ 0268bff0-7ccd-404f-8424-78444885fb4a
md"""
We can take a zonal average as well (using the `plotzonalaverage` function)
"""

# ╔═╡ 9a95c55a-8e44-42d3-a52e-2f18d670811b
md"""
Or a mean profile
"""

# ╔═╡ a6407c69-ac39-4819-9429-239f43a5dc63
md"20 minutes?"

# ╔═╡ 56a1aa18-0134-4e4f-a6cc-3662316284dc
md"""
## Extras for the breakout session! 

### Compose with other packages
"""

# ╔═╡ 3db99f8c-f238-427c-8239-19a7d94401b7
md"""
As an example, let's load the polygons of the 3 major basins and look at plots restriced to these areas
"""

# ╔═╡ e0267345-df37-410b-bcca-218e35388f87
OCEANS = oceanpolygons()

# ╔═╡ 246063c6-d2ca-4796-88c5-fa1489e6381b
masks = (
	ATL = isatlantic(latvec(grd), lonvec(grd), OCEANS),
	PAC = ispacific(latvec(grd), lonvec(grd), OCEANS),
	IND = isindian(latvec(grd), lonvec(grd), OCEANS),
)

# ╔═╡ e14f6526-c5e3-4c17-896e-2983c030f5e9
plotverticalmean(sum(k*mask for (k,mask) in enumerate(masks)), grd; title="Ocean basins: $([(k,m) for (k,m) in enumerate(keys(masks))])", color=cgrad(cgrad(:tab10)[1:length(masks)+1], categorical=true), clim=(-0.5,length(masks)+0.5))

# ╔═╡ f7ad9931-1d92-4d08-a0ee-9393dd524f92
md"""
## Extra settings

Below are cells that are not necessarily needed for this notebook but they facilitate some things.
"""

# ╔═╡ 25f7c2a3-a8a7-441e-8da1-b311fd355867
md"""

##### Always accept downloads

The AIBECS.jl package provides functions to download different ocean circulations using the DataDeps.jl package, which ususally requires the user to accept each download individually. The cell below tells DataDeps.jl to always accept downloads, so that this notebook can be self-contained and fetch the circulation you want automatically.
"""

# ╔═╡ c08308c9-eade-4209-b755-2415542e542a
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

# ╔═╡ cffb9e34-c2ab-40c2-8006-e1fd934fd066
md"""
##### Plot options
"""

# ╔═╡ d5a447b4-2d22-4466-b89f-de859ed9f09e
clim=(0,2500)

# ╔═╡ 9f2ece88-dbc7-495e-9ca5-6518af58caf8
# repeated plot options are colorbartitle, colormap, extra margin for colorbar title
plot_options = (;
	colorbartitle="\nage",  # The new line moves cbar label to avoid overlap but...
	right_margin=3Plots.mm, # but the moved cbar label needs more space!
	color=:viridis, 
	clim)

# ╔═╡ d4521113-4a5c-4154-877f-7228885e0ac2
function myhorizontalslice(x, grd; depth, kwarhs...)
	# Default heatmap
	plothorizontalslice(x, grd; depth, title="Radiocarbon age at $(depth)m", plot_options...)
	# Overlay filled contour and some black contours (for looks)
	plothorizontalslice!(x, grd; depth, plot_options..., levs=0:100:2500, st=:contourf, lw=0)
	plothorizontalslice!(x, grd; depth, plot_options..., levs=0:500:2500, st=:contour, c=:black, clabels=true)
end

# ╔═╡ 75b087ff-fb90-4a79-b95b-678157f4046e
function myzonalaverage(x, grd; kwargs...)
	p = plotzonalaverage(x, grd; plot_options..., yunit=u"km", kwargs...)
	plotzonalaverage!(p, x, grd; plot_options..., yunit=u"km", st=:contour, levs=0:500:3000, c=:black, clabels=true, kwargs...)
	p
end

# ╔═╡ f0d1a772-2b1d-4feb-8d52-b9dbf29b3f45
myzonalaverage(C14age, grd; title="Global zonal average of radiocarbon age")

# ╔═╡ 8f76c559-ef89-46c0-a792-75a1de96f6e5
plot([myzonalaverage(C14age, grd; mask, title="Zonal average ($m)") 
		for (m,mask) in pairs(masks)]..., layout=(length(masks),1), size=(500,800))

# ╔═╡ ee1d8881-4dce-45d3-95fa-9bb33819ebe9
plothorizontalaverage(C14age, grd; plot_options..., xlim=clim, lab="", xlab="Horizontally averaged radiocarbon age")

# ╔═╡ 98508321-fb4b-4076-9c95-0d38f63303b8
let
	p = plot()
	for (color, (m,mask)) in enumerate(pairs(masks))
		plothorizontalaverage!(C14age, grd; mask, plot_options..., color, lab="$m", xlim=(0,NaN), xlab="Basin averaged radiocarbon age")
	end
	p
end

# ╔═╡ cc9fb021-35f9-4183-b992-40d82a5c1d2b
md"""
## Special UI stuff
"""

# ╔═╡ b505808f-6a05-4f2a-b755-201065b9c5f5
function plotclicktracker(p::Plots.Plot; draggable::Bool=false)

	# we need to render the plot before its dimensions are available:
	# plot_render = repr(MIME"image/svg+xml"(),  p)
	plot_render = repr(MIME"image/svg+xml"(),  p)

	# these are the _bounding boxes_ of our plot
	big = bbox(p.layout)
	small = plotarea(p[1])

	# the axis limits
	xl = xlims(p)
	yl = ylims(p)

	# with this information, we can form the linear transformation from 
	# screen coordinate -> plot coordinate

	# this is done on the JS side, to avoid one step in the Julia side
	# we send the linear coefficients:
	r = (
	x_offset = xl[1] - (xl[2] - xl[1]) * small.x0[1] / small.a[1],
	x_scale = (big.a[1] / small.a[1]) * (xl[2] - xl[1]),
	y_offset = (yl[2] - yl[1]) + (small.x0[2] / small.a[2]) * (yl[2] - yl[1]) + yl[1],
	y_scale = -(big.a[2]/ small.a[2]) * (yl[2] - yl[1]),
	x_min = xl[1], # TODO: add margin
	x_max = xl[2],
	y_min = yl[1],
	y_max = yl[2],
	)

	HTML("""<script id="hello">
		const body = $(PlutoRunner.publish_to_js(plot_render))
		const mime = "image/svg+xml"
		const img = this ?? document.createElement("img")
		let url = URL.createObjectURL(new Blob([body], { type: mime }))
		img.type = mime
		img.src = url
		img.draggable = false
		
		const clamp = (x,a,b) => Math.min(Math.max(x, a), b)
		img.transform = f => [
			clamp(f[0] * $(r.x_scale) + $(r.x_offset), $(r.x_min), $(r.x_max)),
			clamp(f[1] * $(r.y_scale) + $(r.y_offset), $(r.y_min), $(r.y_max)),
		]
		img.fired = false
		
		const val = {current: undefined }
		
		if(this == null) {
		Object.defineProperty(img, "value", {
			get: () => val.current,
			set: () => {},
		})
		
		const handle_mouse = (e) => {
			const svgrect = img.getBoundingClientRect()
			const f = [
				(e.clientX - svgrect.left) / svgrect.width, 
				(e.clientY - svgrect.top) / svgrect.height
			]
		console.log(f)
			if(img.fired === false){
				img.fired = true
				val.current = img.transform(f)
				img.dispatchEvent(new CustomEvent("input"), {})
			}
		}
		img.addEventListener("click", onclick)
		img.addEventListener("pointerdown", e => {
			if($(draggable)){
				img.addEventListener("pointermove", handle_mouse);
			}
			handle_mouse(e);
		});
		const mouseup = e => {
			console.log(e)
			img.removeEventListener("pointermove", handle_mouse);
		};
		document.addEventListener("pointerup", mouseup);
		document.addEventListener("pointerleave", mouseup);
		}
		return img
		</script>""")
end

# ╔═╡ 7a209d2b-06b8-4177-a728-14a6b7e364ed
begin
	default_usage_error = :(error("Example usage:\n\n@intially [1,2] @bind x f(x)\n"))
	
	macro initially(::Any)
		default_usage_error
	end
	
	macro initially(default, bind_expr::Expr)
		if bind_expr.head != :macrocall || bind_expr.args[1] != Symbol("@bind")
			return default_usage_error
		end
		
		# warn if the first argument is a @bind
		if default isa Expr && default.head == :macrocall && default.args[1] == Symbol("@bind")
			return default_usage_error
		end
			
		esc(intially_function(default, bind_expr))
	end
	
	
	function intially_function(default, bind_expr)
		sym = bind_expr.args[3]
		@gensym setval bond

		quote
			if !@isdefined($sym)
				$sym = $default
			end

			$setval = $sym


			$bond = @bind $sym $(bind_expr.args[4])
			PlutoRunner.Bond

			if $sym isa Missing
				$sym = $setval
			end

			$bond
		end
	end
end

# ╔═╡ bd6d5348-56e9-45d0-8770-0499d437b0a6
depth_slider = @bind depth Slider(0:50:6000, show_value=true, default=500)

# ╔═╡ d2a72eda-685e-4033-acc4-d384f2816e80
md"""
depth = $(depth_slider)
"""

# ╔═╡ be2300ea-50a6-42b3-bc3b-3e65e01b3203
plothorizontalslice(C14age, grd; depth, st=:heatmap)

# ╔═╡ 25379449-9118-4d18-8688-1bedccd454b9
myhorizontalslice(C14age, grd; depth, title="Radiocarbon age at $(depth)m")

# ╔═╡ 1d63cee3-0d5d-497d-b019-ce9c8865c4d4
md"""
depth = $(depth_slider)
"""

# ╔═╡ 2dbc5948-3c9e-4a93-9919-702fab2a992d
function prettylon(lon)
	lon = round(Int, lon)
	lon == 0 ? "0°" : lon > 0 ? "$(lon)°E" : "$(-lon)°W"
end

# ╔═╡ 1a9b44b4-415f-47aa-acae-dbef7c9af54a
function prettylat(lat)
	lat = round(Int, lat)
	lat == 0 ? "Eq" : lat > 0 ? "$(lat)°N" : "$(-lat)°S"
end

# ╔═╡ dd59eee8-2265-4487-a875-a09aca2d9cb0
function slice_and_profile(click_coordinate)
	p = plothorizontalslice(C14age, grd; depth, title="Radiocarbon age at $(depth)m", color=:viridis, clim, colorbar=false)
	plothorizontalslice!(p, C14age, grd; depth, seriestype=:contourf, levels=0:125:2500, clim, color=:viridis, lw=0)
	plothorizontalslice!(p, C14age, grd; depth, seriestype=:contour, levels=0:250:5000, color=:black, contour_labels=true)
	# profile
	if !isnothing(click_coordinate)
		lon, lat = click_coordinate
		lon2, lat2 = prettylon(lon), prettylat(lat)
		scatter!(p, [lon], [lat], c=:black, lab="")
		annotate!(lon, lat-10, text("$(lon2), $(lat2)", 10, :red))
		vline!(p, [lon], lab="", c=:red, linestyle=:dash)
		hline!(p, [lat], lab="", c=:red, linestyle=:dash)
		# profile
		pro = plotdepthprofile(C14age, grd; lonlat=(lon,lat), xlim=(0,3000), yunit=u"km", lab="", xlab="Age")
		hline!(pro, [depth * u"m"], lab="", c=:red, linestyle=:dash)
	else
		pro = plothorizontalaverage(C14age, grd, xlim=(0,3000), yunit=u"km", lab="", xlab="Age")
	end
	xlims!(p, (0,360))
	ylims!(p, (-90,90))
	# home-made colorbar
	crange = range(clim..., length=100)
	cbar = contourf(crange * u"yr", [0,1], [crange crange]'; color=:viridis, levels=0:125:2500, clim, colorbar=false, yticks=[], xlab="Radiocarbon age", lw=0)
	contour!(cbar, crange * u"yr", [0,1], [crange crange]'; color=:black, levels=0:250:5000) 
	# combine plots
	p1 = plot(p, cbar, layout=grid(2, 1, heights=(0.925, 0.075)))
	plot(p1, pro, layout=grid(1, 2, widths=(0.8, 0.2)), size=(800,400), bottom_margin=3Plots.mm)
end

# ╔═╡ d8ac806d-99b8-4a69-a31d-a92711b274af
mainplot = @initially [160.0,0.0] @bind x0 plotclicktracker(slice_and_profile(x0); draggable=true)

# ╔═╡ e552da86-a408-4b3f-a0b6-eb9e9702ddc5
lon, lat = x0 # a tuple of longitude and latitude sliders

# ╔═╡ ebce3a87-1fd8-4f63-a7ac-c97476c15323
let
	p = plotmeridionalslice(C14age, grd; lon, title="Radiocarbon age at $(prettylon(lon))", plot_options...)
	# Overlay black contours (nicer on the eyes)
	plotmeridionalslice!(C14age, grd; lon, st=:contour, levs=0:250:3000, clabels=true, c=:black, plot_options...)
	# Overlay red lines for latitude of target latitude and depth
	#hline!([depth], lab="", c=:red, linestyle=:dash)
	#vline!([lat], lab="", c=:red, linestyle=:dash)
end

# ╔═╡ 66364ec2-d852-4d02-940f-feafd7deb53d
let
	plotdepthprofile(C14age, grd; lonlat=(lon,lat), title="Radiocarbon age at ($(prettylon(lon)), $(prettylat(lat)))", xlim=(0,NaN), plot_options..., lab="")
end

# ╔═╡ 1022c705-7bb4-4f6f-a5b7-834d1c9aaaa6
md"""
## Notebook environment
"""

# ╔═╡ Cell order:
# ╟─c91da4d9-a3d1-4c38-aadb-c4548b6cc0e3
# ╟─41e825fa-43e2-4b33-9e25-98145cb4eb63
# ╠═d8ac806d-99b8-4a69-a31d-a92711b274af
# ╟─d2a72eda-685e-4033-acc4-d384f2816e80
# ╟─1d5df523-879e-4398-849e-d4b21003dd12
# ╟─fc00184b-9bcc-4833-9fc8-49ef2efbd3af
# ╟─2b531ea8-87eb-4a63-aabc-e7520f18a893
# ╟─59304232-cf0f-47ec-b993-86afcd2cc227
# ╟─5473f343-4dd5-4aef-88c9-4287a876d7f9
# ╟─ba54f82c-9e6c-409d-86bd-a85cc4609357
# ╟─38129859-df5e-4faa-a75c-8304b2bed55b
# ╠═0a889312-d965-11eb-28e0-a72a1379850a
# ╟─56bf3614-092d-4891-b1aa-7b2c5e1ddcbc
# ╟─4bbffe2c-0520-4e7b-bf38-0857454336ef
# ╠═35641e0f-537c-476a-86bd-d690d0989736
# ╟─bcc5fa0f-7595-4a21-a057-1cbabf94a049
# ╟─d76ec9fb-77f4-4a14-afb2-713cbd03ffab
# ╠═ba5fb3fa-b9f8-43c3-881c-0dc26c399b70
# ╟─5862c0fe-ff2c-45c9-aef4-ab0635442ee5
# ╠═ac71b40d-5f79-4881-9086-c71c89077499
# ╟─fe54d279-6d03-4892-9c0f-9dccb42b77a1
# ╠═daef9e51-08b5-4ab3-946a-703052ef1a34
# ╟─3c300825-ef08-4e5d-ab29-9264bb8e1024
# ╠═5b99fcb5-c3e1-4a01-84e6-334b26562301
# ╟─c9f07270-f1ce-4b1b-b8b0-a4ab53440366
# ╠═021a82d2-4a7b-4768-9e49-bdcb2a2b5da4
# ╟─7679b0ba-5a04-4af3-9f58-80644da08f9e
# ╠═202247cd-7a29-4f24-a5f9-831badce334f
# ╟─55292329-200f-4cf3-bb1d-e4cff6def15b
# ╠═e533c789-73e5-4843-b58e-9b00c950aa0b
# ╟─fa2265a6-3bf5-4b29-9fcb-271d1215ea62
# ╠═9a40e63d-bcb0-47d9-9daa-f98e6b2d90d7
# ╟─4e82b6b7-abb8-4248-ad92-65e082c2eedc
# ╠═8875b4e6-95ca-4c31-a6c5-70fd7193f5a5
# ╟─77d30cef-dc91-4db1-886f-15b4787f1ec7
# ╠═fa03c28c-5c2a-4764-9260-6acd9b5ec5a3
# ╟─da45685f-db9b-4c8a-a62b-0bc9cfc04311
# ╠═d3074335-0625-402a-b4d3-6ae65c8cd168
# ╟─5ef6f676-ead9-486d-afdd-b145d48e5a24
# ╠═59dd2d1e-78e2-431d-ad92-a7368b8e54ce
# ╟─fdc77b5c-5a56-43bc-b304-89f1f064c9a5
# ╟─5d8b2518-66a5-4ab2-992b-b80b98d1a498
# ╠═4d57848c-a2e7-4d5a-b848-f3822437da6a
# ╠═be2300ea-50a6-42b3-bc3b-3e65e01b3203
# ╠═25379449-9118-4d18-8688-1bedccd454b9
# ╠═d4521113-4a5c-4154-877f-7228885e0ac2
# ╟─1d63cee3-0d5d-497d-b019-ce9c8865c4d4
# ╠═ebce3a87-1fd8-4f63-a7ac-c97476c15323
# ╠═66364ec2-d852-4d02-940f-feafd7deb53d
# ╟─0268bff0-7ccd-404f-8424-78444885fb4a
# ╠═f0d1a772-2b1d-4feb-8d52-b9dbf29b3f45
# ╠═75b087ff-fb90-4a79-b95b-678157f4046e
# ╟─9a95c55a-8e44-42d3-a52e-2f18d670811b
# ╠═ee1d8881-4dce-45d3-95fa-9bb33819ebe9
# ╟─a6407c69-ac39-4819-9429-239f43a5dc63
# ╟─56a1aa18-0134-4e4f-a6cc-3662316284dc
# ╟─3db99f8c-f238-427c-8239-19a7d94401b7
# ╠═361bbe46-bf9a-4a6c-834f-5a0660f9a0e6
# ╠═e0267345-df37-410b-bcca-218e35388f87
# ╠═246063c6-d2ca-4796-88c5-fa1489e6381b
# ╠═e14f6526-c5e3-4c17-896e-2983c030f5e9
# ╠═98508321-fb4b-4076-9c95-0d38f63303b8
# ╠═8f76c559-ef89-46c0-a792-75a1de96f6e5
# ╟─f7ad9931-1d92-4d08-a0ee-9393dd524f92
# ╟─25f7c2a3-a8a7-441e-8da1-b311fd355867
# ╠═c08308c9-eade-4209-b755-2415542e542a
# ╟─cffb9e34-c2ab-40c2-8006-e1fd934fd066
# ╠═d5a447b4-2d22-4466-b89f-de859ed9f09e
# ╠═9f2ece88-dbc7-495e-9ca5-6518af58caf8
# ╟─cc9fb021-35f9-4183-b992-40d82a5c1d2b
# ╠═325084cb-b233-4a63-8f8c-7ba32b7e8d43
# ╠═b505808f-6a05-4f2a-b755-201065b9c5f5
# ╠═7a209d2b-06b8-4177-a728-14a6b7e364ed
# ╠═dd59eee8-2265-4487-a875-a09aca2d9cb0
# ╟─e552da86-a408-4b3f-a0b6-eb9e9702ddc5
# ╠═bd6d5348-56e9-45d0-8770-0499d437b0a6
# ╠═2dbc5948-3c9e-4a93-9919-702fab2a992d
# ╠═1a9b44b4-415f-47aa-acae-dbef7c9af54a
# ╟─1022c705-7bb4-4f6f-a5b7-834d1c9aaaa6
# ╠═0df28961-4e34-496b-86c0-a4740b368127
