### A Pluto.jl notebook ###
# v0.12.18

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

# ╔═╡ bc2fdc0c-04b2-11eb-26dd-ab2d62794853
# Include package. 
begin 
	using Pkg
	Pkg.activate(joinpath(Pkg.envdir(), "dev-env"))
	using FractalTools 
	using Plots
	plotly()
	default(:label, "")
	default(:lw, 2)
end 

# ╔═╡ c5132b30-04b2-11eb-173c-898862565f53
# Define function
begin 
	f(x) = 1 / 2 - (x - 1 / 2)^2 + 1 
	g(x) = x^2
end

# ╔═╡ 924bbe12-04b3-11eb-3ab4-f1d84faeb86f
begin 
	xi, dx, xf = 0, 0.25,  1
	x = collect(xi : dx : xf) 
	y = f.(x) 
end

# ╔═╡ ee35dc64-04b3-11eb-22c3-5fc92810c3f1

@bind d html"<input type=range min=0 max=0.2 step=0.01>"

# ╔═╡ 30386406-04b4-11eb-2e1b-2d0eee87f5ed
@bind h  html"<input type=range min=0 max=0.2 step=0.01>"

# ╔═╡ 321a63d2-04b4-11eb-283d-9fa41a2534ef
@bind l  html"<input type=range min=0 max=0.2 step=0.01>"

# ╔═╡ 47922bdc-04b4-11eb-239b-411078a8b27b
@bind m  html"<input type=range min=0 max=0.2 step=0.01>"

# ╔═╡ 5f9872c4-04b4-11eb-162f-9f998adb0c9f
@bind niter html"<input type=range min=1 max=30 step=1>"

# ╔═╡ 508bc310-04b4-11eb-3d4c-e78390b8b8a1
# Construct interpolants 
begin 
	z = g.(x)
	dv =  d * ones(length(x) - 1)
	hv =  h * ones(length(x) - 1)
	lv =  l * ones(length(x) - 1)
	mv =  m * ones(length(x) - 1)
	interp = interpolate(x, y, d = dv, niter=niter)
	hiddeninterp = hiddenfif(x, y, z, d = dv, h = hv, l = lv, m = mv, niter=niter)
end

# ╔═╡ 34596672-04b5-11eb-1f41-f781db5f981f
d, h, l, m, niter

# ╔═╡ 99eb66ee-04b3-11eb-34e3-5f2bcf5f383f
begin 
	xt = collect(xi : 0.1dx : xf)
	plt = plot()
	scatter!(x, y, marker=(:circle, 3), label="data")
	plot!(xt, f.(xt), label="f")
	plot!(xt, interp.(xt), label="finterp")
	plot!(xt, getindex.(hiddeninterp.(xt), 1), label="finterphidden")
end

# ╔═╡ Cell order:
# ╠═bc2fdc0c-04b2-11eb-26dd-ab2d62794853
# ╠═c5132b30-04b2-11eb-173c-898862565f53
# ╠═924bbe12-04b3-11eb-3ab4-f1d84faeb86f
# ╠═508bc310-04b4-11eb-3d4c-e78390b8b8a1
# ╠═ee35dc64-04b3-11eb-22c3-5fc92810c3f1
# ╠═30386406-04b4-11eb-2e1b-2d0eee87f5ed
# ╠═321a63d2-04b4-11eb-283d-9fa41a2534ef
# ╠═47922bdc-04b4-11eb-239b-411078a8b27b
# ╠═5f9872c4-04b4-11eb-162f-9f998adb0c9f
# ╠═34596672-04b5-11eb-1f41-f781db5f981f
# ╠═99eb66ee-04b3-11eb-34e3-5f2bcf5f383f
