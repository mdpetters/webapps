### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 37ff1a70-02b5-11eb-1507-41208db7f97c
begin
	import Pkg
        Pkg.activate(Base.current_project())
	
	using DifferentialMobilityAnalyzers
	using Gadfly
	using LinearAlgebra
	using Printf
	using Distributions
	using DataFrames
	using Random
	using LsqFit
	using MLStyle
	using Optim
	using CSV
	using RegularizationTools
	using Lazy
	using Colors
	using Underscores
	using Dates
	using PlutoUI
	using DataStructures
	using Interpolations
	using Observables
	using NumericIO
	using Cairo
	using Fontconfig
	using LinearAlgebra
	using ProgressLogging
	
	import PlutoUI: combine
        Logging.disable_logging(Logging.Warn)
	
end

# ╔═╡ fb64cac1-928f-4555-8cdc-fc11c8adc1a1
module DataLoader

using DifferentialMobilityAnalyzers
using DataFrames
using Dates

struct ParticleDistributions
    timestamp::Array{DateTime,1}
    𝕟::Array{SizeDistribution,1}
end

struct HTDMADistributions
    timestamp::Array{DateTime,1}
	Dd::Array{Float64,1}
    𝕟::Array{SizeDistribution,1}
end

# Accepts a DataFrame, which is materialized from a CSV file upon load
function df_to_psd(df, i)
    Dl = Vector(df[1, i:end])
    Dp = Vector(df[2, i:end])
    Du = Vector(df[3, i:end])
	ts =  Millisecond.(Vector(df[:, 3])) .+ DateTime(0,12,31,0,0,0)

	De = [Dl[1:end];Du[end]]
	ΔlnD = log.(De[2:end]./De[1:end-1])
	ℕ = map(4:length(df[1:end,1])) do j
	N = Vector(df[j, i:end])
		
		𝕟 = SizeDistribution( # note that order must be from high D to low D
			[[]],             # Empty array - not used
			reverse(De),      # Bin edges
			reverse(Dp),      # Bin midpoints
			reverse(ΔlnD),    # Log bin width
			reverse(N./ΔlnD), # Spectral density
			reverse(N),       # Number concentration per bin
			:data             # Label for distribution origin, not used
		)
	end
	
	return ParticleDistributions(ts, ℕ)
end

function df_to_htdma(df, i)
	ℕ = map(1:4:length(df[1:end,1])) do j
	    Dl = Vector(df[j, i:end])
	    Dp = Vector(df[j+1, i:end])
	    Du = Vector(df[j+2, i:end])
	
		De = [Dl[1:end];Du[end]]
		ΔlnD = log.(De[2:end]./De[1:end-1])

		N = Vector(df[j+3, i:end])
		𝕟 = SizeDistribution(
			[[]], 
			reverse(De), 
			reverse(Dp), 
			reverse(ΔlnD), 
			reverse(N./ΔlnD), 
			reverse(N), 
			:data
		)
	end

	ts = Millisecond.(Vector(df[4:4:end, 3])) .+ DateTime(0,12,31,0,0,0) 
	Dd = Vector(df[4:4:end, 4]) .* 1e-9
	
	return HTDMADistributions(ts, Dd, ℕ)
end

end

# ╔═╡ 000bfaa0-3e61-41d3-a460-448eac65f381
module TDMA

using DifferentialMobilityAnalyzers
using Memoize
using Underscores
using Lazy
using RegularizationTools
using LinearAlgebra
using Interpolations
using MLStyle
using LsqFit

function getModel(𝕟ᵢₙ, 𝕣ᵢₙ, Λ₁ᵢₙ, Λ₂ᵢₙ)
    Λ₁, Λ₂, 𝕟, 𝕣 = deepcopy(Λ₁ᵢₙ), deepcopy(Λ₂ᵢₙ), deepcopy(𝕟ᵢₙ), deepcopy(𝕣ᵢₙ)
	δ₁ = setupDMAgridded(Λ₁, (𝕟.De))
	δ₂ = setupDMAgridded(Λ₂, (𝕣.De))

    O(k) = mapfoldl(
		zs -> (δ₂.Ω(Λ₂, δ₁.Z, zs/k, k) .* δ₂.Tl(Λ₂, δ₁.Z, k))', vcat, δ₂.Z)
	T₁(zˢ, k) = δ₁.Ω(Λ₁, δ₁.Z, zˢ / k, k) .* δ₁.Tc(k, δ₁.Dp) .* δ₁.Tl(Λ₁, δ₁.Dp)
	DMA₁(𝕟, zˢ, gf) = @_ map((gf ⋅ (T₁(zˢ, _) * 𝕟)), 1:3)
	DMA₂(𝕟, k) = O(k) * 𝕟
	
	function TDMA(𝕟, zˢ, gf)
		ℕ = @> DMA₁(𝕟, zˢ, gf) 
		r = mapreduce(k -> DMA₂(ℕ[k].N, k), +,  1:3) 
		return SizeDistribution([[]], δ₂.De, δ₂.Dp, δ₂.ΔlnD, r./δ₂.ΔlnD, r, :TDNA)
	end

	model(𝕟, P, Dd, gf) = 
		sum(@_ map(P[_] * TDMA(𝕟, dtoz(Λ₁, Dd), gf[_]), 1:length(P)))

	return model, δ₁, δ₂
end

function getMatrix(𝕟, δ₂, Dd, model)
	mgf = δ₂.Dp ./ (Dd * 1e9)
	λ = @_ map(model(𝕟, 1.0 , Dd, _), mgf) 
	𝐁 = hcat((@_ map(_.N, λ))...)
	return 𝐁
end

RMSE(x,y) = sqrt(sum((x .- y).^2.0)./length(x))

function invertme(p) 
	Ax = Float64[]

	xλ = @match p.method begin
		"L₀DₓB" => begin
			try
				invert(p.𝐁, p.𝕣.N, LₖDₓB(0, 0.001, p.lb, p.ub))
			catch
				zeros(length(p.𝕣.N))
			end
		end
		"LSQ₁"  => begin
			init = map(x->parse(Float64,x), split(p.initial, ","))
			f(_, x) = (p.model(p.𝕟, x[1], p.Dd, x[2])).N
			fit = curve_fit(f, p.𝕣.N, p.𝕣.N, init, lower = zeros(2))
			Ax = fit.param
			kk = argmin(abs.(Ax[2] .- p.gf))
			@> zeros(length(p.gf)) setindex!(Ax[1], kk)
		end
		"LSQ₂"  => begin
			init = map(x->parse(Float64,x), split(p.initial, ","))
			f(_, x) = (p.model(p.𝕟, x[1:2], p.Dd, x[3:4])).N
			fit = curve_fit(f, p.𝕣.N, p.𝕣.N, init, lower = [0.0,0.0,1.0,1.0])
			Ax = fit.param
			kk = argmin(abs.(Ax[3] .- p.gf))
			ll = argmin(abs.(Ax[4] .- p.gf))
			@> zeros(length(p.gf)) setindex!(Ax[1], kk) setindex!(Ax[2], ll)
		end
		_ => throw("Not supported")
	end

	pred = @match p.method begin
		"L₀DₓB" => p.𝐁 * xλ
		"LSQ₁"  => (p.model(p.𝕟, Ax[1], p.Dd, Ax[2])).N 
		"LSQ₂"  => (p.model(p.𝕟, Ax[1:2], p.Dd, Ax[3:4])).N 
	end
	
	r = round(RMSE(p.𝕣.N, pred), digits = 2)
	
	return xλ, Ax, pred, r
end

end

# ╔═╡ 566cb6a0-0332-11eb-36e1-3bd2acf329ab
md"""
$(TableOfContents(depth=5))
# Introduction 

Welcome to this interactive tandem DMA inversion tool. It is to illustrate the inversion method described in [Petters 2021](https://amt.copernicus.org/articles/14/7909/2021/amt-14-7909-2021.pdf).

After each change, the reactive Pluto notebook revaluates all of the cells. This behavior is equivalent to that of a spreadsheet application.  
	
$(Markdown.MD(Markdown.Admonition("note", "Author", [md" 
If you have questions or comments, please send an email to\
	Markus Petters: **markus.petters@ucr.edu**"])))

## Acknowledgements

$(Resource("https://raw.githubusercontent.com/mdpetters/webapps/main/virtualDMA/notebooks/ASR.png", :width => 300)) 

$(Resource("https://raw.githubusercontent.com/mdpetters/webapps/main/virtualDMA/notebooks/nsflogo.jpg", :width => 300)) 

The tandem DMA inversion code was supported by DOE grant SC 0021074. Development of this online tool was supported by NSF grant AGS-2037704. Infrastructure for hosting this notebook online is supported by NSF-AGS-2112978.

## References

Petters, M.D. (2018) A language to simplify computation of differential mobility analyzer response functions Aerosol Science & Technology, 52 (12), 1437-1451, https://doi.org/10.1080/02786826.2018.1530724.

Petters, M. D.: Revisiting matrix-based inversion of scanning mobility particle sizer (SMPS) and humidified tandem differential mobility analyzer (HTDMA) data, Atmos. Meas. Tech., 14, 7909–7928, https://doi.org/10.5194/amt-14-7909-2021, 2021.
"""



# ╔═╡ a817257a-e8aa-4e00-890f-dea76d03de76
md" # Module DataLoader
This module loads the data and places it in the appropriate data structre. Modify this module if you like to work with a different file format or alternatively, format the data as in the provided sample files. 
"

# ╔═╡ 1c253882-0332-11eb-2f0b-43c1b874a6ab
md" # Load Files
"

# ╔═╡ 5de4e17d-af63-476a-9a32-e928ca3e229a
begin
	@bind values confirm( 
		combine() do Child
	
	md"""**Enter Filenames**

	Enter the filename of the SMPS and TDMA Data. See [example csv files](https://drive.google.com/drive/folders/1vV--usvcJgH54bMsFKjmYIPa37kAFJ3J?usp=sharing) for the required file format. The start index is corresponds to the column with the first size bin. Diameters are in [nm]. Columns 5 to n can be used to store metadata about the scan, e.g. RH, temperature etc. Note that timeISO is provided for compatibility/verification reasons. Time int64 is used to avoid parsing ambiguously formated inputs. Press submit to when you change the index. Files load instantaneously.  
	
	| Instrument |  Upload File | Start Index | 
	| ---------- | --------- | ----------  |
	| SMPS       | $(@bind file1 FilePicker()) |  $(Child(NumberField(1:20, default = 4))) |
	| TDMA       | $(@bind file2 FilePicker()) |  $(Child(NumberField(1:20, default = 5))) |
	"""
			
	end
	)
end

# ╔═╡ 15f2f933-d961-4796-8607-a8967da3d023
	begin
		sps = try 
			file1["data"]
		catch
			nothing
		end
		isnothing(sps) || (open(io -> write(io, sps), "smps.csv", "w"))
		
		htd = try 
			file2["data"]
		catch
			nothing
		end
		isnothing(htd) || (open(io -> write(io, htd), "htdma.csv", "w"))

	 	smps = @> CSV.read("smps.csv", DataFrame) DataLoader.df_to_psd(values[1])
	 	htdma = @> CSV.read("htdma.csv", DataFrame) DataLoader.df_to_htdma(values[2])
		
		Markdown.MD(
			Markdown.Admonition("info", "Information", [md"Data are loaded here and then preprocessed using the DataLoader module. Potential error messages during file loading will display here."]))
	end

# ╔═╡ 84c14876-f2fd-4bef-9829-75b2dc702861
md" # Module TDMA 

"

# ╔═╡ ed4f6300-e18d-4821-962b-4e1dc9a535a2
Markdown.MD( Markdown.Admonition("info", "Information", [md" The code applies the concepts from [DifferentialMobilityAnalyzers.jl](https://github.com/mdpetters/DifferentialMobilityAnalyzers.jl). The TDMA method is described in [Petters 2021](https://amt.copernicus.org/articles/14/7909/2021/amt-14-7909-2021.pdf). Modify this module to change the model (e.g. swap charge functions). You can set properties such as flow rates using the interactive controls below."]))

# ╔═╡ 2984464e-0332-11eb-303d-cf9502464cbc
md""" ## Functions
Each cell contains a function that links the control elements/inputs and generates the outputs. The fields in the "Input Conditions" section are queried inside these functions. The main function is invertTDMA. Modify this function if you like to alter the behavior of the app.
"""

# ╔═╡ c7d9e924-02b6-11eb-2dfc-998639d0751f
function plotGFDual(dfl, dfr)
    set_default_plot_size(18cm, 7cm)
    colors = ["black", "darkred", "steelblue3", "darkgoldenrod"]
    xlabels = collect(1:0.5:3)
    p1 = plot(dfl, x = :gf, y = :N, color = :Color, Geom.step,
		Guide.xlabel("Growth Factor (-)"),
		Guide.ylabel("Number (cm⁻³)", orientation = :vertical),
		Guide.xticks(ticks = collect(0.8:0.1:2.5)),
		Scale.color_discrete_manual(colors...),
		Scale.x_continuous(labels = x -> x in xlabels ? @sprintf("%.1f", x) : ""),
		Scale.y_continuous(labels = x -> formatted(x, :SI, ndigits=2)),
		Theme(plot_padding=[1mm,5mm,1mm,1mm]),
		Coord.cartesian(xmin = 0.8, xmax = 2.5))

    p2 = plot(dfr, x = :gf, y = :Frequency, color = :Color, Geom.step,
		Guide.xlabel("Growth Factor (-)"),
		Guide.ylabel("Probability Density (-)", orientation = :vertical),
		Guide.xticks(ticks = collect(0.8:0.1:2.5)),
		Scale.color_discrete_manual(colors...),
		Scale.x_continuous(labels = x -> x in xlabels ? @sprintf("%.1f", x) : ""),
		Theme(plot_padding=[1mm,6mm,1mm,1mm]),
		Coord.cartesian(xmin = 0.8, xmax = 2.5))
	
    hstack(p1,p2)
end

# ╔═╡ fe9e6990-030f-11eb-3896-ed3b8e0a9e16
function plotPSD(dfl, dfr)
	colors = ["black", "darkred"]
    set_default_plot_size(18cm, 7cm)
	
	xt = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500]
	xlabels = log10.([10, 20, 50, 100, 200, 500])
    
	p1 = plot(dfl, x = :Dp, y = :S, color = :Dist, Geom.line,
		Guide.xlabel("Particle diameter (nm)"),
		Guide.ylabel("Number (cm⁻³)", orientation = :vertical),
		Scale.color_discrete_manual(colors...),
        Guide.colorkey(; title = "GF Raw"),
        Guide.xticks(ticks = log10.(xt)),
		Scale.y_continuous(labels = x -> formatted(x, :SI, ndigits=2)),
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
		Theme(plot_padding=[1mm,3mm,1mm,1mm]),
        Coord.cartesian(xmin = log10(10), xmax = log10(500)))
	
    p2 = plot(dfr, x = :Dp, y = :S, color = :Dist, Geom.step,
        Guide.xlabel("Particle diameter (nm)"),
        Guide.ylabel("dN/dlnD (cm⁻³)"),
        Guide.xticks(ticks = log10.(xt)),
        Guide.colorkey(; title = "PSD"),
		Scale.y_continuous(labels = x -> formatted(x, :SI, ndigits=2)),
        Scale.x_log10(labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""),
		Theme(plot_padding=[1mm,5mm,1mm,1mm]),
        Scale.color_discrete_manual(colors...),
        Coord.cartesian(xmin = log10(10), xmax = log10(500)))
	
	hstack(p2, p1)
end

# ╔═╡ 90a47615-784d-427c-ac94-735855369af0
function setupDMA(t, p, qsh, qsa, leff, DMAtype)
	(r₁, r₂, l, geom) = @match DMAtype begin
		"TSI Long"      => (9.37e-3, 1.961e-2, 0.44369, :cylindrical)
		"High Flow DMA" => (0.05, 0.058, 0.6, :cylindrical)
		"RDMA"          => (2.4e-3, 50.4e-3, 10e-3, :radial)
		"Brechtel"      => (6.24e-2, 7.23e-2, 0.339, :cylindrical)
		_               => throw("Error")
	end
	
	return DMAconfig(t, p, qsa, qsh, r₁, r₂, l, leff, :-, 6, geom)
end

# ╔═╡ 49d26bd2-3a3d-4147-aad7-07402a5ca836
function assemble(result, ge)
	gf = sqrt.(ge[2:end] .* ge[1:end-1])
	
	dummy1 = DataFrame(
        Label = "Lower bin bounds",
		timeISO8601 = result.timeISO8601,
		Dd = result.Dd,
		Method = result.Method,
		RMSE = result.RMSE,
		Integral = result.Integral,
		Fit = result.Fit,
	)
    map(i -> dummy1[!, Symbol("g$i")] = [ge[i+1]], 1:length(gf))

    dummy2 = DataFrame(
        Label = "Midpoints",
		timeISO8601 = result.timeISO8601,
		Dd = result.Dd,
		Method = result.Method,
		RMSE = result.RMSE,
		Integral = result.Integral,
		Fit = result.Fit,
	)
    map(i -> dummy2[!, Symbol("g$i")] = [gf[i]], 1:length(gf))

    dummy3 = DataFrame(
        Label = "Upper bin bounds",
		timeISO8601 = result.timeISO8601,
		Dd = result.Dd,
		Method = result.Method,
		RMSE = result.RMSE,
		Integral = result.Integral,
		Fit = result.Fit,
	)
    map(i -> dummy3[!, Symbol("g$i")] = [ge[i]], 1:length(gf))

 	[dummy1; dummy2; dummy3; result]
end

# ╔═╡ 284aaeac-0a52-11eb-1eaf-d59f0ac49aff
begin
	aaa = @bind tc NumberField(15.0:1.0:35.0, default = 22.0)
	bbb = @bind phPa NumberField(800:10.0:1020, default = 1010.0)
	
	ccc = @bind qshl1 NumberField(2.0:0.1:15.0, default = 5.0)
	ddd = @bind qsal1 NumberField(0.3:0.1:2, default = 0.63)
	
	eee = @bind DMAtype1 Select(
		["TSI Long", "High Flow DMA", "Radial DMA", "Brechtel"],
		default = "Brechtel"
	)
	fff = @bind leff1 NumberField(0.0:1:100, default = 0.0)

	hhh = @bind qshl NumberField(2.0:0.1:15.0, default = 5.0)
	iii = @bind qsal NumberField(0.3:0.1:2, default = 1.0)
	jjj = @bind DMAtype Select(
		["TSI Long", "High Flow DMA", "Radial DMA", "Brechtel"],
		default = "Brechtel"
	)
	kkk = @bind leff NumberField(0.0:1:100, default = 0.0)
		
	dd = @bind initial1 TextField((15); default="1.0, 1.3")
	ee = @bind initial2 TextField((30); default="0.5, 0.5, 1.1, 1.6")

	cc = @bind method Select(
		["Pick Best", "L₀DₓB", "LSQ₁", "LSQ₂"]
	)

	ff = @bind hardbound CheckBox(default=true)
	ss = @bind method1 Select(
		[1, 5, length(htdma.timestamp)]
	)
	
	md"""

	## Input Conditions
	
	**Configure DMAs and Size Grid**
	
	Enter the appropriate values for your data. If you run the app locally and wish to specialize the app for a given instrument you can change the defaults by editing this code block. If you need a DMA with different geometry, you can add the label and edit :Setup to register the dimensions.

	Temperature (°C): $(aaa)
	Pressure (hPa): $(bbb)

	**DMA 1**
	
	DMA Sheath Flow (L min⁻¹): $(ccc)
	DMA Sample Flow (L min⁻¹): $(ddd)
	
	DMA Type: $(eee)
	Effective length: $(fff)

	**DMA 2**
	
	DMA Sheath Flow (L min⁻¹): $(hhh)
	DMA Sample Flow (L min⁻¹): $(iii)
	
	DMA Type: $(jjj)
	Effective length: $(kkk)
		
	**Invert Data**	
	Prohibit gf < 1: $(ff)

	LSQ₁ Initial Guess (f, gf): $(dd) 
	
	LSQ₂ Initial Guess (f1, f2, gf1, gf2): $(ee)
	
	Use Method: $(cc)
	Number of spectra to invert: $(ss)
	"""
end

# ╔═╡ 1c0ce412-f589-4e10-8086-a9a7d4462b28
begin
	t = tc + 273.15
	p = phPa*100.0
	init1 = parse.(Float64,split(initial1, ","))
	init2 = parse.(Float64,split(initial2, ","))
	Λ₁ = setupDMA(t, p, qshl1*1.6666666666e-5, qsal1*1.6666666666e-5, leff1, DMAtype1)
	Λ₂ = setupDMA(t, p, qshl*1.6666666666e-5, qsal*1.6666666666e-5, leff, DMAtype)
	:Setup
end

# ╔═╡ 60f16b68-ef53-476c-b140-5666c35ad541
function test()
	i = 1
	𝕟 = smps.𝕟[i]
	𝕣 = htdma.𝕟[i]
	Dd = htdma.Dd[i]
	model, δ₁a, δ₂a = TDMA.getModel(𝕟, 𝕣, Λ₁, Λ₂)
	𝐁 = TDMA.getMatrix(𝕟, δ₂a, Dd, model)
	mgf = δ₂a.Dp ./ (Dd * 1e9)
	
	Normalize(x) = x./sum(x)
	f = @> (pdf.(Normal(1.3,0.1), mgf)) Normalize
	o = model(𝕟, f, Dd, mgf)
	plot(layer(x = mgf, y = 𝐁*f), layer(x = δ₂a.Dp./(Dd * 1e9), y = o.N, Geom.line))
end

# ╔═╡ f1a7d609-3943-41e2-9aae-8e08016c5cae
function invertTDMA(i, method)
	𝕟 = smps.𝕟[i]
	𝕣 = htdma.𝕟[i]
	Dd = htdma.Dd[i]

	model, δ₁, δ₂ = TDMA.getModel(𝕟, 𝕣, Λ₁, Λ₂)
	ge = δ₂.De ./ (Dd * 1e9) 
	gf = δ₂.Dp ./ (Dd * 1e9) 
	Δgf =  ge[1:end-1] .- ge[2:end]

	Normalize(x) = x./sum(x)

	function sanitize(x)
		y = deepcopy(x)
		y[y .< 1e-4] .= 0
		return y
	end
	function LDB()
		𝐁 = TDMA.getMatrix(𝕟, δ₂, Dd, model)
		lb = zeros(length(gf))
		ub = ones(length(gf))
		if hardbound == true
			ii = gf .< 1
			ub[ii] .= 0.0
		end	
		p1 = (method = "L₀DₓB", 𝐁 = 𝐁, 𝕣 = 𝕣, lb = lb, ub = ub)
		xλ, Ax, pred, r = TDMA.invertme(p1)
	end

	function LSQ1()
		p2 = (
			method = "LSQ₁", 
			model = model, 
			gf = gf, 
			𝕣 = 𝕣, 
			Dd = Dd, 
			initial = initial1, 
			𝕟 = 𝕟
		)
		xλ, Ax, pred, r = TDMA.invertme(p2)
	end

	function LSQ2()
		p3 = (
			method = "LSQ₂", 
			model = model, 
			gf = gf, 
			𝕣 = 𝕣, 
			Dd = Dd, 
			initial = initial2, 
			𝕟 = 𝕟
		)
		xλ, Ax, pred, r = TDMA.invertme(p3)
	end

	zz = 1
	xλ, Ax, pred, r = @match method begin
		"L₀DₓB"  =>  LDB()
		"LSQ₁"   =>  LSQ1()
		"LSQ₂"   =>  LSQ2()
		"Pick Best" => begin
			xλ1, Ax1, pred1, r1 = LDB()
			xλ2, Ax2, pred2, r2 = LSQ1()
			xλ3, Ax3, pred3, r3 = LSQ2()
			zz = argmin([r1,r2,r3])
			@match zz begin
				1 => (xλ1, Ax1, pred1, r1)
				2 => (xλ2, Ax2, pred2, r2)
				3 => (xλ3, Ax3, pred3, r3)
			end
		end
	end

	mydict = Dict([(1, "L₀DₓB"), (2, "LSQ₁"), (3,"LSQ₂")])
	(length(Ax) == 0) && (Axstr = "")
	(length(Ax) == 2) && (Axstr = @sprintf("g = %.2f", Ax[2]))
	(length(Ax) == 4) && (Axstr = 
		@sprintf("f = %.2f %.2f\ng = %.2f %.2f", Ax[1], Ax[2], Ax[3], Ax[4]))

	Axstr1 = ""
	(length(Ax) == 2) && (Axstr1 = @sprintf("%2f %.2f", Ax[1], Ax[2]))
	(length(Ax) == 4) && (Axstr1 = 
		@sprintf("%.2f %.2f %.2f %.2f", Ax[1], Ax[2], Ax[3], Ax[4]))

	meth = (method == "Pick Best") ? mydict[zz] : method
	rf = @match meth begin
		"L₀DₓB" => round(sum(xλ), digits = 2) 
		"LSQ₁" => round(Ax[1], digits = 2) 
		"LSQ₂" => round(Ax[1] + Ax[2], digits = 2)
	end
	
	result = DataFrame(
		Label = "invertTDMA",
		timeISO8601 = htdma.timestamp[i],
		Dd = htdma.Dd[i],
		Method = meth,
		RMSE = r,
		Integral = rf,
		Fit = Axstr1,
	)
	theN = @_ map(_ < 0 ? 0 : _, Normalize(xλ./Δgf))
	map(i -> result[!, Symbol("g$i")] = [theN[i]], 1:length(xλ))

	tscn = Dates.Time(smps.timestamp[i])
 	tsht = Dates.Time(htdma.timestamp[i])

 	cn_label = ["$(tscn)" for i = 1:length(𝕟.Dp)]
 	ht_label = ["$(tsht)" for i = 1:length(𝕣.Dp)]
 	Ddnm = @sprintf("%i", Dd*1e9)
 	Dd_label = ["Dd = $(Ddnm) nm" for i = 1:2]
 	dfp1 = DataFrame(Dp = 𝕟.Dp, S = 𝕟.S, Dist = cn_label)
 	dfp2a = DataFrame(Dp = 𝕣.Dp, S = 𝕣.N, Dist = ht_label)
 	dfp2b = DataFrame(Dp = [Dd, Dd]*1e9, S = [0, maximum(𝕣.N)], Dist = Dd_label)
 	plot1 = plotPSD([dfp2a;dfp2b], dfp1)

	freq = (xλ./Δgf)
	freq[isnan.(freq)] .= 0
 	df1 = DataFrame(gf = 𝕣.Dp./(Dd.*1e9), N = 𝕣.N, Color = "raw")
 	df2 = DataFrame(gf = gf, N = pred, Color = "fitted")
	dfr = DataFrame(gf = gf, Frequency = freq, 
		Color = "$(meth)($r) \n $(Axstr)")
	
 	plot2 = plotGFDual([df1;df2], dfr)
	set_default_plot_size(16cm, 12cm)
	return vstack(plot1,plot2), result, ge

end

# ╔═╡ 4652778c-0a41-11eb-08c9-c53d0b312295
begin
	fig, res, ge = [], [], []
	try
		dump1 = read(`sh cull.sh`, String)
	catch
	end

	@progress for i = 1:method1 
		a,b,c = invertTDMA(i, method)
		push!(fig, a)
		push!(res, b)
		push!(ge, c)
	end
	map(i -> draw(PNG("spectrum$i.png", dpi = 300), fig[i]), 1:method1)
	dataDF = invertedData = map(assemble, res, ge) |> x -> vcat(x...)
	dataDF |> CSV.write("report.csv")
	try
		dump2 = read(`sh zip.sh`, String)
	catch
	end
	fig
end

# ╔═╡ a74e4657-9d9d-4ed9-8905-d92fd559ca57
dataDF

# ╔═╡ 2bb1a8c2-19b7-4ccc-a709-324bb1be800d
begin
	dataDF
	filename1 = "report.csv"
	xx = DownloadButton(read(filename1), basename(filename1))

	# filename2 = "figures.zip"
	# yy = DownloadButton(read(filename2), basename(filename2))

	md""" # Download Output

	The report.csv file contains a comma separated file with the data frame shown above. 
	
	| Data |  File |  
	| ---------------- | --------- | 
	| Inverted Data    | $(xx)  |    
	"""		
end

# ╔═╡ Cell order:
# ╟─566cb6a0-0332-11eb-36e1-3bd2acf329ab
# ╟─37ff1a70-02b5-11eb-1507-41208db7f97c
# ╟─a817257a-e8aa-4e00-890f-dea76d03de76
# ╟─fb64cac1-928f-4555-8cdc-fc11c8adc1a1
# ╟─1c253882-0332-11eb-2f0b-43c1b874a6ab
# ╟─5de4e17d-af63-476a-9a32-e928ca3e229a
# ╟─15f2f933-d961-4796-8607-a8967da3d023
# ╟─84c14876-f2fd-4bef-9829-75b2dc702861
# ╟─ed4f6300-e18d-4821-962b-4e1dc9a535a2
# ╟─000bfaa0-3e61-41d3-a460-448eac65f381
# ╟─60f16b68-ef53-476c-b140-5666c35ad541
# ╟─2984464e-0332-11eb-303d-cf9502464cbc
# ╟─c7d9e924-02b6-11eb-2dfc-998639d0751f
# ╟─fe9e6990-030f-11eb-3896-ed3b8e0a9e16
# ╟─90a47615-784d-427c-ac94-735855369af0
# ╟─1c0ce412-f589-4e10-8086-a9a7d4462b28
# ╟─f1a7d609-3943-41e2-9aae-8e08016c5cae
# ╟─49d26bd2-3a3d-4147-aad7-07402a5ca836
# ╟─284aaeac-0a52-11eb-1eaf-d59f0ac49aff
# ╟─4652778c-0a41-11eb-08c9-c53d0b312295
# ╟─a74e4657-9d9d-4ed9-8905-d92fd559ca57
# ╟─2bb1a8c2-19b7-4ccc-a709-324bb1be800d
