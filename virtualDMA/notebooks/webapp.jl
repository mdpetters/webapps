### A Pluto.jl notebook ###
# v0.17.7

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

# ╔═╡ dc1a7d9c-fad9-11eb-3360-4be773edd179
begin    
	import Pkg
    Pkg.activate(Base.current_project())

	using Gadfly
	using PlutoUI
	using Colors
	using DifferentialMobilityAnalyzers
	using Markdown
	using MKL
	using Lazy
	using MLStyle
	using Interpolations
	using Lazy
	using Underscores
	using DataFrames
	using LinearAlgebra
	using Printf
	using Memoize
        using Random
	using Distributions
	using DifferentialEquations
	using ParameterizedFunctions
	import Logging	
	import PlutoUI: combine

	Logging.disable_logging(Logging.Warn)
	md"""
	# Virtual Volatility-Humidity TDMA 

	Welcome to this educational notebook simulating the transfer through a tandem differential mobility analyzer system. This work was supported by the U.S. National Science Foundation grant AGS-2037704. 
	"""
end

# ╔═╡ 0128e0b7-1911-49f0-964c-f0a3f3c90e56
module EvaporationModel

using DifferentialEquations
using ParameterizedFunctions
using Memoize

ODEs = @ode_def_bare begin
	dDp = -G/Dp
end G

@memoize function evapmodel(Dd, G)
	function condition(u, t, integrator)
		out = (u[1] < 5e-9) ? 0.0 : 1.0
	end
	affect!(integrator) = terminate!(integrator)
	cb = ContinuousCallback(condition, affect!)

	prob = ODEProblem(ODEs, [Dd], (0.0, 3600.0), [G])
	sol = solve(prob, RK4(), abstol = 1e-20, reltol = 1e-20, callback = cb)
end

@memoize function getEFfunction(Ds, G)
	sols = map(D -> evapmodel(D, G), Ds)
	function EF(t)
		Dev = map(sols) do s
			out = s(t)[1] < 5e-9 ? 0.0 : s(t)[1]
		end

		return Dev./Ds
	end

	return EF
end

end

# ╔═╡ de2272eb-5476-4112-a94a-cf04847a263c
md"""

$(Resource("https://raw.githubusercontent.com/mdpetters/webapps/main/virtualDMA/notebooks/nsflogo.jpg", :width => 300)) 
		
This notebook is a "reactive" interactive application. After each change, the notebook revaluates all of the cells. This behavior is equivalent to that of a spreadsheet, where all cells are linked. Use show/hide code icon to the left view the code. If you have questions or comments, please send an email to  Markus Petters: [markus.petters@ucr.edu](mailto:markus.petters@ucr.edu).

#### Purpose

The main purpose of this resource is instructional. The instructor guides the students through the basic theory of operation of a tandem differential mobility analyzer system. Changing the inputs allows students to visualize the role of multiply charged particles in shaping the tandem DMA response function.    
"""

# ╔═╡ 93a3dafa-fafb-45d6-9cf0-53984e0bf56b
begin

	md"""
	# Tandem DMA Introduction 

	## DMA Basics
	The cylindrical differential mobility analyzer is an annular capacitor. The column's properties are defined by the radii $r_1$, $r_2$, the length of the aerosol path, $l$. Operation conditions are defined by the the electric potential $v$ applied across the annulus and the four flow rates: sheath flow, $q_{sh}$, polydisperse aerosol flow $q_a$, excess flow, $q_{ex}$, and monodisperse sample flow, $q_{sa}$. Throughout this work it is assumed that the flows are balanced, i.e. $q_{sh} = q_{ex}$ and $q_{sa} = q_a$. The two flows tracked are $q_{sh}$ and $q_{sa}$. Dimensions and flow range for one commonly used model are provided in Table 1 below. 
	$(Resource("https://raw.githubusercontent.com/mdpetters/webapps/main/virtualDMA/notebooks/dma.png", :width => 400))

	 $(Resource("https://raw.githubusercontent.com/mdpetters/webapps/main/virtualDMA/notebooks/table.png", :width => 400))

	"""
end

# ╔═╡ 6ad7bd8b-5c5e-492c-abfd-0703ee7ebdcf
begin
md"""

	!!! note
	
	For the cylindrical DMA and balanced flows, the relationship between voltage, mobility, and diameter is given by Knutson and Whitby (1975)
	
	$z^s = \frac{q_{sh}}{2\pi l v} \ln \left(\frac{r_2}{r_1}\right)$
	

	where $v$ is the potential applied between the inner and out section of the annulus, $r_1$, $r_2$, and $l$ are the dimensions of the cylindrical DMA (Table 1) and $q_{sh}$ is the sheath flow rate. The relationship between diameter $d_p$ and centroid mobility "z star" ($z^s$) is 

	$d_p =  \frac{kec_c}{3\pi \eta z^s}$ 
	
	where $e$ is the elementary charge and $k$ is the number of charges on the particle, $c_c$ is the Cunningham correction factor, and $\eta$ is the viscosity of the fluid.
	"""
end

# ╔═╡ 6e7f331c-76cf-42ee-95e3-1e21cca3f553
begin
	md"""
	## Tandem DMA
	
	In a typical tandem DMA, dried, charge equilibrated particles are classified in DMA 1. The flow is split between a condensation particle counter (CPC) and some conditioner. Conditioners may humidify the particles to grow them, heat the particles to evaporate them (thermodenuder), store the particles to evaporate them (isothermal evaporation), or pass them through an ion field to alter the charge distribution. The changed mobility distribution is measured using the second DMA that is operated in scanning or stepping mode. 

	$(Resource("https://raw.githubusercontent.com/mdpetters/webapps/main/virtualDMA/notebooks/tdma.png", :width => 500))

	## Bipolar charge distribution

	The bipolar charge conditioner (neutralizer) imparts a statistical charge distribution on the particles. Note that multiply charged particles are more prevalent at larger diameters. Also note that the charging of +1 and -1 is asymmetric.

	$(Resource("https://raw.githubusercontent.com/mdpetters/webapps/main/virtualDMA/notebooks/charge.png", :width => 900))

	
	"""
end

# ╔═╡ 963d9fae-2b4e-4067-8800-dc1e04cde1ba
@bind values confirm( 
	combine() do Child
		md"""
		### Input
		Please select the input. Once complete, press the **Submit** button to update the calculations. Note that check boxes evaluate immediately.

		**Aerosol size distribution:** The multi-modal lognormal size distribution is 
		
		 $\frac{dN}{d\ln D_p} = \sum_{i=1}^n \frac{N_{t,i}}{\sqrt{2\pi}\ln\sigma_{g,i}} \exp \left(- \frac{\left[\ln D_p-\ln D_{pg,i}\right]^2}{2\ln \sigma_{g,i}^2}\right)$  
		
		where $\frac{dN}{d\ln D_p}$ is the spectral number density, $N_{t,i}$ is the total number concentration, $\sigma_{g,i}$ is the geometric standard deviation, and $D_{pg,i}$ is the geometric mean diameter of the $i^{th}$ mode, and $n$ is the number of modes. The units are $N_{t,i}$ [$cm^{-3}$], $D_{g,i}$ [$nm$], and $\sigma_{g,i}$ [$-$]. The distribution is a superposition of multiple modes.
		
		| Mode | $N_{t,i}$ | $D_{g,i}$ | $\sigma_{g,i}$ |
		|----------- |----------- |----------- | ----------- | 
		| 1 | $(Child(NumberField(0:100:10000, default = 2000))) | $(Child(NumberField(5:1:1000, default = 200))) | $(Child(NumberField(1:0.05:3, default = 1.6))) |
		| 2 | $(Child(NumberField(0:100:10000, default = 0))) | $(Child(NumberField(5:1:1000, default = 50))) | $(Child(NumberField(1:0.05:3, default = 1.3))) |
		| 3 | $(Child(NumberField(0:100:10000, default = 0))) | $(Child(NumberField(5:1:1000, default = 600))) | $(Child(NumberField(1:0.05:3, default = 1.4))) |

		**Aerosol properties/mixing state:** Hygroscopicity parameter, growth/evaporation parameter, and number fraction for class 1, 2, and 3.

		The hygroscopic growth factor $g$ is computed from $a_w = RH/100\%$ and $\kappa$. The Kelvin effect is ignored.

		$g = \left (1 + \kappa \frac{a_w}{1-a_w}\right)^{1/3}$

		The evaporation is modeled via integration of 

		$\frac{dD_p}{dt} = \frac{G}{D_p}$

		where $G$ subsumes the diffusivity, vapor pressure and accomodation coefficient. Note that this standard growth law has a square root time dependence. Small particle grow/evaporate faster than larger particles.
		
		| Properties |   Class  1 | Class 2    | Class 3    | 
		|----------- |----------- |----------- |----------- |  
		| $\kappa (-)$   |  $(Child(NumberField(0:0.01:2.0, default = 0))) | $(Child(NumberField(0:0.01:2.0, default = 0.5))) | $(Child(NumberField(0:0.01:2.0, default = 1.3))) |
		| $log_{10} G$ | $(Child(NumberField(-20:0.5:-1, default = -20))) |  $(Child(NumberField(-20:0.5:-1, default = -17))) |  $(Child(NumberField(-20:0.5:-1, default = -16))) |
		| $f_i (-)$   | $(Child(NumberField(0:0.05:1, default = 0.0))) | $(Child(NumberField(0:0.05:1, default = 1.0))) | $(Child(NumberField(0:0.05:1, default = 0.0))) |

		**DMA configuration:** The units are $D$ [$nm$], $Q_{sh}$ [$L\;min^{-1}$], $Q_{sa}$ [$L\;min^{-1}$]. The neutralizer is a bipolar charger. The detector is a CPC with sample volume equal to the sample flow rate and sampling at 1 Hz.
		
		|DMA | $D$ | $Q_{sh}$ | $Q_{sa}$ | Neutralizer | Detector |
		|----------- |----------- |----------- | ----------- | ----------- |----------- |
		| 1 | $(Child(NumberField(1:1:1000, default = 100))) |  $(Child(NumberField(1:0.5:10, default = 5))) | $(Child(NumberField(0.5:0.5:5, default = 1))) | *always on* |
		| 2 | *scanning* | $(Child(NumberField(1:0.5:10, default = 5))) | $(Child(NumberField(0.5:0.5:5, default = 1))) | on/off $(@bind neutch CheckBox(default = false)) | noise on/off $(@bind pnoise CheckBox(default = false)) |

		**Instrument configuration**

		|     |            |            |            |
		|----------- |----------- |----------- |----------- |  
		| Num charges | $(Child(NumberField(1:1:10, default = 2)))   | Show all charges? $(@bind showmulti CheckBox(default = false))  | Full x-axis? $(@bind fullax CheckBox(default = false))  | 
		| Humidifier  | on/off $(@bind humch CheckBox(default = false))  | $RH(\%)$ =  $(Child(NumberField(10:5:95, default = 80)))  | | 
		| Volatility  | on/off $(@bind volch CheckBox(default = false))  | $t =$ $(Child(NumberField(1:1:3600, default = 200))) | |
		
		"""
	end
)


# ╔═╡ 4928d6ea-fe33-40e6-86f1-b7d1be5f236a
begin
	@memoize function initializeDMAs(Dd, k, q)
		t, p = 295.15, 1e5
		r₁, r₂, l = 9.37e-3, 1.961e-2, 0.44369
		Λ₁ = DMAconfig(t, p, q.sa1, q.sh1, r₁, r₂, l, 160.0, :-, 10, :cylindrical)
	    Λ₂ = DMAconfig(t, p, q.sa2, q.sh2, r₁, r₂, l, 0.0, :-, 10, :cylindrical)
	    δ₁ = setupDMA(Λ₁, dtoz(Λ₁, 5000e-9), dtoz(Λ₁, 8e-9), 64)   
		δ₂ = setupDMA(Λ₂, dtoz(Λ₂, 5000e-9), dtoz(Λ₂, 8e-9), 256)
	    Λ₁, Λ₂, δ₁, δ₂ 
	end
	
	function poisson_noise(Qcpc, N; seed = seed, t = 2.0)
		Q = Qcpc * 1e6   # Flow in cm3 s-1
	    c = N * Q * t   # number of counts in each bin
	
		map(c) do i
			f = rand(Poisson(i), 1)
			f[1] / (Q * t)
		end
	end
	
	@memoize function TDMAFunction(Dd, k, q, m; neutralizer = false) 
		
		Λ₁, Λ₂, δ₁, δ₂ = initializeDMAs(Dd, k, q)
		
		O(k) = mapfoldl(
			zs -> (δ₂.Ω(Λ₂, δ₂.Z, zs / k, k) .* δ₂.Tl(Λ₂, δ₂.Z, k))', vcat, δ₂.Z
		)
		
		T₁(zˢ, k) = δ₁.Ω(Λ₁, δ₁.Z, zˢ / k, k) .* δ₁.Tc(k, δ₁.Dp) .* δ₁.Tl(Λ₁, δ₁.Dp)
		DMA₁(𝕟, zˢ, gf) = @_ map((gf ⋅ (T₁(zˢ, _) * 𝕟)), 1:m)
		itp(𝕟) = interpolateSizeDistributionOntoδ((𝕟, δ₂))
		
		if neutralizer == false
			DMA₂ = (𝕟, k) -> O(k) * 𝕟
		else
			DMA₂ = (𝕟, k) -> δ₂.𝐀 * 𝕟
		end
		
		return Λ₁, Λ₂, δ₁, δ₂, T₁, DMA₁, itp, DMA₂
	end

	function mainplot(Ax, Dd, q, m)
		Λ₁, Λ₂, δ₁, δ₂, T₁, DMA₁, itp, DMA₂ = 
			TDMAFunction(Dd, 120, q, m; neutralizer = neutch) 
	
		𝕟ᶜⁿ = DMALognormalDistribution(Ax, δ₁)
		zˢ = dtoz(Λ₁, Dd);    
		
		if humch == false
			gf₁ = 1.0
			gf₂ = 1.0
			gf₃ = 1.0
		else
			aw = values[25]/100.0               
			gf₁ = (1 + values[10] * aw/(1-aw))^0.3333
			gf₂ = (1 + values[11] * aw/(1-aw))^0.3333
			gf₃ = (1 + values[12] * aw/(1-aw))^0.3333
		end

		if volch == false
			ef₁ = 1.0
			ef₂ = 1.0
			ef₃ = 1.0
		else
			t = values[26]
			ef₁ = EvaporationModel.getEFfunction(δ₁.Dp*1e-9, 10.0^values[13]) |> 
				x -> x(t)
			ef₂ = EvaporationModel.getEFfunction(δ₁.Dp*1e-9, 10.0^values[14]) |> 
				x -> x(t)
			ef₃ = EvaporationModel.getEFfunction(δ₁.Dp*1e-9, 10.0^values[14]) |> 
				x -> x(t)
			println(ef₁)
		end
		
		f₁ = values[16]
		f₂ = values[17]
		f₃ = values[18]

		ℕ₀ = [0.0*𝕟ᶜⁿ for i = 1:m]	
		ℕ₁ = (f₁ < 0.01) ? ℕ₀ : DMA₁(𝕟ᶜⁿ, zˢ, gf₁.*ef₁)      
	    ℕ₂ = (f₂ < 0.01) ? ℕ₀ : DMA₁(𝕟ᶜⁿ, zˢ, gf₂.*ef₂)
		ℕ₃ = (f₃ < 0.01) ? ℕ₀ : DMA₁(𝕟ᶜⁿ, zˢ, gf₃.*ef₃)
		ℕ = map((𝕟₁,𝕟₂,𝕟₃) -> f₁*𝕟₁ + f₂*𝕟₂ + f₃*𝕟₃, ℕ₁, ℕ₂, ℕ₃)
		
		𝕄ₐ = map(k -> (@> itp(ℕ[k]) DMA₂(k)), 1:m) 

		if pnoise == true
			𝕄 = map(𝕄ₐ) do 𝕟
				N = poisson_noise(q.sa2, 𝕟.N; seed = 703, t = 1.0)
				SizeDistribution(
					[[]], 
					𝕟.De, 
					𝕟.Dp, 
					𝕟.ΔlnD, 
					N./𝕟.ΔlnD, 
					N, 
					:perturbed
				)
				
			end
		else
			𝕄 = 𝕄ₐ
		end

		𝕞ᵗ = sum(𝕄)                               # total response

		mdf(k) = DataFrame(
		    Dp = 𝕄[k].Dp, 
		    S = 𝕄[k].S .* 𝕞ᵗ.ΔlnD, 
		    Dist = ["𝕄[$k]" for i = 1:length(𝕄[k].Dp)]
		)
		
		df1 = mapreduce(mdf, vcat, 1:m)
		df2 = DataFrame(
			Dp = 𝕞ᵗ.Dp, 
			S = 𝕞ᵗ.S .* 𝕞ᵗ.ΔlnD, 
			Dist = ["𝕞ᵗ" for i = 1:length(𝕞ᵗ.Dp)]
		)
		df = showmulti == true ? [df2; df1] : df2
		
		colors = ["black", "darkred", "steelblue3", "darkgoldenrod", "navajowhite"]
		xlabels = log10.([10, 30, 100, 300, 1000])
	
		if fullax == true
			Coordinates = Coord.cartesian(xmin = log10(8), xmax = log10(2000))
			lfun = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""
		else
			Coordinates = Coord.cartesian(
				xmin = log10(Dd*1e9/3), 
				xmax = log10(Dd*1e9*5)
			)
			lfun = x -> (x in xlabels) && (x > log10(Dd*1e9/3)) && 
				(x < log10(Dd*1e9*5)) ? @sprintf("%2i", exp10(x)) : ""
		end
	
		z = dtoz(Λ₂, Dd)
		Ds = map(k->ztod(Λ₂,k,z), 1:m)
		dfx = map(1:1) do x
			layer(
				x = [Ds[x], Ds[x]], 
				y = [0, maximum(𝕞ᵗ.N)], 
				Geom.line, 
				Theme(default_color=colors[x+1])
			)
		end
			
		p1 = plot(
			dfx...,
		    layer(df,
		    x = :Dp,
		    y = df[!,:S],
		    color = :Dist,
		    Geom.line),
		    Guide.xlabel("Apparent mobility diameter (nm)", 
				orientation = :horizontal),
		    Guide.ylabel("N (cm⁻³)"),
			Guide.xticks(ticks = log10.([[8,9];10:10:90;100:100:1000;[2000]])),
		    Guide.colorkey(; title = ""),
		    Scale.color_discrete_manual(colors...),
			Scale.x_log10(labels = lfun),
			Coordinates,
		    Theme(plot_padding = [5mm, 5mm, 0mm, 0mm]),
		)	
		p0 = pdfplot(𝕟ᶜⁿ, Dd, Λ₁,m)
		
		println("Executed Main Plot")

		return hstack(p0, p1)
	end


	function pdfplot(𝕟₁, D, Λ, m)
	    df1 = DataFrame(Dp = 𝕟₁.Dp, S = 𝕟₁.S, Dist = "Size dist.")
	    z = dtoz(Λ, D)
		Ds = map(k->ztod(Λ,k,z), 1:m)
		knots = reverse(𝕟₁.Dp)
		Interpolations.deduplicate_knots!(knots)
		itp = interpolate((knots,), reverse(𝕟₁.S), Gridded(Linear()))
		extp = extrapolate(itp, 0)
		dfx = map(1:m) do x
			y = [0, extp(Ds[x])]
			DataFrame(Dp = [Ds[x], Ds[x]], S = y, Dist = "+$x charge")
		end
		
		df = vcat(df1, dfx...)
	    xlabels = log10.([10, 30, 100, 300, 1000])
	    colors = ["black", "darkred", "steelblue3", "darkgoldenrod", "navajowhite"]
	    
	    return plot(
	        df,
	        x = :Dp,
	        y = :S,
	        color = :Dist,
	        Geom.line,
	        Guide.xlabel("Mobility diameter (nm)"),
	        Guide.ylabel("dN/dlnD (cm⁻³)"),
	        Guide.xticks(ticks = log10.([[8,9];10:10:90;100:100:1000;[2000]])),
	        Guide.colorkey(; title = ""),
	        Scale.x_log10(
				labels = x -> x in xlabels ? @sprintf("%2i", exp10(x)) : ""
			),
	        Scale.color_discrete_manual(colors...),
	        Coord.cartesian(xmin = log10(8), xmax = log10(2000)),
			Theme(plot_padding = [5mm, 0mm, 0mm, 0mm], default_color="black")
	    )
	end

end

# ╔═╡ 188d473e-17e5-443b-bc66-8d8ffa25a878
begin
Dd = values[19]	* 1.0e-9
charges = values[24]

lpm = 1.666e-5

q = (
	sh1 = values[20]*lpm, 
	sa1 = values[21]*lpm, 
	sh2 = values[22]*lpm, 
	sa2 = values[23]*lpm
)
	
Ax = [
	[values[1], values[2], values[3]],
	[values[4], values[5], values[6]],
	[values[7], values[8], values[9]]
]

set_default_plot_size(18cm, 8cm)
aaaa = mainplot(Ax, Dd, q, charges)
end

# ╔═╡ Cell order:
# ╟─dc1a7d9c-fad9-11eb-3360-4be773edd179
# ╟─de2272eb-5476-4112-a94a-cf04847a263c
# ╟─93a3dafa-fafb-45d6-9cf0-53984e0bf56b
# ╟─6ad7bd8b-5c5e-492c-abfd-0703ee7ebdcf
# ╟─6e7f331c-76cf-42ee-95e3-1e21cca3f553
# ╟─963d9fae-2b4e-4067-8800-dc1e04cde1ba
# ╟─188d473e-17e5-443b-bc66-8d8ffa25a878
# ╟─0128e0b7-1911-49f0-964c-f0a3f3c90e56
# ╟─4928d6ea-fe33-40e6-86f1-b7d1be5f236a
