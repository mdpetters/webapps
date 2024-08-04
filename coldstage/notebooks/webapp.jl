### A Pluto.jl notebook ###
# v0.19.9

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

# ╔═╡ cc5fa8b8-4da3-11ee-28fd-7dacab544db3
begin
	using Pkg
	Pkg.activate(Base.current_project())
	using DifferentialEquations
	using Plots
	using Plots.PlotMeasures
	using PlutoUI
	import PlutoUI:combine
	using Logging

	Logging.disable_logging(Logging.Warn)

	#https://mdpetters.github.io/cee200/assets/
	cs_url = "https://mdpetters.github.io/cee200/assets/cold_stage.png"
	tec_url = "https://mdpetters.github.io/cee200/assets/tec_performance.png"
	control1_url = "https://mdpetters.github.io/cee200/assets/control1.svg"
	control2_url = "https://mdpetters.github.io/cee200/assets/control2.svg"
	control3_url = "https://mdpetters.github.io/cee200/assets/control3.svg"
	control4_url = "https://mdpetters.github.io/cee200/assets/control4.svg"
	control5_url = "https://mdpetters.github.io/cee200/assets/TC-36-25_RS232_02.jpg"

	md"""
$(TableOfContents(depth=4))
#  Introduction 

## The Ice Nucleation Cold Stage
Welcome to this interactive introduction into the ice nucleation cold stage control system. The ice nucleation cold stage instrument for measuring the concentration of ice nucleating particles in liquids (Mahant et al., 2023). The technique has been widely used in the ice nucleation literature (e.g. DeMott et al., 2018, 2024, Yadav et al, 2019). A news article describing the instrument can be found [here](https://www.cert.ucr.edu/news/2024/06/24/ce-cert-team-collaborates-west-texas-am-low-temperature-research). 

This page is a Pluto notebook. After each change, the reactive Pluto notebook 		revaluates all of the cells. This behavior is equivalent to that of a spreadsheet 	application. The purpose of the notebook is to explain the temperature control system of the instrument. The notebook is suitable for coursework teaching about control systems, or in coursework discussion ice nucleation instrumentation. 

$(Markdown.MD(Markdown.Admonition("note", "Author", [md" 
If you have questions or comments, please send an email to\
	Markus Petters: **markus.petters@ucr.edu**"])))

## Acknowledgements

$(Resource("https://raw.githubusercontent.com/mdpetters/webapps/main/virtualDMA/notebooks/nsflogo.jpg", :width => 300)) 

The instrument was developed and characterized using funding from the National Science Foundation awards NSF-AGS 1010851 and NSF-AGS 1450690. This notebook and the infrastructure to host the notebook was supported via NSF-AGS 2112978. 

## References

DeMott, PJ, JA Mirrielees, SS Petters, DJ Cziczo, MD Petters, HG Bingemer, TCJ Hill, K Froyd, S Garimella, AG Hallar, EJT Levin, IB McCubbin, AE Perring, CN Rapp, T Schiebel, J Schrod, KH Suski, D Weber, MJ Wolf, M Zawadowicz, J Zenker, O Möhler, SD Brooks, The Fifth International Workshop on Ice Nucleation Phase 3 (FIN-03): Field Intercomparison of Ice Nucleation Measurements, Atmospheric Measurement Techniques - Discussions, 2024.

S Mahant, S Yadav, C Gilbert, ER Kjærgaard, MM Jensen, T Kessler, M Bilde, MD Petters, An open-hardware community ice nucleation cold stage for research and teaching, HardwareX, e00491, https://doi.org/10.1016/j.ohx.2023.e00491, 2023.

Yadav, S., Venezia, R. E., Paerl, R. W., & Petters, M. D.: Characterization of ice-nucleating particles over Northern India. Journal of Geophysical Research: Atmospheres, 124, 10467–10482. https://doi.org/10.1029/2019JD030702, 2019. 

DeMott, P. J., Möhler, O., Cziczo, D. J., Hiranuma, N., Petters, M. D., Petters, S. S., Belosi, F., Bingemer, H. G., Brooks, S. D., Budke, C., Burkert-Kohn, M., Collier, K. N., Danielczok, A., Eppers, O., Felgitsch, L., Garimella, S., Grothe, H., Herenz, P., Hill, T. C. J., Höhler, K., Kanji, Z. A., Kiselev, A., Koop, T., Kristensen, T. B., Krüger, K., Kulkarni, G., Levin, E. J. T., Murray, B. J., Nicosia, A., O'Sullivan, D., Peckhaus, A., Polen, M. J., Price, H. C., Reicher, N., Rothenberg, D. A., Rudich, Y., Santachiara, G., Schiebel, T., Schrod, J., Seifried, T. M., Stratmann, F., Sullivan, R. C., Suski, K. J., Szakáll, M., Taylor, H. P., Ullrich, R., Vergara-Temprado, J., Wagner, R., Whale, T. F., Weber, D., Welti, A., Wilson, T. W., Wolf, M. J., and Zenker, J.: The Fifth International Workshop on Ice Nucleation phase 2 (FIN-02): laboratory intercomparison of ice nucleation measurements, Atmos. Meas. Tech., 11, 6231–6257, https://doi.org/10.5194/amt-11-6231-2018, 2018. 
"""
end

# ╔═╡ aabae4cb-3912-408d-bacf-267f2440c047
Markdown.MD(
	Markdown.Admonition("warning", "Key Concept", [md"
A control system consists of components and circuits that work together to maintain the process at a desired operating point. Every home or an industrial plant has a temperature control that maintains the temperature at the thermostat setting. In industry, a control system may be used to regulate some aspect of production of parts or to maintain the speed of a motor at a desired level.
"]))

# ╔═╡ 9db0e27f-32f5-4b03-9d70-5dea9f079ecb
md"""
# Thermolectric Cold-Stage

## Introduction
This module will develop the design and analysis of simple control systems using a simple example application: an enclolsed thermoelectrically cooled/heated plate surface whose temperature can be varied between −45 °C and +90 °C.  

$(Resource(cs_url, :width => 2500px))

**Figure.** A  thermoelectric cold plate with an enclousre cell. 1: heat sink, 2: barb connectors, 3: TEC module, 4: base block, 5: cold plate, 6: spacers, 7: M2.5 screw, 8: M4 screw, 9: glass lid, T: thermistor opening. Two thermistor openings are located on the left front facing side of the baseblock. (Source: Mahant et al., 2023).

The purpose of the heat sink is the maintain one side of the TEC module at a constant temperature. 

## Instrument Model
The temperature at the metal surface of the instrument can be modeled using [Newton's law of cooling](https://en.wikipedia.org/wiki/Newton%27s_law_of_cooling).

```math
\frac{dT}{dt} = k(T_{env}-T) - \frac{Q(V,\Delta T)}{mc}
```

where ``k`` is the coefficient of heat transfer, ``T_{env}`` is the temperature of the environment, ``T`` is the temperature of the place, ``Q`` is the heat transfer across the TEC, ``m`` is the mass of the metal stage, and ``c`` is the heat capacity of the metal. 

Performance of a [thermoelectric element](https://en.wikipedia.org/wiki/Thermoelectric_generator) can be modeled via

```math
Q(V,\Delta T) = a_1\sqrt{V} + a_2\Delta T + a_3 
```

where  ``V`` is the voltage applied to the TEC module, ``\Delta T`` is the temperature difference between the hot and cold side, and ``a_i`` are fitted coefficients that determine the performance for a specific TEC model.  

$(Resource(tec_url, :height => 900px))
**Figure**. Example relationship via heat removed vs. temperature gradient for different applied voltages. Symbols show data and the line corresponds to the multilinear model.

"""

# ╔═╡ ba749db8-c015-4137-ac7e-92fc0bfbe7ce
Markdown.MD(
	Markdown.Admonition("info", "Exercise", [md"
Based on the performance characteristic of the TEC module. Predict the coldest temperature that the stage can reach assuming that the heat sink temperature is 280K. 
	"]))

# ╔═╡ a25d1d5d-abc8-4ec4-80b2-ecdaea755a1c
md"""
## State Space Representation

The instrument performance can alternatively described by the canonical [state-space model](https://en.wikipedia.org/wiki/State-space_representation) representation:

```math
\begin{eqnarray}
 & \quad \vec{q}'(t) &= \mathbf{A}\vec{q}(t) + \mathbf{B} \vec{f}(t) \\
 & \quad \vec{y}(t) &=  \mathbf{C}\vec{q}(t) + \mathbf{D} \vec{f}(t) \\
\end{eqnarray}
```

where ``\vec{q} = [T]``, ``\quad \vec{y} = [T]``, and  ``\quad\vec{f} = 
\begin{bmatrix}
\sqrt{V} & T_{hot} & 1 & T_{env}
\end{bmatrix}
``

And the coefficients:

``\quad A = [k + \frac{a_2}{mc}]``,  ``\quad B = \begin{bmatrix} -\frac{a_1}{mc} \\ -\frac{a_2}{mc} \\ -\frac{a_3}{mc} \\ k \end{bmatrix}``, ``\quad C = [1]``, and ``D = \begin{bmatrix} 0 \\ 0 \\ 0 \\ 0 \end{bmatrix}``

Thus we can use a generic state space solver to model the system. 
"""

# ╔═╡ be1a604c-6cf6-405f-9fa0-86e3b454315e
function state_space_solver(A, B, C, D, q0, f, tspan)
	h(q, p, t) = A*q + B*f(t)  # canonical equation 1
	y(q, t) = C*q + D*f(t)     # canonical equation 2
	
	problem = ODEProblem(h, q0, tspan, Float64[]) 
	solution = solve(problem, RK4(), reltol = 1e-12, abstol = 1e-12)
	return map(y, solution.u, solution.t), solution.t
end

# ╔═╡ 759f6f53-866c-4cc8-9bb5-41be3b43f1dc
md"""
## Solution for Constant V

The function below intializes the problem with known coefficients for a constant input voltage V. The stage cools and reaches an equilibrium temperature that is controlled by the heat sink temperature and the heat loss rate to the environment. 
"""

# ╔═╡ bdc5a009-c6c4-434d-a619-63081c6af32d
function solve_const_v(V)
	m = 0.3                          # kg
	c = 502.0                        # J kg-1 K-1
	k = 300.0*(0.025*0.025)/(m*c)    # s^-1
	a1 =  11.07 				     # W V^-0.5 
	a2 =  -0.352                     # W K^-1
	a3 =  -8.29                      # W

	u0 = [300.0]
	tspan  = [0.0, 3000.0]                   # 0 to 1000 s
	f(t) = [sqrt(V), 280.0, 1.0, 300.0] # V = 5V, Th = 280K, 1, Tenv = 300K
	
	A = -k + a2/(m*c)
	B = [-a1/(m*c) -a2/(m*c) -a3/(m*c) k]
	C = 1.0
	D = [0.0 0.0 0.0 0.0]
	
	state_space_solver(A, B, C, D, u0, f, tspan)	
end

# ╔═╡ b8fe4424-9eac-4011-85ba-f1bf5e72c718
Markdown.MD(
	Markdown.Admonition("info", "Exercise", [md"
- Compare the model solution against your earlier prediction. 
- Does the equilibration time depend on the voltage?
- Using the same equipment, could you imaging a faster path to cool to the equilibrium temperature?
	"]))

# ╔═╡ 1b1dbf06-c67b-4447-9049-4ab1bf82ee8a
@bind myV combine() do Child
	md"""
	Voltage:  $(
		Child(Slider(0:1:12, default = 6, show_value = true))
	) 
	"""
end

# ╔═╡ f4c5f8dd-a66f-4b1a-8a6f-3e4ad9615999
outvec = solve_const_v(myV[1])

# ╔═╡ cf63fd31-4c00-46e5-8648-713b59080eba
begin
	pos3 = map(x -> x[1], outvec[1])
	ts3 = outvec[2]
	plot(ts3/60.0, pos3, color = :black, label = "RK4 Solution", lw = 1, 
		xlabel = "Time (min)", ylabel = "T (K)", size = (700, 300), bottom_margin = 20px, left_margin = 10px, minorgrid = :true, framestyle = :box)
end

# ╔═╡ 451a1910-4ce9-4e15-bca3-c516e20e5bd5
md"""
## Solution for Time Varying V

The function below intializes the problem with known coefficients for a time-varying input voltage V. Unsurprisingly, the changing voltage is mirrored in the changing temperature around some average temperature.
"""

# ╔═╡ c42f51a3-a02d-4c86-beee-7feedd378ad0
function solve_variable_v(fv)
	m = 0.3                          # kg
	c = 502.0                        # J kg-1 K-1
	k = 300.0*(0.025*0.025)/(m*c)    # s^-1
	a1 =  11.07 				     # W V^-0.5 
	a2 =  -0.352                     # W K^-1
	a3 =  -8.29                      # W

	
	u0 = [300.0]
	tspan  = [0.0, 3000.0]                   # 0 to 1000 s
	f(t) = [fv(t), 280.0, 1.0, 300.0]    
	
	A = -k + a2/(m*c)
	B = [-a1/(m*c) -a2/(m*c) -a3/(m*c) k]
	C = 1.0
	D = [0.0 0.0 0.0 0.0]
	
	state_space_solver(A, B, C, D, u0, f, tspan)	
end

# ╔═╡ 538e7d53-5cc3-4447-b1c7-3d02a0946f09
begin
	fex(t) = (sin(0.01*t) + 1.0)/2.0*sqrt(3) + 1
	outvec1 = solve_variable_v(fex)
	ts4 = outvec1[2]
	p1 = plot(ts4./60, fex.(ts4), color = :black, label = :none,
		ylabel = "V")
	pos4 = map(x -> x[1], outvec1[1])
	
	p2 = plot(ts4/60.0, pos4, color = :black, label = "RK4 Solution", 
		lw = 1, xlabel = "Time (min)", ylabel = "T (K)")
	plot(p1, p2, layout = grid(2,1))
end

# ╔═╡ 6369c1ac-e52c-4d81-869d-77cbd483fd8b
md"""
## Block Diagram Representation

Block diagrams consist of a single block or a combination of blocks. These are used to represent the control systems in pictorial form.

Recall the state-space equations:

```math
\begin{eqnarray}
 & \quad \vec{q}'(t) &= \mathbf{A}\vec{q}(t) + \mathbf{B} \vec{f}(t) \\
 & \quad \vec{y}(t) &=  \mathbf{C}\vec{q}(t) + \mathbf{D} \vec{f}(t) \\
\end{eqnarray}
```

We have the input vector ``\vec{f}``, the output vector ``\vec{y}``. The system is subject to initial conditions ``\vec{q_0}``. The `process` box solves for the evolution of the system over time ``t``.

$(Resource(control1_url, :width => 1500px))

For the cold stage, the inputs are ``\quad\vec{f} = 
\begin{bmatrix}
\sqrt{V} & T_{hot} & 1 & T_{env}
\end{bmatrix}
``, the process state variables are ``\vec{q} = [T]``, and the outputs are ``\quad \vec{y} = [T]``. Therefore the intial input ``\vec{q}_0`` is the intial temperature of the stage ``T_0``.
"""

# ╔═╡ 03a0fc10-9875-4248-8711-8a3887db30a5
md"""
# Discrete Event System Simulation

## Introduction
A discrete-event simulation (DES) models the operation of a system as a (discrete) sequence of events in time. Each event occurs at a particular instant in time and marks a change of state in the system. Between events, the system is allowed to evolve according to some set of rules. Here we simulate the system as a data acquisition system would perceive it. At each update event (determined by the sampling frequency), the data acquisition system reads the state using sensors, here the temperature of the cold stage, and perhaps other states such as the environmental temperature or the heat sink temperature. A typical data acquisition system would then write the data to file. The user may also change inputs to the system, for example the voltage supplied to the TEC module. 

## Block Diagram

Instead of solving for the evolution over the entire time domain, the discrete simulation soves for the evolution from time ``t \rightarrow t + \Delta t``.

$(Resource(control2_url, :width => 1500px))

To simulate the evolution of the system over longer times, the output of ``\vec{y}(t + \Delta t)`` is used to re-initialize the system. For the cold stage, the output is ``T`` and maps directly to the initial condition ``T_0``. This introduces: 

- External Observer: Every ``\Delta t`` we observe the state of the system.
- External Controller: Every ``\Delta t`` we can, in principle, update the system forcing vector ``\vec{f}``. 

"""

# ╔═╡ 8aa24500-d236-4e24-bf06-4db94e612a50
md"""

## Implementation

### ODE Solution 
The `update_system` function incrementally solves the differential equation in `dt` intervals, where `dt` is the time step. 
"""


# ╔═╡ 54a7e106-4537-4857-a9e7-fdc3e9bf4ba1
function update_system(T, V, dt; Thot = 280.0, Tenv = 300.0)
	m = 0.3                          # kg
	c = 502.0                        # J kg-1 K-1
	k = 200.0*(0.025*0.025)/(m*c)    # s^-1
	a1 =  11.07 				     # W V^-0.5 
	a2 =  -0.352                     # W K^-1
	a3 =  -8.29                      # W

	u0 = [T]
	tspan  = [0.0, dt]               # 0 to 1 s
	sig = sign(V)                    # Needed for negative V
	f(t) = [sig.*sqrt(abs(V)), Thot, 1.0, Tenv]    
	
	A = -k + a2/(m*c)
	B = [-a1/(m*c) -a2/(m*c) -a3/(m*c) k]
	C = 1.0
	D = [0.0 0.0 0.0 0.0]
	
	out = state_space_solver(A, B, C, D, u0, f, tspan)
	return out[1][end][]
end

# ╔═╡ f0ecf625-53f3-4595-8764-ab0a95c5b6e3
# example: new temperature for -3V after 1 second
update_system(300.0, -3.0, 1.0)

# ╔═╡ 17c3d4a5-048a-44c1-a5f0-5278a2322f60
md"""

### Simulation 
We can then perform a simulation where we step in 1s increments for n seconds. After each 1s increment, the temperature is updated. Note that `update_system` still used the ODE integration under the hood which may be at a higher time resolution. The trick here is that we observe the system every second and have the opportunity to intervene, for example by changing the voltage. Without intervention, the result should be the same as running the ODE solver for the same time interval. 

Note that the simulation is generic. It calls a black box "update_system", that can be replaced with any system imaginable.
"""

# ╔═╡ 25b351fa-69e6-482e-9543-12e1dde9e438
function simulation(T, V, n; dt = 1.0)
	ts = Float64[]
	Ts = Float64[]
	Vs = Float64[]
	push!(Ts, T)
	push!(Vs, V)
	push!(ts, 0.0)
	for i = 1:n
		theT = Ts[end]
		theV = Vs[end]
		newT = update_system(theT, theV, dt)
		push!(Ts, newT)
		push!(Vs, theV)
		push!(ts, ts[end] + dt)
	end

	return ts, Ts, Vs
end

# ╔═╡ 1163c920-d2ca-42e7-9b9e-1e19d494bb46
# Example simulation
ts, Ts, Vs = simulation(300.0, 2.0, 360; dt = 10)

# ╔═╡ 43c38357-4d59-457d-9b90-94ae1f267212
plot(ts./60, Ts, color = :black, label = "0.1 Hz Discrete Event Simulation", 
		xlabel = "Time (min)", ylabel = "T (K)", size = (700, 300), bottom_margin = 20px, left_margin = 10px, minorgrid = :true, framestyle = :box)

# ╔═╡ 7d030c04-ac3b-4832-a3e8-6053fc9c14e7
md"""

# Controlling the Temperature

## Controller Block Diagram

$(Resource(control3_url, :width => 2800px))

A controller can be used to control the signal. The setpoint of the controller is ``r``. The output of the process is subtracted from the setpoint to produce the error ``e(t)``. The error function is then passed into the controller, which outputs the input ``f(t)`` for the process to be controlled. The output of the controller is also sometimes referred to as the *manipulated variable*.

In the subsequent section we will focus on single-input and single-output (SISO) systems. Such systems a single-variable control system with one input and one output. Hence the vector notation for ``f(t)`` and ``y(t)`` has been dropped for simplicity. Further note that when employing for the discrete event simulation approach, the process initial condition is updated each time step. In a real-world system this step is of course implicit. For this reason, the initial condition arrow is marked as a dashed line. The controller has a few import properties. It has 

- knowledge about the *current error*
- knowledge about the *past error* 
- **no knowledge** about the process or *future error*
"""

# ╔═╡ 1ae8b353-df38-430f-b220-401adddac1d1
md"""
## On-Off Controller 
### Description
On-Off control is the simplest form of feedback control. An on-off controller simply drives the manipulated variable from fully closed to fully open depending on the position of the controlled variable relative to the setpoint. A common example of on-off control is the temperature control in a domestic heating system. When the temperature is below the thermostat setpoint the heating system is switched on and when the temperature is above the setpoint the heating switches off.

Unfortunately, if the heating switches on and off the instant the measured temperature crossed the setpoint then the system would chatter – repeatedly switch on and off at very high frequency. This, in turn would substantially reduce the lifetime of the unit. To avoid chattering, practical on-off controllers usually have a deadband around the setpoint. When the measured value lies within this dead-band the controller does nothing. Only when the value moves outside the deadband, the controller switches. The effect of this is to introduce continuous oscillation in the value of the controlled variable – the large the dead-band the higher the amplitude and lower the frequency.

### Implementation

An example implementation of the on-off controller is as follows. 

```math
e(t) = r-y(t)
```
```math
f(t) = 
\begin{cases}
 +1  & (e(t)-d) < 0 \land (e(t-\Delta t) - d) > 0 \\ 
 -1 & (e(t)+d) > 0 \land (e(t-\Delta t) + d) < 0
\end{cases}
```

where ``d`` is the deadband value. The output is either ``-1`` (full cooling) or ``+1`` full heating. These values may need to scaled by the output range of the input control, for the TEC to ``\pm 12V``.

For example, if the current output ``y = 15``, the setpoint ``r = 10`` we need to reduce ``y`` and keep the process outut at full cooling. Only when the ``y < 9``, i.e. outside the deadband the signal needs to shift to warming. However, we also only want to switch if we crossed the critical line, hence we need to reach into the past and make sure that the condition changed. The same logic is applied on the other end. We need to heat until we cross the ``y = 11`` before switching back. The function below implements an off control in the context of the simulation. 
"""

# ╔═╡ f85d4b96-e5d8-4d5c-aaed-df30c8c06084
function simulation_onoff(Tenv, n; Tset = 280.0, Tdead = 2.0, Thot = 280.0)
	dt = 1.0 
	# Need to know if we need heat or cool initially
	Vstart = (Tset < Tenv) ? 12.0 : -12.0
	ts = Float64[]  # Memory array for t
	Ts = Float64[]  # Memory array for T
	Vs = Float64[]  # Memory array for V
	
	push!(Ts, Tenv)    # Need to initialize two values for the control loop to work
	push!(Ts, Tenv)    # Need to initialize two values for the control loop to work
	
	push!(Vs, Vstart)
	push!(Vs, Vstart)
	
	push!(ts, 0.0)
	push!(ts, 1.0)

	# Control loop
	for i = 1:n-1
		theT = Ts[end]
		theV = Vs[end]

		et = Tset - Ts[end]     # current error
		etp = Tset - Ts[end-1]  # previous error

		if (et + Tdead) < 0 && (etp + Tdead) > 0
			theV = +12.0 # + 1 is scaled to + 12 V
		end
		
		if (et - Tdead) > 0 && (etp - Tdead) < 0
			theV = -12.0 # - 1 is scaled to + 12 V
		end
		
		newT = update_system(theT, theV, dt; Thot = Thot, Tenv = Tenv)
		push!(Ts, newT)
		push!(Vs, theV)
		push!(ts, ts[end] + dt)
	end

	return ts, Ts, Vs
end

# ╔═╡ 0c6f2477-1252-4485-9596-2748a2508e14
@bind on_off combine() do Child
	md"""
	``T_{set}\;(K)`` $(
		Child(Slider(200:10:400, default = 250, show_value = true))
	) ``T_{env}\;(K)`` $(
		Child(Slider(280:10:320, default = 300, show_value = true))
	) 
	
	``T_{hot}\; (K)`` $(
		Child(Slider(260:10:320, default = 280, show_value = true))
	)  ``T_{dead}\; (K)`` $( 
		Child(Slider(0:0.1:4, default = 2, show_value = true))
	)
	"""
end

# ╔═╡ 77c08e89-b5ca-45d3-8916-3ac8cf9453cb
let 
	Tset = on_off[1]
	Tenv = on_off[2]
	Thot = on_off[3]
	Tdead = on_off[4]
	ts, Ts, Vs = simulation_onoff(Tenv, 1000; Tset = Tset, Thot = Thot, Tdead = Tdead)
	p1 = plot(ts./60.0, Vs, label = :none, color = :black, ylabel = "V")
	p2 = plot(ts./60.0, Ts, color = :black, ylabel = "T (K)", label = "T system")
	plot!([0, ts[end]./60.0], [Tset, Tset], label = "T set", color = :darkred)
	plot!([0, ts[end]./60.0], [Tset-Tdead, Tset-Tdead], l = :dash, label = "±T dead", 
		color = :darkred)
	plot!([0, ts[end]./60.0], [Tset+Tdead, Tset+Tdead], l = :dash, label = :none, 
		color = :darkred)
	plot(p1, p2, layout = grid(2,1))
end

# ╔═╡ 6647a86f-7f6d-43d7-892e-724d3a191ca9
md"""

## Proportional (P) Control

$(Resource(control3_url, :width => 2800px))

### Description

The controller output is proportional to the error signal, which is the scaled difference between the setpoint and the process variable. 

```math
g(t) = K_p \frac{r - y(t)}{s} + p_0 
```

where ``K_p`` is proportional gain, ``r`` is the setpoint, ``y(t)`` is the output variable, ``s`` is the span, and ``p_0`` is the controller output with zero error. The purpose of the span is to set the output to maximum when the difference between the setpoint and the process variable exceeds the span. Furthermore, in practical application the output needs to be constrained between ``-1`` and ``1``.

```math
f(t) = \begin{cases}
-1 &  g(t) < -1 \\
1 & g(t) > 1 \\
g(t) & else \\
\end{cases}
```

The case statement ensures that output ranges from ``-1`` (here maximum heating) to ``+1`` (here maximum cooling). Finally the output is scaled to ``\pm 12V``, which corresponds to the operating voltages of the TEC element. 

### Implementation
"""

# ╔═╡ b40aefaa-4600-42f7-ab1d-d92b54fb6401
function simulation_P(T, n; Tset = 280.0, Kp = 1.4, span = 20.0, p0 = 0.0)
	dt = 1.0 
	ts = Float64[]
	Ts = Float64[]
	Vs = Float64[]
	push!(Ts, T)
	push!(Vs, 12.0)
	push!(ts, 0.0)
	for i = 1:n
		theT = Ts[end]
		theV = Vs[end]
		g =  Kp*(theT - Tset)/span + p0
		
		if g < -1
			MV = -1
		elseif g > 1
			MV = 1
		else
			MV = g
		end
		
		theV = 12.0*MV # The voltage changes between -12 and 12  
		newT = update_system(theT, theV, dt)
		push!(Ts, newT)
		push!(Vs, theV)
		push!(ts, ts[end] + dt)
	end

	return ts, Ts, Vs
end

# ╔═╡ a20f135e-f061-4287-9e11-3ffad5cc64b8
@bind P combine() do Child
	md"""
	``T_{set}`` $(
		Child(Slider(200:10:400, default = 250, show_value = true))
	) 
	
	``T_{span}`` $( 
		Child(Slider(1:1:40, default = 20, show_value = true))
	) ``K_p`` $( 
		Child(Slider(0:0.1:4, default = 2, show_value = true))
	) ``p_0`` $( 
		Child(Slider(-1:0.1:1, default = 0, show_value = true))
	) 

	
	"""
end

# ╔═╡ 177d0f02-9045-49d5-84cd-bc7b6ebb58fc
let 
	Tset = P[1]
	span = P[2]
	Kp = P[3]
	p0 = P[4]
	ts, Ts, Vs = simulation_P(300.0, 3600; Tset = Tset, Kp = Kp, span = span, p0)
	p1 = plot(ts./60.0, Vs, label = :none, color = :black, ylabel = "V")
	p2 = plot(ts./60.0, Ts, color = :black, ylabel = "T (K)", label = "T")
	plot!([0, ts[end]./60.0], [Tset, Tset], label = "Tset")
	plot(p1, p2, layout = grid(2,1), size= (700, 400))
end

# ╔═╡ 2a2f1d95-60d1-41b5-9d94-e69e7c5bba78
md"""

### Advantages

The proportional contontroller gradually changes the input to converge on the set point, which eliminates chatter. There are also no oscillations around the final value. 

### Disadvanges

The P-controller suffers from offset error. Offset error is the difference between the desired value and the actual value. Over a range of operating conditions, proportional control alone is unable to eliminate offset error, as it requires an error to generate an output adjustment. Although a proportional controller may be tuned by adjusting ``p_0`` for a specific operating condition, it cannot be eliminated over the entire state space. 
"""

# ╔═╡ 9e35dc47-ff08-4254-96dd-fe27c0c7366d
md"""
## Proportional-Integral (PI) Control

$(Resource(control3_url, :width => 2800px))

### Description

Including an integral term increases action in relation not only to the error but also the time for which it has persisted. The integral term can be used to eliminate the bias in the P-controller. 



```math
e(t) = \frac{r - y(t)}{s}
```

```math
g(t) = K_p e(t) + K_i \int_0^t e(\tau) d\tau  
```

where ``K_p`` is proportional gain, ``r`` is the setpoint, ``y(t)`` is the output variable, ``s`` is the span (these terms are the same as for the P-controller), and ``K_i`` is the the integral gain.

As with the P-controller, the purpose of the span is to set the output to maximum when the difference between the setpoint and the process variable exceeds the span. Furthermore, in practical application the output needs to be constrained between ``-1`` and ``1``.

```math
f(t) = \begin{cases}
-1 &  g(t) < -1 \\
1 & g(t) > 1 \\
g(t) & else \\
\end{cases}
```

The case statement ensures that output ranges from ``-1`` (here maximum heating) to ``+1`` (here maximum cooling). Finally the output is scaled to ``\pm 12V``, which corresponds to the operating voltages of the TEC element. 

### Implementation
"""

# ╔═╡ 6e2bc9b4-6b61-49c8-ae72-0bc94e72543a
function simulation_PI(T, n; Tset = 280.0, Kp = 1.4, span = 20.0, Ki = 1.0)
	dt = 1.0 
	ts = Float64[]
	Ts = Float64[]
	Vs = Float64[]
	integral = Float64[]
	push!(Ts, T)
	push!(Vs, 12.0)
	push!(ts, 0.0)
	push!(integral, 0.0)
	for i = 1:n
		theT = Ts[end]
		theV = Vs[end]
		error = (theT - Tset)/span
		theIntegral = integral[end] + error*Ki*dt
		MV = Kp*error + Ki*theIntegral
		
		if MV < -1
			MV = -1
		elseif MV > 1
			MV = 1
		end
	
		theV = 12.0*MV
		newT = update_system(theT, theV, dt)
		push!(Ts, newT)
		push!(Vs, theV)
		push!(integral, theIntegral)
		push!(ts, ts[end] + dt)
	end

	return ts, Ts, Vs, integral
end

# ╔═╡ bbbb2744-a2d9-412a-9938-e0a97aa2620b
@bind PI combine() do Child
	md"""
	``T_{set}`` $(
		Child(Slider(200:10:400, default = 250, show_value = true))
	) 
	
	``T_{span}`` $( 
		Child(Slider(1:1:40, default = 20, show_value = true))
	) ``K_p`` $( 
		Child(Slider(0:0.1:4, default = 2, show_value = true))
	) ``K_i`` $( 
		Child(Slider(0:0.1:4, default = 1.0, show_value = true))
	) 

	
	"""
end

# ╔═╡ fed350c9-aaa1-4a1a-9338-7f244b5e738b
let 
	Tset = PI[1]
	span = PI[2]
	Kp = PI[3]
	Ki = PI[4]
	ts, Ts, Vs, integral = simulation_PI(300.0, 3600; Tset = Tset, 
		span = span, Kp = Kp, Ki = Ki)
	p1 = plot(ts./60.0, Vs, label = :none, color = :black, ylabel = "V")
	p2 = plot(ts./60.0, Ts, color = :black, ylabel = "T (K)", label = "T")
	p2 = plot!([0, ts[end]./60.0], [Tset, Tset], label = "Tset")
	p3 = plot(ts./60.0, integral, color = :black, label = "Integral", 
		xlabel = "Time (min)")
	plot(p1, p2, p3, layout = grid(3,1), size = (700,500))
end

# ╔═╡ f4a47b15-3c2c-4f08-b611-08d57971b04f
md"""
### Advantages
The integral term can clearly eliminate the bias from the P-controller.

### Disadvantages
The I term can introduce osillations around the setpoint. The ``K_i`` value needs to be tuned for the system.
"""

# ╔═╡ 3fba7ef0-a994-4f7e-934a-cb475ceea441
md"""
## Proportional-Integral-Derivative (PID) Control
$(Resource(control3_url, :width => 2800px))

### Description

The derivative of the process error is calculated by determining the slope of the error over time and multiplying this rate of change by the derivative gain ``K_d``. Derivative action predicts system behavior and thus improves settling time and stability of the system.

The full PID controller equation is given by 

```math
e(t) = \frac{r - y(t)}{s}
```

```math
g(t) = K_p e(t)  + K_i \int_0^t e(\tau) d\tau + K_d \frac{de(t)}{dt}
```

where ``K_p`` is proportional gain, ``r`` is the setpoint, ``y(t)`` is the output variable, ``s`` is the span, ``K_i`` is the the integral gain (these terms are the same as for the PI-controller), and ``K_d`` is the derivative gain.

As with the PI-controller, the purpose of the span is to set the output to maximum when the difference between the setpoint and the process variable exceeds the span. Furthermore, in practical application the output needs to be constrained between ``-1`` and ``1``.

```math
f(t) = \begin{cases}
-1 &  g(t) < -1 \\
1 & g(t) > 1 \\
g(t) & else \\
\end{cases}
```

The case statement ensures that output ranges from ``-1`` (here maximum heating) to ``+1`` (here maximum cooling). Finally the output is scaled to ``\pm 12V``, which corresponds to the operating voltages of the TEC element. 

### Implementation
"""

# ╔═╡ 9c1588f1-b49d-470e-b71c-9e70bebafc25
function simulation_PID(T, n; Tset = 280.0, Kp = 1.4, span = 20.0, Ki = 1.0, Kd = 1.0)
	dt = 1.0 
	ts = Float64[]
	Ts = Float64[]
	Vs = Float64[]
	integral = Float64[]
	derivative = Float64[]
	push!(Ts, T)
	push!(Ts, T)
	push!(Vs, 12.0)
	push!(Vs, 12.0)
	push!(ts, 0.0)
	push!(ts, 1.0)
	push!(integral, 0.0)
	push!(integral, 0.0)
	push!(derivative, 0.0)
	push!(derivative, 0.0)
	for i = 1:n
		theT = Ts[end]
		theV = Vs[end]
		error = (Ts[end] - Tset)/span
		error_past = (Ts[end-1] - Tset)/span
		
		theIntegral = integral[end] + error*Ki*dt
		theDerivative = (error_past - error)/dt
		MV = Kp*error + Ki*theIntegral + Kd*theDerivative
		
		if MV < -1
			MV = -1
		elseif MV > 1
			MV = 1
		end

		theV = 12.0*MV
		newT = update_system(theT, theV, dt)
		push!(Ts, newT)
		push!(Vs, theV)
		push!(integral, theIntegral)
		push!(derivative, theDerivative)
		push!(ts, ts[end] + dt)
	end

	return ts, Ts, Vs, integral, derivative
end

# ╔═╡ 5fe29e57-fe90-4ae8-9beb-3662c51d9c16
@bind PID combine() do Child
	md"""
	``T_{set}`` $(
		Child(Slider(200:10:400, default = 250, show_value = true))
	) ``T_{span}`` $( 
		Child(Slider(1:1:40, default = 20, show_value = true))
	) 
	
	``K_p`` $( 
		Child(Slider(0:0.1:4, default = 2, show_value = true))
	) ``K_i`` $( 
		Child(Slider(0:0.1:4, default = 1.0, show_value = true))
	) ``K_d`` $( 
		Child(Slider(0:1.0:40, default = 0.0, show_value = true))
	) 

	
	"""
end

# ╔═╡ 537f07ab-f141-4a38-9dcd-be8de62fc3e5
let 
	Tset = PID[1]
	span = PID[2]
	Kp = PID[3]
	Ki = PID[4]
	Kd = PID[5]
	ts, Ts, Vs, integral, derivative = simulation_PID(300.0, 3600; Tset = Tset, 
		span = span, Kp = Kp, Ki = Ki, Kd = Kd)
	p1 = plot(ts./60.0, Vs, label = :none, color = :black, ylabel = "V")
	p2 = plot(ts./60.0, Ts, color = :black, ylabel = "T (K)", label = "T")
	p2 = plot!([0, ts[end]./60.0], [Tset, Tset], label = "Tset")
	p3 = plot(ts./60.0, integral, color = :black, label = "Integral")
	p4 = plot(ts./60.0, derivative, color = :black, label = "Derivative", 
	 	xlabel = "Time (min)")
	
	plot(p1, p2, p3, p4, layout = grid(4,1), size = (700,500))
end

# ╔═╡ 0a42d713-ecc8-41cb-90f0-30f3c598c1af
md"""
### Advantages
The derivative term can improve convergence rates to the setpoint.

### Disadvantages
Noise in derivative calculations can lead to error amplification. In many applications PI control is preferred to PID control.
"""

# ╔═╡ 4267d65e-1cfd-4fc1-ba61-3f15ba5885c6
md"""
## Trajactory Control
$(Resource(control4_url, :width => 2800px))

### Description

The standard PID controller operates with a fixed set point ``r``. However, the set point may itself be a function of time, ``r(t)``, thus forcing a temperature trajectory, rather than a fixed-point control. 

It is possible to try to force the system along the trajectory by updating ``r(t)`` at each time step. The example below tries to force the system with a linear cooling rate of 1 K min⁻¹.

"""

# ╔═╡ 06007dd6-94dc-4ed6-bf53-18a44db17c06
function simulation_PID_trajectory(
	T, n; Kp = 1.4, span = 20.0, Ki = 1.0, Kd = 1.0, cr = 1.0)

	dt = 1.0 
	ts = Float64[]
	Ts = Float64[]
	Vs = Float64[]
	Tset = Float64[]
	integral = Float64[]
	derivative = Float64[]
	push!(Ts, T)
	push!(Ts, T)
	push!(Vs, 12.0)
	push!(Vs, 12.0)
	push!(ts, 0.0)
	push!(ts, 1.0)
	push!(integral, 0.0)
	push!(integral, 0.0)
	push!(derivative, 0.0)
	push!(derivative, 0.0)
	push!(Tset, 300.0)
	push!(Tset, 300.0)
	for i = 1:n
		theT = Ts[end]
		theV = Vs[end]
		error = (Ts[end] - Tset[end])/span
		error_past = (Ts[end-1] - Tset[end-1])/span
		
		theIntegral = integral[end] + error*Ki*dt
		theDerivative = (error_past - error)/dt
		MV = Kp*error + Ki*theIntegral + Kd*theDerivative
		
		if MV < -1
			MV = -1
		elseif MV > 1
			MV = 1
		end

		theV = 12.0*MV
		newT = update_system(theT, theV, dt)
		newTset = Tset[end] - cr/60.0*dt
		push!(Ts, newT)
		push!(Vs, theV)
		push!(integral, theIntegral)
		push!(derivative, theDerivative)
		push!(Tset, newTset)
		push!(ts, ts[end] + dt)
	end

	return ts, Ts, Vs, integral, derivative, Tset
end

# ╔═╡ 6fc7e520-e85e-4b96-ae6a-1d1312e9e6c8
@bind PID_traj combine() do Child
	md"""
	``c_{r}\;[K\;min^{-1}]`` $( 
		Child(Slider(0:0.1:2, default = 1, show_value = true))
	) 
	
	``T_{span}`` $( 
		Child(Slider(1:1:40, default = 20, show_value = true))
	) 
	
	``K_p`` $( 
		Child(Slider(0:0.1:4, default = 2, show_value = true))
	) ``K_i`` $( 
		Child(Slider(0:0.1:4, default = 1.0, show_value = true))
	) ``K_d`` $( 
		Child(Slider(0:1.0:40, default = 0.0, show_value = true))
	) 

	
	"""
end

# ╔═╡ e2698160-39e6-4a76-ac17-6ccd56fb151f
let 
	cr = PID_traj[1]
	span = PID_traj[2]
	Kp = PID_traj[3]
	Ki = PID_traj[4]
	Kd = PID_traj[5]
	
	ts, Ts, Vs, integral, derivative, Tset = simulation_PID_trajectory(
		300.0, 3600; span = span, Kp = Kp, Ki = Ki, Kd = Kd, cr = cr)
	p1 = plot(ts./60.0, Vs, label = :none, color = :black, ylabel = "V")
	p2 = plot(ts./60.0, Ts, color = :black, ylabel = "T (K)", label = "T")
	p2 = plot!(ts./60.0, Tset, label = "Tset")
	p3 = plot(ts./60.0, integral, color = :black, label = "Integral")
	p4 = plot(ts./60.0, derivative, color = :black, label = "Derivative", 
	 	xlabel = "Time (min)")
	
	plot(p1, p2, p3, p4, layout = grid(4,1), size = (700,500))
end

# ╔═╡ 2ae4f6e6-a1f0-4618-b0a2-ca3e4c73ee7f
md"""
For this particular system, forcing a linear profile is easy. However, oscillations may occur for ``K_i`` and ``K_d`` settings, particularly for small span values. Tuned PID parameters (see more below) may be different from the constant set-point case. 
"""

# ╔═╡ d3780235-46e6-472c-9f46-18ac5cf8ee40
md"""
## Hardware Controllers

Many hardware controllers have on-board PID capabilities. The picture shows a TE Technology TC-36-25-RS232 device, which is a bi-polar proportional-integral-derivative temperature controller that can modulate power input to a thermolectric device. It communicates through an RS 232 port. 

$(Resource(control5_url, :width => 2800px))

The controller stores the ``T_{set}``, ``T_{span}``, ``K_p``, ``K_i``, and ``K_d`` values in the unit's memory. The values can be changed through serial commands passed to the controller. The unit takes fixed voltage power as input provides regulated output (``\pm 100\%`` of the input) to a TEC module. The advantage of hardware controllers is that they can operate offline as part of an instrument. 
"""

# ╔═╡ 945aba6b-ad67-44e7-ae6b-50cff0ac014a
md"""
## PID Tuning

PID parameters can be tuned either on the physical system, or by simulating the system as we did above. Manual tuning is common. One tuning method is to first set ``K_{i}`` and ``K_{d}`` values to zero. Next, increase the ``K_{p}`` value until oscillations are observed. Then set ``K_{p}`` to approximately half that value. Next, increase ``K_{i}`` until any offset is corrected in sufficient time for the process, but not until too great a value causes instability. Finally, increase ``K_{d}``, if required, until the loop is acceptably quick to reach its reference after a load disturbance. 
"""

# ╔═╡ 57d40abd-2a35-468a-9334-f63d962ba604
md"""
# Some Key Concepts in Systems Control 

## Mutli-Input Control System

SISO - Single Input, Single Output

These systems use data/input from one sensor to control one output. These are the simplest to design since they correspond one sensor to one actuator. For example, temperature (TC) is used to control the valve state of v1 through a PID controller.

SIMO - Single Input, Multiple Output

These systems use data/input from one sensor to control multiple outputs. For example, temperature (TC) is used to control the valve state of v1 and v2 through PID controllers.

MISO - Multiple Input, Single Output

These systems use data/input from multiple sensors to control one ouput. For example, a cascade controller can be considered MISO. Temperature (TC) is used in a PID controller (#1) to determine a flow rate set point i.e. FCset. With the FCset and FC controller, they are used to control the valve state of v1 through a PID controller (#2).

MIMO - Multiple Input, Multiple Output

These systems use data/input from multiple sensors to control multiple outputs. These are usually the hardest to design since multiple sensor data is integrated to coordinate multiple actuators. For example, flow rate (FC) and temperature (TC) are used to control multiple valves (v1, v2, and v3). Often, MIMO systems are not PID controllers but rather designed for a specific situation.

## Classification of Systems

- Continuous Linear Time-Invariant Systems (LTI Systems)

```math
\begin{eqnarray}
 & \quad \vec{q}'(t) &= \mathbf{A}\vec{q}(t) + \mathbf{B} \vec{f}(t) \\
 & \quad \vec{y}(t) &=  \mathbf{C}\vec{q}(t) + \mathbf{D} \vec{f}(t) \\
\end{eqnarray}
```

- Continuous Linear Time-Variant (LTV Systems)

```math
\begin{eqnarray}
 & \quad \vec{q}'(t) &= \mathbf{A}(t)\vec{q}(t) + \mathbf{B}(t) \vec{f}(t) \\
 & \quad \vec{y}(t) &=  \mathbf{C}(t)\vec{q}(t) + \mathbf{D}(t) \vec{f}(t) \\
\end{eqnarray}
```

- Non-linear systems

## Controllability and Observability
Controllability and observability represent two major concepts of modern control
system theory. They can be roughly defined as follows.
"""

# ╔═╡ 0d40bd72-6be9-4127-9c63-6041473a3ccc
Markdown.MD(
	Markdown.Admonition("warning", "Key Concepts", [md"
**Controllability:** In order to be able to do whatever we want with the given
dynamic system under control input, the system must be controllable.

**Observability:** In order to see what is going on inside the system under observation, the system must be observable.
	"]))

# ╔═╡ b9f569a7-3854-40c8-8e64-de13f5eb19cf
md"""
For linear systems, controllability and observability are determined by the matrices ``\mathbf{A}`` and ``\mathbf{B}``.

The controllability matrix for LTI systems is given by

```math
 \mathbf{R}={\begin{bmatrix}\mathbf{B}&\mathbf{AB}&\mathbf{A^{{2}}B}&...&\mathbf{A}^{{n-1}}\mathbf{B}\end{bmatrix}}
```

The system is controllable if the controllability matrix has full row rank. The solution is more difficult for LTV systems, but remains closely related to the matrices ``\mathbf{A}`` and ``\mathbf{B}``
"""

# ╔═╡ 4af334c0-c8bd-416f-956e-6c702a6b5f01
md"""

## Optimal Control

Optimal control theory is a branch of control theory that deals with finding a control for a dynamical system over a period of time such that an objective function is optimized. For example, the dynamical system might be a spacecraft with controls corresponding to rocket thrusters, and the objective might be to reach the Moon with minimum fuel expenditure

"""

# ╔═╡ Cell order:
# ╟─cc5fa8b8-4da3-11ee-28fd-7dacab544db3
# ╟─aabae4cb-3912-408d-bacf-267f2440c047
# ╟─9db0e27f-32f5-4b03-9d70-5dea9f079ecb
# ╟─ba749db8-c015-4137-ac7e-92fc0bfbe7ce
# ╟─a25d1d5d-abc8-4ec4-80b2-ecdaea755a1c
# ╠═be1a604c-6cf6-405f-9fa0-86e3b454315e
# ╟─759f6f53-866c-4cc8-9bb5-41be3b43f1dc
# ╠═bdc5a009-c6c4-434d-a619-63081c6af32d
# ╠═f4c5f8dd-a66f-4b1a-8a6f-3e4ad9615999
# ╟─b8fe4424-9eac-4011-85ba-f1bf5e72c718
# ╟─1b1dbf06-c67b-4447-9049-4ab1bf82ee8a
# ╠═cf63fd31-4c00-46e5-8648-713b59080eba
# ╟─451a1910-4ce9-4e15-bca3-c516e20e5bd5
# ╠═c42f51a3-a02d-4c86-beee-7feedd378ad0
# ╠═538e7d53-5cc3-4447-b1c7-3d02a0946f09
# ╟─6369c1ac-e52c-4d81-869d-77cbd483fd8b
# ╟─03a0fc10-9875-4248-8711-8a3887db30a5
# ╟─8aa24500-d236-4e24-bf06-4db94e612a50
# ╠═54a7e106-4537-4857-a9e7-fdc3e9bf4ba1
# ╠═f0ecf625-53f3-4595-8764-ab0a95c5b6e3
# ╟─17c3d4a5-048a-44c1-a5f0-5278a2322f60
# ╠═25b351fa-69e6-482e-9543-12e1dde9e438
# ╠═1163c920-d2ca-42e7-9b9e-1e19d494bb46
# ╠═43c38357-4d59-457d-9b90-94ae1f267212
# ╟─7d030c04-ac3b-4832-a3e8-6053fc9c14e7
# ╟─1ae8b353-df38-430f-b220-401adddac1d1
# ╠═f85d4b96-e5d8-4d5c-aaed-df30c8c06084
# ╟─0c6f2477-1252-4485-9596-2748a2508e14
# ╠═77c08e89-b5ca-45d3-8916-3ac8cf9453cb
# ╟─6647a86f-7f6d-43d7-892e-724d3a191ca9
# ╠═b40aefaa-4600-42f7-ab1d-d92b54fb6401
# ╟─a20f135e-f061-4287-9e11-3ffad5cc64b8
# ╠═177d0f02-9045-49d5-84cd-bc7b6ebb58fc
# ╟─2a2f1d95-60d1-41b5-9d94-e69e7c5bba78
# ╟─9e35dc47-ff08-4254-96dd-fe27c0c7366d
# ╠═6e2bc9b4-6b61-49c8-ae72-0bc94e72543a
# ╟─bbbb2744-a2d9-412a-9938-e0a97aa2620b
# ╟─fed350c9-aaa1-4a1a-9338-7f244b5e738b
# ╟─f4a47b15-3c2c-4f08-b611-08d57971b04f
# ╟─3fba7ef0-a994-4f7e-934a-cb475ceea441
# ╠═9c1588f1-b49d-470e-b71c-9e70bebafc25
# ╟─5fe29e57-fe90-4ae8-9beb-3662c51d9c16
# ╟─537f07ab-f141-4a38-9dcd-be8de62fc3e5
# ╟─0a42d713-ecc8-41cb-90f0-30f3c598c1af
# ╟─4267d65e-1cfd-4fc1-ba61-3f15ba5885c6
# ╟─06007dd6-94dc-4ed6-bf53-18a44db17c06
# ╟─6fc7e520-e85e-4b96-ae6a-1d1312e9e6c8
# ╟─e2698160-39e6-4a76-ac17-6ccd56fb151f
# ╟─2ae4f6e6-a1f0-4618-b0a2-ca3e4c73ee7f
# ╟─d3780235-46e6-472c-9f46-18ac5cf8ee40
# ╟─945aba6b-ad67-44e7-ae6b-50cff0ac014a
# ╟─57d40abd-2a35-468a-9334-f63d962ba604
# ╟─0d40bd72-6be9-4127-9c63-6041473a3ccc
# ╟─b9f569a7-3854-40c8-8e64-de13f5eb19cf
# ╟─4af334c0-c8bd-416f-956e-6c702a6b5f01
