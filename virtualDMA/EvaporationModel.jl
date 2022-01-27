module EvaporationModel

using DifferentialEquations
using ParameterizedFunctions
using Interpolations

ODEs = @ode_def_bare begin
	dDp = -G/Dp
end G

function evapmodel(Dd, G)
	function condition(u, t, integrator)
		out = (u[1] < 5e-9) ? 0.0 : 1.0
	end
	affect!(integrator) = terminate!(integrator)
	cb = ContinuousCallback(condition, affect!)

	prob = ODEProblem(ODEs, [Dd], (0.0, 3600.0), [G])
	sol = solve(prob, RK4(), abstol = 1e-20, reltol = 1e-20, callback = cb)
end

function getEFfunction(Ds, G)
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

Ds = exp10.(range(-8, stop = -6, length = 240))

@time f = EvaporationModel.getEFfunction(Ds, 1e-18)
f(20)
#@time f = evaluate_model(Ds, 1e-17)
#@time f = evaluate_model(Ds, 1e-16)

