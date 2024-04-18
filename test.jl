using OrdinaryDiffEq
using DiffEqCallbacks
using Plots

function f(u, p, t)
    a, d, e = p
    return u * (a - d) + e
end

u0 = 1.0
p = [0.0, 1.0, 1.0]
tspan = (0.0, 100.0)
prob = ODEProblem(f, u0, tspan, p)

function make_events(start_time, end_time, pulse_width, period, strength=2.0)
    function go_up!(integrator)
        integrator.p[1] = strength
        integrator.p[3] = 0.0
    end

    function go_down!(integrator)
        integrator.p[1] = 0.0
        integrator.p[3] = 1.0
    end

    upevents = PresetTimeCallback(start_time:period:end_time, go_up!)
    downevents = PresetTimeCallback(start_time+pulse_width:period:end_time, go_down!)
    return CallbackSet(upevents, downevents)
end

alg = Tsit5()
events = make_events(tspan[begin], tspan[end], 1.0, 10.0)

sol = solve(prob, alg, callback=events)

plot(sol)
