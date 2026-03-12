# # Parameter fitting
# Fitting CaMKII parameters to experimental decay rates
using ModelingToolkit
using OrdinaryDiffEq
using CurveFit
using Optimization
using OptimizationBBO
using CaMKIIModel
using CaMKIIModel: second

# Setup simulation
tend = 210second
@time "Build system" @mtkcompile sys = build_neonatal_ecc_sys()
@time "Build problem" prob = ODEProblem(sys, [], tend)
callback15 = build_stim_callbacks(Istim, stimstart + 15second; period=1second, starttime=stimstart)
sol15 = solve(prob, alg; callback=callback15)
callback30 = build_stim_callbacks(Istim, stimstart + 30second; period=1second, starttime=stimstart)
sol30 = solve(prob, alg; callback=callback30)
callback60 = build_stim_callbacks(Istim, stimstart + 60second; period=1second, starttime=stimstart)
sol60 = solve(prob, alg; callback=callback60)
callback90 = build_stim_callbacks(Istim, stimstart + 90second; period=1second, starttime=stimstart)
sol90 = solve(prob, alg; callback=callback90)

# Decay times for 15sec, 30 sec, 60 sec, and 90 sec pacing
# Fitted in `Pacing durations`

# Loss function
function loss(p; alg = FBDF(), stimstart = 30second)
    const experimental_durations = [16.48, 16.73, 17.65, 18.08]
    simulation_durations = map((callback15, callback30, callback60, callback90)) do callback
        prob = remake(prob, p=p) ## TODO: remake with p only for the parameters we want to fit
        sol = solve(prob, alg; callback=callback)
        ysim = sol(stimstart+15second:5second:stimstart+15second+50second; idxs=sys.CaMKAct * 100).u
        ts = collect(range(0.0, stop=50.0, step=5.0))
        fit_sim = solve(CurveFitProblem(ts, ysim), ExpSumFitAlgorithm(n=1, withconst=true))
        tau_sim = inv(-fit_sim.u.λ[])
    end
    return sum(abs2, experimental_durations .- simulation_durations)
end
