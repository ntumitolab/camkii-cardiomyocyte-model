# # CaMK parameter fitting
# Fitting CaMKII parameters to experimental decay rates.
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using Optimization
using OptimizationOptimJL
using OptimizationLBFGSB
using ADTypes
using ForwardDiff
using DiffEqCallbacks
using CurveFit
using Plots
using LinearAlgebra
using Model
using Model: second

# ## Experimental data
# CaMKII activity decay times after 1Hz pacing for 15.0, 30.0, 60.0, 90.0 seconds.
pacing_durations = [15.0, 30.0, 60.0, 90.0]
experimental_taus = [16.48, 16.73, 17.65, 18.08]

# ## Simplified calcium transient
# We use a fitted model for calcium transients to speed up CaMKII parameter fitting.
stimstart = 30.0second
stimend = 130.0second
tend = 205second
alg = KenCarp47()
@time "Building ODE system" sys = Model.DEFAULT_SYS
@time "Building ODE problem" prob = ODEProblem(sys, [sys.k_P1_P2=>0], tend)

@unpack Istim = sys
pace15 = Model.build_stim_callbacks(Istim, stimstart + 15 * second; period=1second, starttime=stimstart)
pace30 = Model.build_stim_callbacks(Istim, stimstart + 30 * second; period=1second, starttime=stimstart)
pace60 = Model.build_stim_callbacks(Istim, stimstart + 60 * second; period=1second, starttime=stimstart)
pace90 = Model.build_stim_callbacks(Istim, stimstart + 90 * second; period=1second, starttime=stimstart);

# ## Test run
@time "Solve problems" sols = map([pace15, pace30, pace60, pace90]) do cb
    solve(prob, alg, callback = cb)
end

plot()
for (sol, t) in zip(sols, pacing_durations)
    plot!(sol, idxs=CaMKAct, label="Paced at $(t) seconds")
end
plot!(xlabel="Time (ms)", ylabel="CaMKII Activity", title="Simulated calcium transient", legend=:topright)

# ## Loss function
# Changing `kdeph_CaMK` and `k_P1_P2` to fit decay rate in the experiments.
@unpack kphos_CaMK, kdeph_CaMK, kb_CaMKP, CaMKAct = sys
kphos_dephos_ratio = prob.ps[kphos_CaMK] / prob.ps[kdeph_CaMK]
data = (
    prob=prob,
    cbs=[pace15, pace30, pace60, pace90],
    experimental_taus=experimental_taus,
    pacing_durations=pacing_durations,
    kphos_dephos_ratio=kphos_dephos_ratio,
    stimstart = stimstart,
);

function loss(theta, data)
    @unpack prob, cbs, experimental_taus, pacing_durations, kphos_dephos_ratio, stimstart = data
    @unpack kdeph_CaMK, kphos_CaMK, kb_CaMKP, CaMKAct = prob.f.sys
    dephos_rate = exp10(theta[1])
    kb2_rate = exp10(theta[2])
    ## Parallel ensemble simulation
    function prob_func(prob, i, repeat)
        remake(prob, p=[
            kdeph_CaMK => dephos_rate,
            kphos_CaMK => kphos_dephos_ratio * dephos_rate,
            kb_CaMKP => kb2_rate
            ],
            callback=cbs[i])
    end

    ## Calculate loss in the output function
    function output_func(sol, i)
        SciMLBase.successful_retcode(sol) || return (Inf, false)
        stimend = stimstart + pacing_durations[i] * second
        ts = collect(range(0.0, stop=50.0, step=5.0))
        ysim = sol(stimend .+ ts .* second, idxs=CaMKAct).u
        fit = solve(CurveFitProblem(ts, ysim), ExpSumFitAlgorithm(n=1, withconst=true))
        tau = inv(-fit.u.λ[])
        tauexpected = experimental_taus[i]
        error = (tau - tauexpected)^2
        return (error, false)
    end

    ensemble_prob = EnsembleProblem(prob; prob_func, output_func)
    sim = solve(ensemble_prob, KenCarp47(); trajectories=length(pacing_durations), maxiters=10000, verbose=false)
    return sum(sim)
end

# Test the loss function.
theta = [log10(prob.ps[kdeph_CaMK]), log10(prob.ps[kb_CaMKP])]
@time loss(theta, data)

# ## Optimization
optf = OptimizationFunction(loss)
optprob = OptimizationProblem(optf, theta, data, lb=[-1, -1] + theta, ub=[1, 1] + theta)
optalg = Optim.SAMIN()
@time fitted_dephos = solve(optprob, optalg; maxiters=4000, maxtime=3600)
@show fitted_dephos.objective
fitted_dephos.stats

#---
params = Dict(kdeph_CaMK => exp10(fitted_dephos.u[1]), kphos_CaMK => data.kphos_dephos_ratio * exp10(fitted_dephos.u[1]), kb_CaMKP => exp10(fitted_dephos.u[2]))

println("Fitted parameters:")
println("Dephosphorylation time: " , 1e-3 / params[kdeph_CaMK], " seconds.")
println("Autophosphorylation rate of CaMKB: " , params[kphos_CaMK] * 1000, " Hz")
println("Dephosphorylation time of CaMKP: " , 1e-3 / params[kb_CaMKP], " seconds.")

# ## Test fitted parameters
newprob = remake(prob, p=[kdeph_CaMK => params[kdeph_CaMK], kphos_CaMK => params[kphos_CaMK], kb_CaMKP => params[kb_CaMKP]])

sols = map([pace15, pace30, pace60, pace90]) do cb
    solve(newprob, alg, callback = cb)
end;

plot()
for (sol, t) in zip(sols, [15.0, 30.0, 60.0, 90.0])
    plot!(sol, idxs=CaMKAct, label="Paced at $(t) seconds")
end
plot!(xlabel="Time (ms)", ylabel="CaMKII Activity", title="Simulated calcium transient", legend=:topright)

# ## Decay rates
# Fit data from simulations against an exponential decay model.
# Record 50 seconds after pacing ends.
ts = collect(range(0.0, stop=50.0, step=5.0)) ## in seconds
stimstart = 30.0second
ysim_15 = sols[1](stimstart+15second:5second:stimstart+15second+50second ; idxs=sys.CaMKAct * 100).u
ysim_30 = sols[2](stimstart+30second:5second:stimstart+30second+50second ; idxs=sys.CaMKAct * 100).u
ysim_60 = sols[3](stimstart+60second:5second:stimstart+60second+50second ; idxs=sys.CaMKAct * 100).u
ysim_90 = sols[4](stimstart+90second:5second:stimstart+90second+50second ; idxs=sys.CaMKAct * 100).u

fit_sim_15 = solve(CurveFitProblem(ts, ysim_15), ExpSumFitAlgorithm(n=1, withconst=true))
fit_sim_30 = solve(CurveFitProblem(ts, ysim_30), ExpSumFitAlgorithm(n=1, withconst=true))
fit_sim_60 = solve(CurveFitProblem(ts, ysim_60), ExpSumFitAlgorithm(n=1, withconst=true))
fit_sim_90 = solve(CurveFitProblem(ts, ysim_90), ExpSumFitAlgorithm(n=1, withconst=true))

# ### Fitting results (simulations)
p1s = plot(ts, ysim_15, label="Sim 15 sec")
plot!(p1s, ts, predict(fit_sim_15), label="Fit", linestyle=:dash)
p2s = plot(ts, ysim_30, label="Sim 30 sec")
plot!(p2s, ts, predict(fit_sim_30), label="Fit", linestyle=:dash)
p3s = plot(ts, ysim_60, label="Sim 60 sec")
plot!(p3s, ts, predict(fit_sim_60), label="Fit", linestyle=:dash)
p4s = plot(ts, ysim_90, label="Sim 90 sec")
plot!(p4s, ts, predict(fit_sim_90), label="Fit", linestyle=:dash)
plot(p1s, p2s, p3s, p4s, layout=(2,2), xlabel="Time (s)", ylabel="CaMKII activity (%)")

# ### Decay time scales (tau)
tau_sim_15 = inv(-fit_sim_15.u.λ[])
tau_sim_30 = inv(-fit_sim_30.u.λ[])
tau_sim_60 = inv(-fit_sim_60.u.λ[])
tau_sim_90 = inv(-fit_sim_90.u.λ[])

println("The time scales for simulations: ")
for (tau, dur) in zip((tau_sim_15, tau_sim_30, tau_sim_60, tau_sim_90), (15, 30, 60, 90))
    println("$dur sec pacing is $(round(tau; digits=2)) seconds.")
end
