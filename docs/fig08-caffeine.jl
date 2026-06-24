# # Caffeine Effects
# Caffeine increase RyR opening sensitivity to luminal and subspace calcium.
# In this model, we decrease the mid saturation sub-SR calcium concentration for the opening rate.
using CSV
using DataFrames
using DiffEqCallbacks
using DifferentialEquations
using ModelingToolkit
using OrdinaryDiffEqSDIRK
using Plots
using SteadyStateDiffEq
using Model
using Model: second
Plots.default(lw=1.5)

# ## Experimental data of CaMKII activity with caffeine treatment.
chemicaldf = CSV.read(joinpath(@__DIR__, "data/CaMKAR-chemical-normalized.csv"), DataFrame)
ts = 0:5:205
ctl = chemicaldf[!, "Ctrl Mean"]
ctl_error = chemicaldf[!, "Ctrl SD"] ./ sqrt.(chemicaldf[!, "Ctrl N"])

caf = chemicaldf[!, "caffeine 20mM Mean"]
caf_error = chemicaldf[!, "caffeine 20mM SD"] ./ sqrt.(chemicaldf[!, "caffeine 20mM N"])

fig8a = plot(ts, ctl, yerr=ctl_error, lab="Ctl", color=:blue, markerstrokecolor=:blue)
plot!(fig8a, ts, caf, yerr=caf_error, lab="Caf 20mM", color=:red, markerstrokecolor=:red)
plot!(fig8a, xlabel="Time (s)", ylabel="CaMKAR (R/R0)", title= "A", titlelocation=:left)

# ## Caffeine and electrophysiology
# - Adding caffeine in the beginning of the simulation.
# - Adding caffeine and nifedipine in the beginning of the simulation (nifedipine blocks 90% of LCCs).
@time "Build system" sys = Model.DEFAULT_SYS
tend = 205second
stimstart = 30second
stimend = 120second
alg = FBDF()
@unpack Istim = sys
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
@time prob = ODEProblem(sys, [], tend)
@time prob_caf = ODEProblem(sys, [sys.RyRsensitivity => 10], tend)
gCaL = prob.ps[sys.GCaL]
@time prob_nif_caf = ODEProblem(sys, [sys.RyRsensitivity => 10, sys.GCaL => 0.1 * gCaL], tend)

ssalg = DynamicSS(alg)
sprob_caf = SteadyStateProblem(prob_caf)
sssol = solve(sprob_caf, ssalg; abstol=1e-10, reltol=1e-10)
@time sol = solve(prob, alg; callback)
@time sol_caf = solve(remake(prob_caf, u0=sssol.u), alg; callback)
@time sol_nif_caf = solve(remake(prob_nif_caf, u0=sssol.u), alg; callback)

# ### CaMKII activities
i = (sys.t / 1000, sys.CaMKAct)
fig8b = plot(sol, idxs=i,lab="Ctl", color=:blue)
plot!(fig8b, sol_caf, idxs=i, lab="Caf", color=:green)
plot!(fig8b, sol_nif_caf, idxs=i, lab="Caf + Nif", color=:red, ylabel="Active CaMKII fraction", xlabel="Time (s)", title="B", titlelocation=:left)

# ### Action potential
i = (sys.t / 1000, sys.vm)
tspan = (100second, 102second)
fig8c = plot(sol, idxs=i, lab="Ctl", color=:blue; tspan)
plot!(fig8c, sol_caf, idxs=i, lab="Caf", color=:green; tspan)
plot!(fig8c, sol_nif_caf, idxs=i, lab="Caf + Nif", color=:red, tspan=tspan, ylabel="Voltage (mV)", xlabel="Time (s)", title="C", titlelocation=:left)

# ### Calcium transient
i = (sys.t / 1000, sys.Cai_mean * 1000)
fig8d = plot(sol, idxs=i, lab="Ctl", color=:blue; tspan)
plot!(fig8d, sol_caf, idxs=i, lab="Caf", color=:green; tspan)
plot!(fig8d, sol_nif_caf, idxs=i, lab="Caf + Nif", color=:red, tspan=tspan, ylabel="Mean intracellular Ca (nM)", xlabel="Time (s)", title="D", titlelocation=:left)

# ## Comparing with experiments
fig8 = plot(fig8a, fig8b, fig8c, fig8d, layout=(2, 2), size=(800, 600))
savefig("fig8.pdf")
savefig("fig8.png")

# ## Adding caffeine halfway
# at t = 200 second
@time "Build system" sys = Model.DEFAULT_SYS
tend = 500second
@time "Build problem" prob = ODEProblem(sys, [], tend)
stimstart = 100second
stimend = 300second
alg = FBDF()
function add_coffee_affect!(integrator)
    integrator.ps[sys.RyRsensitivity] = 10
end

#---
@unpack Istim = sys
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
# Add caffeine at t = 200 econd
callback_caf = CallbackSet(build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart), PresetTimeCallback(200.0second, add_coffee_affect!));

# ### Single-dose caffeine
@time sol = solve(prob, alg; callback)
@time sol_caf = solve(prob, alg; callback=callback_caf)

i = (sys.t / 1000, sys.vm)
plot(sol, idxs=i, title="Action potential", lab="Ctl", tspan=(198second, 205second))
plot!(sol_caf, idxs=i, lab="Caf", tspan=(198second, 205second), ylabel="Voltage (mV)", xlabel="Time (s)")

#---
i = (sys.t / 1000, sys.Cai_sub_SR * 1000)
plot(sol, idxs=i, title="Calcium transient (During caffeine addition)", lab="Ctl", tspan=(198second, 205second))
plot!(sol_caf, idxs=i, tspan=(198second, 205second), lab="Caf", ylabel="Subspace calcium (nM)", xlabel="Time (s)")

#---
i = (sys.t / 1000, sys.PO1RyR)
plot(sol, idxs=i, title="RyR open (During caffeine addition)", lab="Ctl", tspan=(198second, 205second))
plot!(sol_caf, idxs=i, tspan=(198second, 205second), lab="Caf", ylabel="Open probability", ylims=(0, 1), xlabel="Time (s)")

#---
i = (sys.t / 1000, sys.Cai_sub_SR * 1000)
plot(sol, idxs=i, title="Calcium transient (After caffeine addition)", lab="Ctl", ylabel="Subspace calcium (nM)", tspan=(198second, 205second))
plot!(sol_caf, idxs=i, lab="Caf", xlabel="Time (s)", tspan=(198second, 205second))

#---
i = (sys.t / 1000, sys.CaJSR)
plot(sol, idxs=i, title="SR Calcium (During caffeine addition)", lab="Ctl", ylabel="SR calcium (μM)", tspan=(198second, 205second))
plot!(sol_caf, idxs=i, tspan=(198second, 205second), lab="Caf", ylims=(0, 850), xlabel="Time (s)")

#---
i = (sys.t / 1000, sys.Jrel)
plot(sol, idxs=sys.Jrel, title="Ca flux", lab="Ctl  (Jrel)", tspan=(198second, 205second))
plot!(sol_caf, idxs=sys.Jrel, lab="Caf (Jrel)", tspan=(198second, 205second), ylabel="μM/ms", xlabel="Time (s)")

#---
i = (sys.t / 1000, sys.CaMKAct)
plot(sol, idxs=i, title="Active CaMKII", lab="Ctl")
plot!(sol_caf, idxs=i, lab="Caf", ylabel="Active CaMKII fraction", xlabel="Time (s)")

#---
i = (sys.t / 1000, sys.Cai_sub_SR * 1000)
tspan = (100second, 102second)
plot(sol, idxs=i, title="Calcium transient", lab="Ctl", color=:blue; tspan)
plot!(sol_caf, idxs=i, lab="Caf", color=:green; tspan)
plot!(sol_nif_caf, idxs=i, lab="Caf + Nif", color=:red, ylabel="Subspace calcium (nM)", xlabel="Time (s)"; tspan)

#---
i = (sys.t / 1000, sys.CaJSR)
tspan = (100second, 102second)
plot(sol, idxs=i, title="SR Calcium", lab="Ctl", ylabel="SR calcium (μM)", color=:blue; tspan)
plot!(sol_caf, idxs=i, lab="Caf", color=:green; tspan)
plot!(sol_nif_caf, idxs=i, lab="Caf + Nif", color=:red, ylims=(0, 850), xlabel="Time (s)"; tspan)
