# # Pacing frequency
using ModelingToolkit
using OrdinaryDiffEq
using Plots
using CSV
using DataFrames
using CaMKIIModel
Plots.default(lw=1.5)

# ## Setup the ODE system
# Electrical stimulation starts at `t`=100 seconds and ends at `t`=300 seconds
sys = build_neonatal_ecc_sys(simplify=true)
tend = 500.0
prob = ODEProblem(sys, [], tend)
stimstart = 100.0
stimend = 300.0
@unpack Istim = sys
alg = FBDF()

# ## 1Hz
callback = build_stim_callbacks(Istim, stimend; period=1, starttime=stimstart)
@time sol = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

#---
plot(sol, idxs=sys.vm*1000, title="Action potential")

#---
plot(sol, idxs=sys.vm*1000, title="Action potential", tspan=(299, 300))

#---
plot(sol, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299, 300), title="Calcium transient")

#---
plot(sol, idxs=sys.CaMKAct, title="Active CaMKII")

# ## 2Hz
callback = build_stim_callbacks(Istim, stimend; period=1/2, starttime=stimstart)
@time sol2 = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

#---
plot(sol2, idxs=sys.vm*1000, title="Action potential")

#---
plot(sol2, idxs=sys.vm*1000, title="Action potential", tspan=(299, 300))

#---
plot(sol2, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299, 300), title="Calcium transient")

#---
plot(sol2, idxs=sys.CaMKAct, title="Active CaMKII")

# ## 3Hz
callback = build_stim_callbacks(Istim, stimend; period=1/3, starttime=stimstart)
@time sol3 = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

#---
plot(sol3, idxs=sys.vm*1000, title="Action potential")

#---
plot(sol3, idxs=sys.vm*1000, title="Action potential", tspan=(299, 300))

#---
plot(sol3, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299, 300), title="Calcium transient")

#---
plot(sol3, idxs=sys.CaMKAct, title="Active CaMKII")

# ## Comparing 1-3 Hz
# Action potential
plot(sol, idxs=sys.vm*1000, title="Action potential", lab="1Hz")
plot!(sol2, idxs=sys.vm*1000, lab="2Hz")
plot!(sol3, idxs=sys.vm*1000, lab="3Hz", tspan=(299, 300), xlabel="Time (sec.)", ylabel="Voltage (mV)")

# CaMKII activity
plot(sol, idxs=sys.CaMKAct, title="Active CaMKII", lab="1Hz")
plot!(sol2, idxs=sys.CaMKAct, lab="2Hz")
plot!(sol3, idxs=sys.CaMKAct, lab="3Hz", ylim=(0.0, 0.8))

# ## Data fitting
# ### Pacing duration and Ca transient
ts = durationdf[!, "Time(sec)"]
durationdf = CSV.read("data/CaMKAR-duration.csv", DataFrame)
fifteen = durationdf[!, "1Hz 15sec (Mean)"]
fifteen_error = durationdf[!, "1Hz 15sec (SD)"] ./ sqrt.(durationdf[!, "1Hz 15sec (N)"])
thirty = durationdf[!, "1Hz 30sec (Mean)"]
thirty_error = durationdf[!, "1Hz 30sec (SD)"] ./ sqrt.(durationdf[!, "1Hz 30sec (N)"])
sixty = durationdf[!, "1Hz 60sec (Mean)"]
sixty_error = durationdf[!, "1Hz 60sec (SD)"] ./ sqrt.(durationdf[!, "1Hz 60sec (N)"])
ninety = durationdf[!, "1Hz 90sec (Mean)"]
ninety_error = durationdf[!, "1Hz 90sec (SD)"] ./ sqrt.(durationdf[!, "1Hz 90sec (N)"])

plot(ts, fifteen, yerr=fifteen_error, lab="15 sec", color = :blue, markerstrokecolor=:blue)
plot!(ts, thirty, yerr=thirty_error, lab="30 sec", color = :red, markerstrokecolor=:red)
plot!(ts, sixty, yerr=sixty_error, lab="60 sec", color = :orange, markerstrokecolor=:orange)
plot!(ts, ninety, yerr=ninety_error, lab="90 sec", color = :green, markerstrokecolor=:green)
plot!(title="Pacing duration", xlabel="Time (sec.)", ylabel="Calcium (AU)")

# Resolution is reduced to 1Hz
callback15 = build_stim_callbacks(Istim, stimstart+15; period=1, starttime=stimstart)
sol15 = solve(prob, alg; callback=callback15, abstol=1e-6, reltol=1e-6, saveat=1.0)
callback30 = build_stim_callbacks(Istim, stimstart+30; period=1, starttime=stimstart)
sol30 = solve(prob, alg; callback=callback30, abstol=1e-6, reltol=1e-6, saveat=1.0)
callback60 = build_stim_callbacks(Istim, stimstart+60; period=1, starttime=stimstart)
sol60 = solve(prob, alg; callback=callback60, abstol=1e-6, reltol=1e-6, saveat=1.0)
callback90 = build_stim_callbacks(Istim, stimstart+90; period=1, starttime=stimstart)
sol90 = solve(prob, alg; callback=callback90, abstol=1e-6, reltol=1e-6, saveat=1.0)

idx = sys.Cai_mean*1000
plot(sol15, idxs=idx, tspan=(100, 300), lab="15 sec")
plot!(sol30, idxs=idx, tspan=(100, 300), lab="30 sec")
plot!(sol60, idxs=idx, tspan=(100, 300), lab="60 sec")
plot!(sol90, idxs=idx, tspan=(100, 300), lab="90 sec")
plot!(title="Pacing duration", xlabel="Time (sec.)", ylabel="Calcium (μM)")
