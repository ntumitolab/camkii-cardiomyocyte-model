# # Pacing frequency
using ModelingToolkit
using OrdinaryDiffEq
using Plots
using CSV
using DataFrames
using CaMKIIModel
Plots.default(lw=1.5)

# ## Setup the ODE system
# Electrical stimulation starts at `t`=100 seconds and ends at `t`=300 seconds.
sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true)
tend = 500.0
prob = ODEProblem(sys, [], tend)
stimstart = 100.0
stimend = 300.0
@unpack Istim = sys
alg = FBDF()

# ## 1Hz
callback = build_stim_callbacks(Istim, stimend; period=1, starttime=stimstart)
sol = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6)

#---
plot(sol, idxs=sys.vm*1000, title="Action potential")

#---
plot(sol, idxs=sys.vm*1000, title="Action potential", tspan=(299, 300))

#---
plot(sol, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299, 300), title="Calcium transient")

# Not RyR problem, ICa problem
prob2 = remake(prob, p = [sys.kRyR => 3600])

sol2 = solve(prob2, alg; callback, abstol=1e-6, reltol=1e-6)

plot(sol2, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(105, 107), title="Calcium transient")

# Imbalance in the SR
plot(sol2, idxs=[sys.CaNSR, sys.CaJSR], tspan=(299, 300), title="Calcium")

#---
plot(sol, idxs=sys.CaMKAct, title="Active CaMKII")

# ## 2Hz
callback = build_stim_callbacks(Istim, stimend; period=1/2, starttime=stimstart)
sol2 = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6)

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
sol3 = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6)

#---
plot(sol3, idxs=sys.vm*1000, title="Action potential")

#---
plot(sol3, idxs=sys.vm*1000, title="Action potential", tspan=(299, 300))

#---
plot(sol3, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299, 300), title="Calcium transient")

#---
plot(sol3, idxs=sys.CaMKAct, title="Active CaMKII")

# ## Comparing 1-3 Hz
plot(sol, idxs=sys.vm*1000, title="Action potential", lab="1Hz")
plot!(sol2, idxs=sys.vm*1000, lab="2Hz")
plot!(sol3, idxs=sys.vm*1000, lab="3Hz", tspan=(299, 300), xlabel="Time (sec.)", ylabel="Voltage (mV)")

#---
plot(sol, idxs=sys.CaMKAct, title="CaMKII", lab="1Hz")
plot!(sol2, idxs=sys.CaMKAct, lab="2Hz")
plot!(sol3, idxs=sys.CaMKAct, lab="3Hz", ylim=(0.0, 1.0), xlabel="Time (sec.)", ylabel="Act. fraction (AU)")

# ## Data fitting
### Pacing duration and CaMKII activity
durationdf = CSV.read(joinpath(@__DIR__, "data/CaMKAR-duration.csv"), DataFrame)
ts = durationdf[!, "Time(sec)"]
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
plot!(title="Pacing duration", xlabel="Time (sec.)", ylabel="CaMKII activity (AU)")

# Simulation resolution is reduced to 1Hz.
stimstart=30.0
callback15 = build_stim_callbacks(Istim, stimstart+15; period=1,starttime=stimstart)
sol15 = solve(prob, alg; callback=callback15, abstol=1e-6, reltol=1e-6)
callback30 = build_stim_callbacks(Istim, stimstart+30; period=1, starttime=stimstart)
sol30 = solve(prob, alg; callback=callback30, abstol=1e-6, reltol=1e-6)
callback60 = build_stim_callbacks(Istim, stimstart+60; period=1, starttime=stimstart)
sol60 = solve(prob, alg; callback=callback60, abstol=1e-6, reltol=1e-6)
callback90 = build_stim_callbacks(Istim, stimstart+90; period=1, starttime=stimstart)
sol90 = solve(prob, alg; callback=callback90, abstol=1e-6, reltol=1e-6)
idx = sys.CaMKAct

plot(sol15, idxs=idx, tspan=(0, 205), lab="15 sec", color = :blue)
plot!(sol30, idxs=idx, tspan=(0, 205), lab="30 sec", color = :red)
plot!(sol60, idxs=idx, tspan=(0, 205), lab="60 sec", color = :orange)
plot!(sol90, idxs=idx, tspan=(0, 205), lab="90 sec", color = :green)
plot!(title="Pacing duration", xlabel="Time (sec.)", ylabel="CaMKII activity (AU)")

# ### Pacing frequency and CaMKII activity
freqdf = CSV.read(joinpath(@__DIR__,"data/CaMKAR-freq.csv"), DataFrame)
ts = 0:5:205
onehz = freqdf[!, "1Hz (Mean)"]
onehz_error = freqdf[!, "1Hz (SD)"] ./ sqrt.(freqdf[!, "1Hz (N)"])
twohz = freqdf[!, "2Hz (Mean)"]
twohz_error = freqdf[!, "2Hz (SD)"] ./ sqrt.(freqdf[!, "2Hz (N)"])

plot(ts, onehz, yerr=onehz_error, lab="1 Hz", color = :blue, markerstrokecolor=:blue)
plot!(ts, twohz, yerr=twohz_error, lab="2 Hz", color = :red, markerstrokecolor=:red)
plot!(title="Pacing frequency", xlabel="Time (sec.)", ylabel="CaMKII activity (AU)")

#---
tend = 205.0
prob = ODEProblem(sys, [], tend)
stimstart = 30.0
stimend = 120.0
callback = build_stim_callbacks(Istim, stimend; period=1, starttime=stimstart)
sol1 = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6)

callback2 = build_stim_callbacks(Istim, stimend; period=0.5, starttime=stimstart)
sol2 = solve(prob, alg; callback=callback2, abstol=1e-6, reltol=1e-6)
idx = sys.CaMKAct

plot(sol1, idxs=idx, lab="1 Hz", color = :blue)
plot!(sol2, idxs=idx, lab="2 Hz", color = :red)
plot!(title="Pacing frequency", xlabel="Time (sec.)", ylabel="CaMKII activity (AU)", ylims=(0.0, 0.9))
