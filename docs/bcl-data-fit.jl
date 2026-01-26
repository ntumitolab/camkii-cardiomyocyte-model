# # Pacing data
# Experiments vs simulations.
using ModelingToolkit
using OrdinaryDiffEq, SteadyStateDiffEq, DiffEqCallbacks
using Plots
using CSV
using DataFrames
using LsqFit
using CaMKIIModel
using CaMKIIModel: second
Plots.default(lw=1.5)

#---
@time "Build system" sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true, reduce_camk=true)
tend = 500second
@time "Build problem" prob = ODEProblem(sys, [], tend)
stimstart = 100second
stimend = 300second
@unpack Istim = sys
alg = KenCarp47()

#===
## Pacing duration and CaMKII activity

### Experiments

30 seconds resting + N seconds 1Hz pacing + resting.
===#
