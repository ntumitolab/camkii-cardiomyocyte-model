using ProgressLogging
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM
Plots.default(lw=1.5)

# ## Setup the ODE system
sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true)
tend = 300.0
prob = ODEProblem(sys, [], tend)
@unpack Istim = sys
alg = FBDF()

callback = build_stim_callbacks(Istim, tend; period=1)

sol = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6, progress=true)

# Calcium transient is smaller than paper (550 nM vs 800 nM)
plot(sol, idxs=sys.Cai_sub_SL*1E6, tspan=(100, 102), ylabel="nM", lab="Ca_i")

# ICaL and INaCa are OK
@unpack INaCa, ICaL, ICaT, ICab = sys
plot(sol, idxs=[INaCa, ICaL, ICaT, ICab], tspan=(100, 102), ylabel="uA/uF")

plot(sol, idxs=sys.JCa_SL, tspan=(100, 102), ylabel="uM/ms")
