using ProgressLogging
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM, Hz, μm
Plots.default(lw=1.5)

# ## Setup the ODE system
sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true)
tend = 300.0
prob = ODEProblem(sys, [], tend, [sys.ktrCaSR => 5])
@unpack Istim = sys
alg = FBDF()

callback = build_stim_callbacks(Istim, tend; period=1)

sol = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6, progress=true)

# Calcium transient is smaller than paper (550 nM vs 800 nM)
plot(sol, idxs=[sys.Cai_sub_SL*1E6, sys.Cai_sub_SR*1E6], tspan=(100, 101), ylabel="nM", lab=["Ca (sub SL)" "Ca (sub SR)"])

# ICaL and INaCa are OK
# RyR flux is tiny in NRVM
@unpack INaCa, ICaL, ICaT, ICab, JCa_SL, JCa_SR = sys
plot(sol, idxs=[INaCa, ICaL, JCa_SL, JCa_SR], tspan=(100, 101), ylabel="uA/uF or uM/ms")

@unpack Jup, Jrel, Jtr = sys
plot(sol, idxs=[Jup, Jtr], tspan=(100, 102), ylabel="uM/ms")

# SR Ca release is reletively small (see JSR)
plot(sol, idxs=[sys.CaJSR, sys.CaNSR], tspan=(100, 104), ylabel="mM")

plot(sol, idxs=[sys.betaSR], tspan=(100, 104))

prob.ps[sys.ACAP_F]
