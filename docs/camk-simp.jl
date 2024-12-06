# # Test: rapid binding of Ca and calmodulin
# Assuming rapid equilibrium for CaM0 <--> CaM2 <--> Ca4CaM
using ModelingToolkit
using OrdinaryDiffEq
using SteadyStateDiffEq
using DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM, hil, Hz, second
Plots.default(lw=1.5)

#---
sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true)
sys_simp = build_neonatal_ecc_sys(simplify=true, reduce_iso=true, reduce_camk=true)

println("# of state variable in the original model: ", length(unknowns(sys)))
println("# of state variable in the simplified model: ", length(unknowns(sys_simp)))

#---
tend = 500.0second
stimstart = 100.0second
stimend = 300.0second
alg = FBDF()
@unpack Istim = sys
callback = build_stim_callbacks(Istim, stimend; period=1second, starttime=stimstart)
prob = ODEProblem(sys, [], tend)
prob_simp = ODEProblem(sys_simp, [], tend)

#---
@time sol = solve(prob, alg; callback)
@time sol_simp = solve(prob_simp, alg; callback)

#---
plot(sol, idxs=sys.vm, label= "Full", title="Action potential")
plot!(sol_simp, idxs=sys.vm, tspan=(295second, 300second), label= "Simplified", line=:dash, xlabel="Time (ms)")

# Something is wrong
plot(sol, idxs=sys.CaMKAct, label= "Full model", title="Active CaMKII")
plot!(sol_simp, idxs=sys_simp.CaMKAct, label= "Simplified model", line=:dash, xlabel="Time (ms)")

# Phosphorylated CaMKII too high
@unpack CaM, KCaM, PCaM, CaMKP = sys_simp
@unpack CaM0, CaM2C, CaM2N, CaM4, CaM0_CaMK, CaM2C_CaMK, CaM2N_CaMK, CaM4_CaMK, CaM0_CaMKP, CaM2C_CaMKP, CaM2N_CaMKP, CaM4_CaMKP = sys_simp
plot(sol_simp, idxs=[CaMKP, KCaM, PCaM, CaMK], tspan=(100second, 105second))
