# # Test: rapid binding of Ca and calmodulin
# Assuming rapid equilibrium for CaM0 <--> CaM2 <--> Ca4CaM
using ModelingToolkit
using OrdinaryDiffEq
using SteadyStateDiffEq
using DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM, hil, Hz
Plots.default(lw=1.5)

#---
sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true)
sys_simp = build_neonatal_ecc_sys(simplify=true, reduce_iso=true, reduce_camk=true)

println("# of state variable in the original model: ", length(unknowns(sys)))
println("# of state variable in the simplified model: ", length(unknowns(sys_simp)))

#---
tend = 500.0
stimstart = 100.0
stimend = 300.0
alg = FBDF()
@unpack Istim = sys
callback = build_stim_callbacks(Istim, stimend; period=1, starttime=stimstart)
prob = ODEProblem(sys, [], tend)
prob_simp = ODEProblem(sys_simp, [], tend)

#---
@time sol = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6)
@time sol_simp = solve(prob_simp, alg; callback, abstol=1e-6, reltol=1e-6)

#---
plot(sol, idxs=sys.vm*1000, tspan=(295, 300), label= "Full", title="Action potential")
plot!(sol_simp, idxs=sys.vm*1000, tspan=(295, 300), label= "Simplified", line=:dash)

# Something is wrong
plot(sol, idxs=sys.CaMKAct, label= "Full model", title="Active CaMKII")
plot!(sol_simp, idxs=sys_simp.CaMKAct, label= "Simplified model", line=:dash)

# PCaM too high
@unpack CaM, KCaM, PCaM, CaMK = sys_simp
@unpack CaM0, CaM2C, CaM2N, CaM4, CaM0_CaMK, CaM2C_CaMK, CaM2N_CaMK, CaM4_CaMK, CaM0_CaMKP, CaM2C_CaMKP, CaM2N_CaMKP, CaM4_CaMKP = sys_simp
plot(sol_simp, idxs=[CaM, KCaM, PCaM, CaMK], tspan=(100, 105))

@unpack CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM = sys
plot(sol, idxs=[CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM], tspan=(100, 105))

#---
@parameters Ca = 0 ROS=0
ksys = CaMKIIModel.get_camkii_fast_ca_binding_sys(Ca; ROS)

for eq in equations(ksys)
    println(eq)
end
