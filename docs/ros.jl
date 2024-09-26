# # ROS effects
using ModelingToolkit
using DifferentialEquations
using Plots
using CaMKIIModel

#---
sys = build_neonatal_ecc_sys(simplify=true)
tend = 300.0
prob = ODEProblem(sys, [], tend)

# ## No ROS
@unpack Istim = sys
callback = build_stim_callbacks(Istim, tend; period=1)
alg = FBDF()
sol = solve(prob, alg; callback, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

#---
plot(sol, idxs=sys.vm, tspan=(298, 300), title="Action potential")

#---
plot(sol, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(299, 300), title="Calcium transcient")

#---
plot(sol, idxs=sys.CaMKAct, title="Active CaMKII")

# ## ROS 1uM
prob2 = remake(prob, p=[sys.ROS => 1e-4])
sol2 = solve(prob2, alg; callback, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

#---
plot(sol2, idxs=sys.vm, tspan=(298, 300), title="Action potential")

#---
plot(sol2, idxs=[sys.Cai_sub_SR, sys.Cai_sub_SL, sys.Cai_mean], tspan=(298, 300), title="Calcium transcient")

#---
plot(sol2, idxs=sys.CaMKAct, title="Active CaMKII")

# ## Comparisons
plot(sol, idxs=sys.CaMKAct, title="Active CaMKII", lab="ROS=0")
plot!(sol2, idxs=sys.CaMKAct, lab="ROS=1uM")

#---
plot(sol, idxs=sys.vm, title="Action potential", lab="ROS=0", tspan=(298, 300))
plot!(sol2, idxs=sys.vm, lab="ROS=1uM", tspan=(298, 300))

#---
plot(sol, idxs=sys.Cai_mean, title="Calcium", lab="ROS=0", tspan=(298, 300))
plot!(sol2, idxs=sys.Cai_mean, lab="ROS=1uM", tspan=(298, 300))
