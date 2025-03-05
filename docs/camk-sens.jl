# # CaMKII system analysis
using ModelingToolkit
using OrdinaryDiffEq, SteadyStateDiffEq, DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM, hil, second, Hz
Plots.default(lw=1.5)

# ## CaMKII sensitivity to calcium
# Model CaM/Calcium binding only. No phosphorylation or oxidation.
@parameters Ca = 0μM ROS = 0μM
sys = get_camkii_sys(Ca; ROS, simplify=true)

#---
prob = SteadyStateProblem(sys, [sys.k_phosCaM => 0])
alg = DynamicSS(Rodas5P())
ca = logrange(0.03μM, 10μM, length=1001)
prob_func = (prob, i, repeat) -> begin
    remake(prob, p=[Ca => ca[i]])
end
trajectories = length(ca)
sol0 = solve(prob, alg; abstol=1e-10, reltol=1e-10) ## warmup
sim = solve(EnsembleProblem(prob; prob_func, safetycopy=false), alg; trajectories, abstol=1e-10, reltol=1e-10)

#---
"""Extract values from ensemble simulations by a symbol"""
extract(sim, k) = map(s -> s[k], sim)
"""Calculate Root Mean Square Error (RMSE)"""
rmse(fit) = sqrt(sum(abs2, fit.resid) / length(fit.resid))

# Components
xopts = (xlabel="Ca (μM)", xscale=:log10, minorgrid=true, xlims=(ca[1], ca[end]))
plot(ca, extract(sim, sys.CaM0), lab="CaM0", ylabel="Conc. (μM)"; xopts...)
plot!(ca, extract(sim, sys.Ca2CaM_C), lab="Ca2CaM_C")
plot!(ca, extract(sim, sys.Ca2CaM_N), lab="Ca2CaM_N")
plot!(ca, extract(sim, sys.Ca4CaM), lab="Ca4CaM")
plot!(ca, extract(sim, sys.CaMK), lab="CaMK")
plot!(ca, extract(sim, sys.CaM0_CaMK), lab="CaM0_CaMK")
plot!(ca, extract(sim, sys.Ca2CaM_C_CaMK), lab="Ca2CaM_C_CaMK")
plot!(ca, extract(sim, sys.Ca2CaM_N_CaMK), lab="Ca2CaM_N_CaMK")
plot!(ca, extract(sim, sys.Ca4CaM_CaMK), lab="Ca4CaM_CaMK", legend=:topleft)

# We excluded CaM0_CaMK (bound but lost all calcium) from the active CaMKII fraction
CaMKAct = 1 - (sys.CaMK + sys.CaM0_CaMK) / sys.CAMKII_T
println("Basal activity without ca is ", sol0[CaMKAct][end])
xdata = ca
ydata = extract(sim, CaMKAct)
plot(xdata, ydata, label=false, title="Active CaMKII", ylims=(0, 0.5); xopts...)

# ## Least-square fitting of steady state CaMKII activity
@. model_camk(x, p) = p[1] * hil(x, p[2], p[4]) + p[3]
p0 = [0.4, 1μM, 0.0, 2.0]
lb = [0.0, 0.001μM, 0.0, 1.0]
fit = curve_fit(model_camk, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)

# Parameters
pestim = coef(fit)

#---
println("Basal activity: ", pestim[3])
println("Maximal activated activity: ", pestim[1])
println("Half-saturation Ca concentration: ", pestim[2], " μM")
println("Hill coefficient: ", pestim[4])
println("RMSE: ", rmse(fit))

# Fit result and the original model
yestim = model_camk.(xdata, Ref(pestim))

p1 = plot(xdata, [ydata yestim], lab=["Full model" "Fitted"], line=[:dash :dot], title="CaMKII activity", legend=:topleft; xopts...)

# ## The simplified CaMKII system
simpsys = CaMKIIModel.get_camkii_simp_sys(Ca; ROS, simplify=true)

# ## Comparing original and simplified model
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
prob = ODEProblem(sys, [], tend);
prob_simp = ODEProblem(sys_simp, [sys_simp.kphos_CaMK => 5Hz], tend);

#---
@time sol = solve(prob, alg; callback)
@time sol_simp = solve(prob_simp, alg; callback)

plot(sol, idxs=sys.CaMKAct, label= "Full model", title="Active CaMKII")
plot!(sol_simp, idxs=sys_simp.CaMKAct, label= "Simplified model", xlabel="Time (ms)")

#---
@unpack CaMKB, CaMKP, CaMKA, CaMKA2, CaMK = sys_simp
plot(sol_simp, idxs=[CaMKB, CaMKP, CaMKA, CaMKA2, CaMK], legend=:right)

#---
plot(sol, idxs=sys.CaMKAct, label= "Full model", title="Active CaMKII")
plot!(sol_simp, idxs=sys_simp.CaMKAct, label= "Simplified model", xlabel="Time (ms)", tspan=(295second,300.0second))
