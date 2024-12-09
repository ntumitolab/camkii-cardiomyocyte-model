# # Sensitivity of CaMKII system to calcium
using ModelingToolkit
using OrdinaryDiffEq
using SteadyStateDiffEq
using DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM, hil, second
Plots.default(lw=1.5)

# ## Setup CaMKII system
# Model CaM/Calcium binding only. NO phosphorylation or oxidation.
@parameters Ca = 0μM ROS = 0μM
sys = get_camkii_sys(Ca; ROS, simplify=true)
equations(sys)

#---
prob = SteadyStateProblem(sys, [sys.k_phosCaM => 0])
alg = DynamicSS(Rodas5P())
ca = logrange(0.1μM, 10μM, length=1001)
prob_func = (prob, i, repeat) -> begin
    remake(prob, p=[Ca => ca[i]])
end
trajectories = length(ca)
sol0 = solve(prob, alg) ## warmup
sim = solve(EnsembleProblem(prob; prob_func, safetycopy=false), alg; trajectories)

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
plot!(ca, extract(sim, sys.Ca4CaM_CaMK), lab="Ca4CaM_CaMK", legend=:topleft, size=(600, 600))

#---
println("Basal activity without ca is ", sol0[sys.CaMKAct][end])
plot(ca, extract(sim, sys.CaMKAct), label=false, title="Active CaMKII", ylims=(0, 0.5); xopts...)

# ## Least-square fitting of steady state CaMKII activity
@. model_camk(x, p) = p[1] * hil(x, p[2], p[4]) + p[3]
xdata = ca
ydata = extract(sim, sys.CaMKAct)
p0 = [0.2, 1μM, 0.02, 2.2]
lb = [0.0, 1E-6μM, 0.0, 1.0]
fit = curve_fit(model_camk, xdata, ydata, inv.(ydata), p0; lower=lb, autodiff=:forwarddiff)

# Parameters
pestim = coef(fit)

#---
println("Basal activity: ", pestim[3])
println("Maximal activated activity: ", pestim[1])
println("Half-saturation Ca concentration: ", pestim[2], " μM")
println("Hill coefficient: ", pestim[4])

# Loss value
rmse(fit)

# Fit result
yestim = model_camk.(xdata, Ref(pestim))

p1 = plot(xdata, [ydata yestim], lab=["Full model" "Fitted"], line=[:dash :dot], title="CaMKII activity", legend=:topleft; xopts...)
p2 = plot(xdata, (yestim .- ydata) ./ ydata .* 100, title="Relative error (%)", xlabel="Ca (μM)", lab=false; xopts...)
plot(p1, p2, layout=(2, 1), size=(600, 600))

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
prob = ODEProblem(sys, [], tend)
prob_simp = ODEProblem(sys_simp, [], tend)

#---
@time sol = solve(prob, alg; callback)
@time sol_simp = solve(prob_simp, alg; callback)

#---
plot(sol, idxs=sys.CaMKAct, label= "Full model", title="Active CaMKII")
plot!(sol_simp, idxs=sys_simp.CaMKAct, label= "Simplified model", xlabel="Time (ms)")

#---
@unpack CaMKB, CaMKP, CaMKA, CaMKA2, CaMK = sys_simp
plot(sol_simp, idxs=[CaMKB, CaMKP, CaMKA, CaMKA2, CaMK])

#---
@unpack CaMK, CaMKP, CaMKP2, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP = sys
plot(sol, idxs=[sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK]), sum([CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP]), CaMKP, CaMKP2, CaMK])

#---
plot(sol, idxs=sys.CaMKAct, label= "Full model", title="Active CaMKII")
plot!(sol_simp, idxs=sys_simp.CaMKAct, label= "Simplified model", xlabel="Time (ms)", tspan=(295second,300.0second))
