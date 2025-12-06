# # CaMKII system simplification and sensitivity analysis
using ModelingToolkit
using OrdinaryDiffEq, SteadyStateDiffEq, DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM, hil, second, Hz
Plots.default(lw=1.5)

# ## CaMKII sensitivity to calcium
# Model CaM-Calcium binding to CaMKII only. No phosphorylation or oxidation reactions.
@parameters Ca = 0μM ROS = 0μM
sys = get_camkii_sys(; Ca, ROS, simplify=true)

#---
prob = SteadyStateProblem(sys, [sys.k_phosCaM => 0])
alg = DynamicSS(FBDF())

# Calcium ranges from 30nM to 10μM
ca = exp10.(range(log10(0.03μM), log10(10μM), 1001))
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

# Active CaMKII fraction is defined as CaMKII bound to CaM except CaM0 (no calcium)
CaMKAct = 1 - (sys.CaMK + sys.CaM0_CaMK) / sys.CAMKII_T
println("Basal activity with 30nM Ca is ", sol0[CaMKAct][end])
xdata = ca
ydata = extract(sim, CaMKAct)
plot(xdata, ydata, label=false, title="Active CaMKII", ylims=(0, 0.5); xopts...)

# ## Least-square fitting of steady state CaMKII activity against calcium concentration
@. model_camk(x, p) = p[1] * hil(x, p[2], 2) + p[3] * hil(x, p[4], 4) + p[5]
p0 = [0.4, 1μM, 0.0, 1μM, 0.0]
lb = [0.0, 0.001μM, 0.0, 0.001μM, 0.0]
fit = curve_fit(model_camk, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)

# Parameters
pestim = coef(fit)

#---
println("Basal activity: ", pestim[5])
println("Maximal activity by 2 Ca binding): ", pestim[1])
println("Half saturation concentration for 2 Ca binding: ", pestim[2], " μM")
println("Maximal activity by 4 Ca binding): ", pestim[3])
println("Half saturation concentration for 4 Ca binding: ", pestim[4], " μM")
println("RMSE: ", rmse(fit))

# Fit result and the original model
yestim = model_camk.(xdata, Ref(pestim))

p1 = plot(xdata, [ydata yestim], lab=["Full model" "Fitted"], line=[:dash :dot], title="CaMKII activity", legend=:topleft; xopts...)
