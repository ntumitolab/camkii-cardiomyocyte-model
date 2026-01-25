# # CaMKII system simplification and sensitivity analysis
using ModelingToolkit
using OrdinaryDiffEq, SteadyStateDiffEq, DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM, hil, second, Hz
Plots.default(lw=1.5)

# ## CaMKII sensitivity to calcium
# Modeling CaM-Cax binding to CaMKII only. No phosphorylation or oxidation reactions.
@parameters Ca = 0μM ROS = 0μM
@time "Build system" sys = get_camkii_sys(; Ca, ROS, simplify=true)
@time "Build problem" prob = SteadyStateProblem(sys, [sys.k_phosCaM => 0])

# Physiological cytosolic calcium concentrations ranges from 30nM to 10μM.
ca = exp10.(range(log10(0.03μM), log10(10μM), 1001))
prob_func = (prob, i, repeat) -> begin
    remake(prob, p=[Ca => ca[i]])
end
trajectories = length(ca)
alg = DynamicSS(Rodas5P())
sol0 = solve(prob, alg; abstol=1e-10, reltol=1e-10) ## warmup
@time "Solve problem" sim = solve(EnsembleProblem(prob; prob_func, safetycopy=false), alg; trajectories, abstol=1e-10, reltol=1e-10)

#---
"""Extract values from ensemble simulations by a symbol"""
extract(sim, k) = map(s -> s[k], sim)
"""Calculate Root Mean Square Error (RMSE)"""
rmse(fit) = sqrt(sum(abs2, fit.resid) / length(fit.resid))

# Status of the CaMKII system across a range of calcium concentrations.
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

#---
savefig("camkii_cam_binding.png")
savefig("camkii_cam_binding.pdf")

# Active CaMKII is defined as CaMKII bound to CaM except CaM0 (apo CaM)
CaMKAct = 1 - (sys.CaMK + sys.CaM0_CaMK) / sys.CAMKII_T
println("Basal activity with 30nM Ca is ", sol0[CaMKAct][end])
xdata = ca
ydata = extract(sim, CaMKAct)
plot(xdata, ydata, label=false, title="Active CaMKII", ylims=(0, 0.5); xopts...)

# ## Least-square fitting
# Least-square fitting of steady state CaMKII activities against calcium concentrations.

#===
### Old model

Using a single Hill function to fit steady state CaMKII activities.

$$
\text{CaMKII}_{act} = p_1 \frac{c^n}{c^n + p_2^n}+ p_3
$$
===#
@. model_camk_single_hill(x, p) = p[1] * hil(x, p[2], p[4]) + p[3]
p0 = [0.4, 1μM, 0.0, 2.0]
lb = [0.0, 0.001μM, 0.0, 1.0]
@time fit = curve_fit(model_camk_single_hill, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pestim = coef(fit)

#---
println("Basal activity: ", pestim[3])
println("Maximal activated activity: ", pestim[1])
println("Half-saturation Ca concentration: ", pestim[2], " μM")
println("Hill coefficient: ", pestim[4])
println("RMSE: ", rmse(fit))

# Fit result and the original model
yestim = model_camk_single_hill.(xdata, Ref(pestim))
p1 = plot(xdata, [ydata yestim], lab=["Full model" "Fitted"], line=[:dash :dot], title="Single Hill function fit", legend=:topleft, ylabel="Bound CaMKII"; xopts...)

#===
### New model

Using two Hill functions to fit steady state CaMKII activities.

$$
\text{CaMKII}_{act} = p_1 \frac{c^2}{c^2 + p_2^2} + p_3 \frac{c^4}{c^4 + p_4^4} + p_5
$$
===#
@. model_camk(x, p) = p[1] * hil(x, p[2], 2) + p[3] * hil(x, p[4], 4) + p[5]
p0 = [0.4, 1μM, 0.0, 1μM, 0.0]
lb = [0.0, 0.001μM, 0.0, 0.001μM, 0.0]
@time fit = curve_fit(model_camk, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pestim = coef(fit)

#---
println("Basal activity: ", pestim[5])
println("Maximal activity by 2 Ca binding): ", pestim[1])
println("Half saturation concentration for 2 Ca binding: ", pestim[2], " μM")
println("Maximal activity by 4 Ca binding): ", pestim[3])
println("Half saturation concentration for 4 Ca binding: ", pestim[4], " μM")
println("RMSE: ", rmse(fit))

# Visualize fitting
yestim = model_camk.(xdata, Ref(pestim))
p1 = plot(xdata, [ydata yestim], lab=["Full model" "Fitted"], line=[:dash :dot], title="Dual Hill function fit", legend=:topleft, ylabel="Bound CaMKII"; xopts...)

#---
savefig("camkii_act_fit.png")
savefig("camkii_act_fit.pdf")
