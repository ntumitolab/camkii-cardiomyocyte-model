# # Sensitivity of CaMKII system to calcium
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM, hil

# ## Setup CaMKII system
# Model CaM/Calcium binding only. NO phosphorylation or oxidation.
@parameters Ca=0μM ROS=0μM
sys = get_camkii_sys(Ca; ROS, phospho_rate=0, simplify=true)

#---
tend=100.0
prob = ODEProblem(sys, [], tend)
alg = Rodas5()
callback = TerminateSteadyState()
ca = exp10.(range(log10(1e-2μM), log10(10μM), length=10001))
prob_func = (prob, i, repeat) -> begin
    remake(prob, p=[Ca => ca[i]])
end
trajectories = length(ca)
sol0 = solve(remake(prob, p=[Ca => 0.0]), alg; callback) ## warmup
sim = solve(EnsembleProblem(prob; prob_func, safetycopy=false), alg; trajectories, callback)

"""Extract values from ensemble simulations by a symbol"""
extract(sim, k) = map(s->s[k][end], sim)


# Components
plot(ca .* 1000, extract(sim, sys.CaM0), lab="CaM0", ylabel="Conc. (mM)", xlabel="Ca (μM)", xscale=:log10)
plot!(ca .* 1000, extract(sim, sys.Ca2CaM_C), lab="Ca2CaM_C")
plot!(ca .* 1000, extract(sim, sys.Ca2CaM_N), lab="Ca2CaM_N")
plot!(ca .* 1000, extract(sim, sys.Ca4CaM), lab="Ca4CaM")
plot!(ca .* 1000, extract(sim, sys.CaMK), lab="CaMK")
plot!(ca .* 1000, extract(sim, sys.CaM0_CaMK), lab="CaM0_CaMK")
plot!(ca .* 1000, extract(sim, sys.Ca2CaM_C_CaMK), lab="Ca2CaM_C_CaMK")
plot!(ca .* 1000, extract(sim, sys.Ca2CaM_N_CaMK), lab="Ca2CaM_N_CaMK")
plot!(ca .* 1000, extract(sim, sys.Ca4CaM_CaMK), lab="Ca4CaM_CaMK", legend=:topleft, size=(800, 800))

# Active CaMKII
println("Basal activity when Ca=0 is ", sol0[sys.CaMKAct][end])
plot(ca .* 1000, extract(sim, sys.CaMKAct), lab="Active CaMKII", xlabel="Ca (μM)", xscale=:log10)

# ## Least-square fitting of CaMKII activity
hil_s(x, k, n) = sign(x*k) * (abs(x)^n) / (abs(x)^n + abs(k)^n)
@. model_camk(x, p) = p[1] * hil_s(x, p[2], p[3]) + 0.021
xdata = ca
ydata = extract(sim, sys.CaMKAct)

p0 = [0.4, 1e-3, 2.5]
lb = [0.0, 0.0, 0.1]

fit = curve_fit(model_camk, xdata, ydata, p0; lower=lb)
pestim = coef(fit)

#---
confidence_inter = confint(fit; level=0.95)

#---
yestim = model_camk.(xdata, Ref(pestim))
plot(xdata, [ydata yestim], lab=["Data" "Fitted"], line=[:dash :dot], title="CaMKII activity", xscale=:log10, legend=:topleft)

#---
plot(xdata, (yestim .- ydata) ./ ydata .* 100, title="Error", xscale=:log10)
