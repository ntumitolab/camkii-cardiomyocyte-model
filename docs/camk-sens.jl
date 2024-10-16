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
tend=10000.0
prob = ODEProblem(sys, [], tend)
alg = Rodas5()
callback = TerminateSteadyState()
ca = exp10.(range(log10(1e-2μM), log10(10μM), length=1001))
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

# ## Least-square fitting of steady state CaMKII activity
hil_s(x, k, n) = sign(x*k) * (abs(x)^n) / (abs(x)^n + abs(k)^n)
@. model_camk(x, p) = p[1] * hil_s(x, p[2], p[3]) + 0.021
xdata = ca
ydata = extract(sim, sys.CaMKAct)

p0 = [0.4, 1e-3, 2.5]
lb = [0.0, 1e-9, 1.0]

fit = curve_fit(model_camk, xdata, ydata, p0; lower=lb)
pestim = coef(fit)

#---
confidence_inter = confint(fit; level=0.95)

#---
yestim = model_camk.(xdata, Ref(pestim))
plot(xdata, [ydata yestim], lab=["Full model" "Fitted"], line=[:dash :dot], title="CaMKII activity", xscale=:log10, legend=:topleft)

#---
plot(xdata, (yestim .- ydata) ./ ydata .* 100, title="Error", xscale=:log10)

# ## Least-square fitting of kinetic CaMKII activity
# CaM + CaMK + 4Ca <---> CaMCa4_CaMKII
extract(sim, k) = map(s->s[k][end], sim)
hil_s(x, k, n) = sign(x*k) * (abs(x)^n) / (abs(x)^n + abs(k)^n)
sqrt_s(x) = sign(x) * sqrt(abs(x))

xdata = ca
ydata = extract(sim, sys.CaMKAct)

function model(x, p)
    A = 30μM
    B = 70μM
    kmca = p[1]
    nca = p[2]
    kfb = p[3]
    kca = p[4]
    k0 = p[5]
    y = map(x) do ca
        keq = kfb * hil_s(ca, kmca, nca) + kca * ca + k0
        xterm = A + B + 1/keq
        camkcam4 = 0.5 * ( xterm - sqrt_s(xterm^2 - 4 * A * B))
        y = camkcam4 / B
    end
end

p0 = [1μM, 2.0, 1e7, 1000.0, 0.1]
lb = [0.0, 1.0, 0.0, 0.0, 0.0]
fit = curve_fit(model, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pestim = coef(fit)

# Fit result
yestim = model.(xdata, Ref(pestim))

plot(xdata .* 1000, [ydata yestim], lab=["Full model" "Fitted"], line=[:dash :dot], title="CaMKII activity", xlabel="Ca (μM)", xscale=:log10)

# 15% error
plot(xdata .* 1000, (yestim .- ydata) ./ ydata, xscale=:log10, xlabel="Ca (μM)", title="Relative error")
