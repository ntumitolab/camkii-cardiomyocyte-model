# # Sensitivity of CaMKII system to calcium
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM, hil
Plots.default(lw=1.5)

# ## Setup CaMKII system
# Model CaM/Calcium binding only. NO phosphorylation or oxidation.
@parameters Ca = 0μM ROS = 0μM
sys = get_camkii_sys(Ca; ROS, phospho_rate=0, simplify=true)

#---
tend = 10000.0
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
extract(sim, k) = map(s -> s[k][end], sim)

# Components
plot(ca .* 1000, extract(sim, sys.CaM0), lab="CaM0", ylabel="Conc. (mM)", xlabel="Ca (μM)", xscale=:log10, minorgrid=true, xlims=(1e-2, 10))
plot!(ca .* 1000, extract(sim, sys.Ca2CaM_C), lab="Ca2CaM_C")
plot!(ca .* 1000, extract(sim, sys.Ca2CaM_N), lab="Ca2CaM_N")
plot!(ca .* 1000, extract(sim, sys.Ca4CaM), lab="Ca4CaM")
plot!(ca .* 1000, extract(sim, sys.CaMK), lab="CaMK")
plot!(ca .* 1000, extract(sim, sys.CaM0_CaMK), lab="CaM0_CaMK")
plot!(ca .* 1000, extract(sim, sys.Ca2CaM_C_CaMK), lab="Ca2CaM_C_CaMK")
plot!(ca .* 1000, extract(sim, sys.Ca2CaM_N_CaMK), lab="Ca2CaM_N_CaMK")
plot!(ca .* 1000, extract(sim, sys.Ca4CaM_CaMK), lab="Ca4CaM_CaMK", legend=:topleft, size=(800, 600))

# Active CaMKII
println("Basal activity without ca is ", sol0[sys.CaMKAct][end])
plot(ca .* 1000, extract(sim, sys.CaMKAct), lab="Active CaMKII", xlabel="Ca (μM)", xscale=:log10, minorgrid=true, xlims=(1e-2, 10))

# ## Least-square fitting of steady state CaMKII activity
hil_s(x, k, n) = sign(x * k) * (abs(x)^n) / (abs(x)^n + abs(k)^n)
@. model_camk(x, p) = p[1] * hil(x, p[2], p[4]) + p[3]
xdata = ca
ydata = extract(sim, sys.CaMKAct)

p0 = [0.4, 1e-3, 0.019, 2.2]
lb = [0.0, 1e-9, 0.0, 1.0]

fit = curve_fit(model_camk, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pestim = coef(fit)

# Loss value
sum(abs2, fit.resid)

#---
yestim = model_camk.(xdata, Ref(pestim))
plot(xdata .* 1000, [ydata yestim], lab=["Full model" "Fitted"], line=[:dash :dot], title="CaMKII activity", xscale=:log10, legend=:topleft, minorgrid=true, xlims=(1e-2, 10), xlabel="Ca (μM)")

#---
plot(xdata .* 1000, (yestim .- ydata) ./ ydata .* 100, title="Relative error (%)", xscale=:log10, minorgrid=true, xlims=(1e-2, 10), xlabel="Ca (μM)", lab=false)

# ## Least-square fitting of kinetic CaMKII activity
# Finding the Keq for the reaction CaM + CaMK + 4Ca <---> CaMCa4_CaMKII
extract(sim, k) = map(s -> s[k][end], sim)
hil_s(x, k, n) = sign(x * k) * (abs(x)^n) / (abs(x)^n + abs(k)^n)
sqrt_s(x) = sign(x) * sqrt(abs(x))

xdata = ca
ydata = extract(sim, sys.CaMKAct)

function model(x, p)
    A = 30μM
    B = 70μM
    hcoef = p[5]
    kmca = p[1]
    khill = p[2]
    kca = p[3]
    k0 = p[4]
    y = map(x) do ca
        keq = khill * hil(ca, kmca, hcoef) + kca * ca + k0
        xterm = A + B + 1 / keq
        camkcam4 = 0.5 * (xterm - sqrt(xterm^2 - 4 * A * B))
        y = camkcam4 / B
    end
end

p0 = [1μM, 1e7, 1000.0, 0.1, 2.7]
lb = [0.0, 0.0, 0.0, 0.0, 1.0]

fit = curve_fit(model, xdata, ydata, inv.(xdata), p0; lower=lb, autodiff=:forwarddiff)
pestim = coef(fit)

# Loss value
sum(abs2, fit.resid)

# Fit result
yestim = model.(xdata, Ref(pestim))

plot(xdata .* 1000, [ydata yestim], lab=["Full model" "Fitted"], line=[:dash :dot], title="CaMKII activity", xlabel="Ca (μM)", xscale=:log10, minorgrid=true, xlims=(1e-2, 10))

# < 6% error
maximum((yestim .- ydata) ./ ydata .* 100)

#---
plot(xdata .* 1000, (yestim .- ydata) ./ ydata .* 100, xscale=:log10, minorgrid=true, xlabel="Ca (μM)", xlims=(1e-2, 10), title="Relative error (%)", yticks=-6:6, lab=false)
