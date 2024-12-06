# # Sensitivity of CaMKII system to calcium
using ModelingToolkit
using OrdinaryDiffEq
using SteadyStateDiffEq
using DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM, hil
Plots.default(lw=1.5)

# ## Setup CaMKII system
# Model CaM/Calcium binding only. NO phosphorylation or oxidation.
@parameters Ca = 0μM ROS = 0μM
sys = get_camkii_sys(Ca; ROS, simplify=true)
equations(sys)

#---
prob = SteadyStateProblem(sys, [sys.k_phosCaM => 0])
alg = DynamicSS(Rodas5P())
ca = logrange(0.1μM, 10μM, length=101)
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

# Loss value
rmse(fit)

# Fit result
yestim = model_camk.(xdata, Ref(pestim))

p1 = plot(xdata, [ydata yestim], lab=["Full model" "Fitted"], line=[:dash :dot], title="CaMKII activity", legend=:topleft; xopts...)
p2 = plot(xdata, (yestim .- ydata) ./ ydata .* 100, title="Relative error (%)", xlabel="Ca (μM)", lab=false; xopts...)
plot(p1, p2, layout=(2, 1), size=(600, 600))

# ## Least-square fitting of kinetic CaMKII activity
# Finding the Keq for the reaction: CaM + CaMK + 4Ca <---> CaMCa4_CaMKII
xdata = ca
ydata = extract(sim, sys.CaMKAct)

function model(x, p)
    A = 30μM
    B = 70μM
    k0 = p[4]
    k1 = p[2]
    kmca = p[1]
    nca = p[3]
    y = map(x) do ca
        keq = k1 * hil(ca, kmca, nca) + k0
        xterm = A + B + 1 / keq
        camkcam4 = 0.5 * (xterm - sqrt(xterm^2 - 4 * A * B))
        y = camkcam4 / B
    end
end

#---
p0 = [1μM, 1e4, 1.0, 0.5]
lb = [0.1μM, 0.0, 0.5, 0.0]
ub= [1000μM, Inf, 4.0, 1000.0]
fit = curve_fit(model, xdata, ydata, inv.(ydata), p0; lower=lb, upper=ub, autodiff=:forwarddiff)

#---
pestim = coef(fit)

# Loss value
rmse(fit)

# Fit result
yestim = model.(xdata, Ref(pestim))
p1 = plot(xdata, [ydata yestim], lab=["Full model" "Fitted"], line=[:dash :dot], title="CaMKII activity"; xopts...)
p2 = plot(xdata, (yestim .- ydata) ./ ydata .* 100, xlabel="Ca (μM)", title="Relative error (%)", lab=false; xopts...)
plot(p1, p2, layout=(2, 1), size=(600, 600))
