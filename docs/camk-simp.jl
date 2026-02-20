# # CaMKII system simplification
using ModelingToolkit
using OrdinaryDiffEq, SteadyStateDiffEq, DiffEqCallbacks
using Plots
using CurveFit
using CaMKIIModel
using CaMKIIModel: μM, hil, second, Hz
Plots.default(lw=1.5)

#===
## CaMKII sensitivity to calcium

Modeling CaM-Cax binding to CaMKII only. No phosphorylation or oxidation reactions.
===#
@parameters Ca = 0μM ROS = 0μM
@time "Build system" sys = get_camkii_sys(; Ca, ROS, simplify=true)
@time "Build problem" prob = SteadyStateProblem(sys, [sys.k_phosCaM => 0])

# Physiological cytosolic calcium concentrations ranges from 30nM to 10μM.
ca = exp10.(range(log10(0.03μM), log10(10μM), length=1001))
prob_func = (prob, i, repeat) -> (prob.ps[Ca] = ca[i]; prob)
alg = DynamicSS(KenCarp47())
sol0 = solve(prob, alg) ## warmup
@time "Solve problem" sim = solve(EnsembleProblem(prob; prob_func), alg; trajectories=length(ca))

#---
"""Extract values from ensemble simulations by a symbol"""
extract(sim, k) = map(s -> s[k], sim)

# Status of the CaMKII system across a range of calcium concentrations.
xopts = (xlabel="Ca (μM)", xscale=:log10, minorgrid=true, xlims=(ca[1], ca[end]))
let
    plot(ca, extract(sim, sys.CaM0), lab="CaM0", ylabel="Conc. (μM)"; xopts...)
    plot!(ca, extract(sim, sys.Ca2CaM_C), lab="Ca2CaM_C")
    plot!(ca, extract(sim, sys.Ca2CaM_N), lab="Ca2CaM_N")
    plot!(ca, extract(sim, sys.Ca4CaM), lab="Ca4CaM")
    plot!(ca, extract(sim, sys.CaMK), lab="CaMK")
    plot!(ca, extract(sim, sys.CaM0_CaMK), lab="CaM0_CaMK")
    plot!(ca, extract(sim, sys.Ca2CaM_C_CaMK), lab="Ca2CaM_C_CaMK")
    plot!(ca, extract(sim, sys.Ca2CaM_N_CaMK), lab="Ca2CaM_N_CaMK")
    plot!(ca, extract(sim, sys.Ca4CaM_CaMK), lab="Ca4CaM_CaMK", legend=:topleft)
end

#---
savefig("camkii_cam_binding.pdf")

# Active CaMKII is defined as CaMKII bound to CaMCax
CaMKAct = 1 - (sys.CaMK) / sys.CAMKII_T
println("Basal activity with 30nM Ca is ", sol0[CaMKAct][end])
xdata = ca
ydata = extract(sim, CaMKAct)
plot(xdata, ydata, label=false, title="Active CaMKII", ylims=(0, 0.5); xopts...)

#===
## Least-square fitting

Least-square fitting of steady state CaMKII activities against calcium concentrations.

### Old model

Using a single Hill function to fit steady state CaMKII activities.

$$
\text{CaMKII}_{act} = p_1 \frac{c^n}{c^n + p_2^n}+ p_3
$$
===#
model_camk_single_hill(p, x) = @. p[1] * hil(x, p[2], p[4]) + p[3]
p0 = [0.4, 1μM, 0.0, 2.0]
prob = NonlinearCurveFitProblem(model_camk_single_hill, p0, xdata, ydata)
@time fit = solve(prob)

println("Basal activity: ", fit.u[3])
println("Maximal activated activity: ", fit.u[1])
println("Half-saturation Ca concentration: ", fit.u[2], " μM")
println("Hill coefficient: ", fit.u[4])
println("RMSE: ", mse(fit) |> sqrt)

# Fit result and the original
plot(xdata, [ydata predict(fit)], lab=["Full model" "Fitted"], line=[:dash :dot], title="Single Hill function fit", legend=:topleft, ylabel="Bound CaMKII"; xopts...)

#===
### New model

Using two Hill functions to fit steady state CaMKII activities.

$$
\text{CaMKII}_{act} = p_1 \frac{c^2}{c^2 + p_2^2} + p_3 \frac{c^4}{c^4 + p_4^4} + p_5
$$
===#
model_camk(p, x) = @. p[1] * hil(x, p[2], 2) + p[3] * hil(x, p[4], 4) + p[5]
p0 = [0.4, 1μM, 0.2, 1μM, 0.0]
prob = NonlinearCurveFitProblem(model_camk, p0, xdata, ydata)
@time fit = solve(prob)

println("Basal activity: ", fit.u[5])
println("Maximal activity by 2 Ca binding): ", fit.u[1])
println("Half saturation concentration for 2 Ca binding: ", fit.u[2], " μM")
println("Maximal activity by 4 Ca binding): ", fit.u[3])
println("Half saturation concentration for 4 Ca binding: ", fit.u[4], " μM")
println("RMSE: ", mse(fit) |> sqrt)

#---
plot(xdata, [ydata predict(fit)], lab=["Full model" "Fitted"], line=[:dash :dot], title="Dual Hill function fit", legend=:topleft, ylabel="Bound CaMKII"; xopts...)

#---
savefig("camkii_act_fit.pdf")

# ### Polynomial fitting
prob = CurveFitProblem(xdata, ydata)
alg = RationalPolynomialFitAlgorithm(num_degree=2, den_degree=2)
@time sol = solve(prob, alg)

println("Numerator coefficients: ", sol.u[1:3])
println("Denominator coefficients: ", vcat(1.0, sol.u[4:5]))
println("RMSE: ", mse(sol) |> sqrt)

#---
plot(xdata, [ydata sol.(xdata)], lab=["Full model" "Fitted"], line=[:dash :dot], title="Rational polynomial fit", legend=:topleft, ylabel="Bound CaMKII"; xopts...)
