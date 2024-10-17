# # Sensitivity to ISO
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM, hil, Hz, hilr
Plots.default(lw=1.5)

# ## Setup b1AR system
@parameters ATP = 5000μM ISO = 0μM
sys = get_bar_sys(ATP, ISO; simplify=true)

#---
tend = 10000.0
prob = ODEProblem(sys, [], tend)
alg = Rodas5()
callback = TerminateSteadyState()
# Log scale for ISO concentration
iso = logrange(1e-6μM, 10μM, 1001)

prob_func = (prob, i, repeat) -> begin
    remake(prob, p=[ISO => iso[i]])
end
trajectories = length(iso)
sol = solve(prob, alg; callback) ## warmup
sim = solve(EnsembleProblem(prob; prob_func, safetycopy=false), alg; trajectories, callback)

# ## time series
# It took over 1000 seconds to reach equilibrium.
plot(sol, idxs=sys.PKACI / sys.RItot)

#---
"""Extract values from ensemble simulations by a symbol"""
extract(sim, k) = map(s -> s[k][end], sim)

#---
xopts = (xlims=(iso[begin] * 1000, iso[end] * 1000), minorgrid=true, xscale=:log10, xlabel="ISO (μM)",)

plot(iso .* 1000, extract(sim, sys.cAMP); lab="cAMP", ylabel="Conc. (mM)", legend=:topleft, xopts...)

#---
plot(iso .* 1000, extract(sim, sys.PKACI / sys.RItot); lab="PKACI", ylabel="Activation fraction")
plot!(iso .* 1000, extract(sim, sys.PKACII / sys.RIItot), lab="PKACII", legend=:topleft; xopts...)

# ## Least-square fitting of PKACI activity
@. model(x, p) = p[1] * x / (x + p[2]) + p[3]
xdata = iso
ydata = extract(sim, sys.PKACI / sys.RItot)
p0 = [0.3, 1E-5, 0.08]
lb = [0.0, 0.0, 0.0]
pkac1_fit = curve_fit(model, xdata, ydata, p0; lower=lb)
pkac1_coef = coef(pkac1_fit)

println("PKACI")
println("Basal activity: ", pkac1_coef[3])
println("Maximal activity: ", pkac1_coef[1] + pkac1_coef[3])
println("Michaelis constant: ", pkac1_coef[2] * 1000, " μM")
confidence_inter = confint(pkac1_fit; level=0.95)

#---
ypred = model.(xdata, Ref(pkac1_coef))
plot(xdata .* 1000, ydata, lab="Full model", line=:dash, title="PKACI")
plot!(xdata .* 1000, ypred; lab="Fitted", line=:dot, legend=:topleft, xopts...)

#---
plot(xdata .* 1000, (ypred .- ydata) ./ ydata .* 100; title="PKACI error (%)", lab=false, xopts...)

# ## Least-square fitting of PKACII activity
xdata = iso
ydata = extract(sim, sys.PKACII / sys.RIItot)
p0 = [0.4, 0.01μM, 0.2]
lb = [0.0, 0.0, 0.0]
pkac2_fit = curve_fit(model, xdata, ydata, p0; lower=lb)
pkac2_coef = coef(pkac2_fit)

println("PKACII")
println("Basal activity: ", pkac2_coef[3])
println("Maximal activity: ", pkac2_coef[1] + pkac2_coef[3])
println("Michaelis constant: ", pkac2_coef[2] * 1000, " μM")
confidence_inter = confint(pkac2_fit; level=0.95)

#---
ypred = model.(xdata, Ref(pkac2_coef))
plot(xdata .* 1000, ydata; lab="Full model", line=:dash, title="PKACII")
plot!(xdata .* 1000, ypred; lab="Fitted", line=:dot, legend=:topleft, xopts...)

#---
plot(xdata .* 1000, (ypred .- ydata) ./ ydata .* 100; title="PKACII error (%)", lab=false, xopts...)

# ## Least-square fitting of PP1 activity
@. model_pp1(x, p) = p[1] * p[2] / (x + p[2]) + p[3]
xdata = iso
ydata = extract(sim, sys.PP1 / sys.PP1tot)
p0 = [0.1, 8e-6, 0.8]
lb = [0.0, 0.0, 0.0]
pp1_fit = curve_fit(model_pp1, xdata, ydata, p0; lower=lb)
pp1_coef = coef(pp1_fit)

println("PP1")
println("Basal activity: ", pp1_coef[1] + pp1_coef[3])
println("Minimal activity: ", pp1_coef[3])
println("Repressive Michaelis constant: ", pp1_coef[2] * 1000, " μM")
confidence_inter = confint(pp1_fit; level=0.95)

#---
ypred = model_pp1.(xdata, Ref(pp1_coef))
plot(xdata .* 1000, ydata, lab="Data", line=:dash, title="PP1")
plot!(xdata .* 1000, ypred, lab="Fitted", line=:dot; xopts...)

#---
plot(xdata .* 1000, (ypred .- ydata) ./ ydata, title="PP1 Error", xscale=:log10, lab=false; xopts...)

# ## Least-square fitting of PLBp
xdata = iso
ydata = extract(sim, sys.PLBp / sys.PLBtot)

plot(xdata .* 1000, ydata, title="PLBp activity", lab=false; xopts...)

# Hill function

@. model_plb(x, p) = p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-5, 1.0, 0.1]
lb = [0.5, 1e-9, 1.0, 0.0]
fit = curve_fit(model_plb, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pestim = coef(fit)

#---
confidence_inter = confint(fit; level=0.95)

#---
ypred = model_plb.(xdata, Ref(pestim))
plot(xdata .* 1000, ydata, lab="Full model", line=:dash, title="PLBp")
plot!(xdata .* 1000, ypred, lab="Fitted", line=:dot, legend=:topleft; xopts...)

#---
plot(xdata .* 1000, (ypred .- ydata) ./ ydata .* 100; lab="PLBp error (%)", xopts...)

# Analytic PLB by solving quadratic equation
function plbp_analytic(iso)
    PLBtot = 106μM
    PKACItot = 1.18μM
    PP1tot = 0.89μM
    k_PKA_PLB = 54Hz
    Km_PKA_PLB = 21μM
    k_PP1_PLB = 8.5Hz
    Km_PP1_PLB = 7.0μM
    PKACI_basal = 0.0831  ## basal activity
    PKACI_activated = 0.25603
    PKACI_KM = 0.0144μM
    PP1_basal = 0.82365
    PP1_activated = 0.1025
    PP1_KI = 0.008465μM

    ## Solve for Vf * x / (x + k1) = Vr * (1 - x) / (1 - x + k2)
    PKACI = PKACItot * (PKACI_basal + PKACI_activated * hil(iso, PKACI_KM))
    PP1 = PP1tot * (PP1_basal + PP1_activated * hilr(iso, PP1_KI))

    Vf = k_PKA_PLB * PKACI
    k1 = Km_PKA_PLB / PLBtot
    Vr = k_PP1_PLB * PP1
    k2 = Km_PP1_PLB / PLBtot
    A = 1 - (Vf / Vr)
    B = (Vf / Vr) * (1 + k2) + (k1 - 1)
    C = -k1
    plb = ifelse(A ≈ 0.0, k1 / (k1 + k2), (-B + sqrt(B^2 - 4 * A * C)) / 2A)
    plbp = 1 - plb
end

# quadratic prediction looks better
ypred = plbp_analytic.(xdata)
plot(xdata .* 1000, ydata, lab="Full model", line=:dash, title="PLBp")
plot!(xdata .* 1000, ypred, lab="Fitted (quadratic)", line=:dot, legend=:topleft; xopts...)

# 1.5% error
plot(xdata .* 1000, (ypred .- ydata) ./ ydata .* 100; lab="PLBp error (%)", xopts...)
