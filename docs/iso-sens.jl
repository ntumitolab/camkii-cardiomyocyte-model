# # Sensitivity to ISO
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM, hil, Hz, hilr, second
Plots.default(lw=1.5)

# ## Setup b1AR system
@parameters ATP = 5000μM ISO = 0μM
sys = get_bar_sys(ATP, ISO; simplify=true)

#---
tend = 5000second
prob = ODEProblem(sys, [], tend)
alg = Rodas5P()

# Log scale for ISO concentration
iso = logrange(1e-4μM, 3μM, length=1001)

#---
prob_func = (prob, i, repeat) -> begin
    remake(prob, p=[ISO => iso[i]])
end
trajectories = length(iso)
callback = TerminateSteadyState(1e-10, 1e-10)
sol = solve(prob, alg; callback) ## warmup
sim = solve(EnsembleProblem(prob; prob_func, safetycopy=false), alg; trajectories, callback)

#---
"""Extract values from ensemble simulations by a symbol"""
extract(sim, k) = map(s -> s[k][end], sim)
"""Calculate Root Mean Square Error (RMSE)"""
rmse(fit) = sqrt(sum(abs2, fit.resid) / length(fit.resid))

#---
xopts = (xlims=(iso[begin], iso[end]), minorgrid=true, xscale=:log10, xlabel="ISO (μM)",)
plot(iso, extract(sim, sys.cAMP); lab="cAMP", ylabel="Conc. (μM)", legend=:topleft, xopts...)

#---
plot(iso, [extract(sim, sys.PKACI) extract(sim, sys.PKACII)], legend=:topleft, lab=["PKACI" "PKACII"]; xopts...)

#---
plot(iso, [extract(sim, sys.I1) extract(sim, sys.I1p) extract(sim, sys.I1p_PP1)], legend=:topleft, lab=["I1" "I1p" "I1p_PP1"]; xopts...)

#---
plot(iso, extract(sim, sys.PKACI / sys.RItot); lab="PKACI", ylabel="Activation fraction")
plot!(iso, extract(sim, sys.PKACII / sys.RIItot), lab="PKACII")
plot!(iso, extract(sim, sys.PP1 / sys.PP1totBA), lab="PP1", legend=:topleft; xopts...)

# ## Fitting active PKACI
@. model(x, p) = p[1] * x / (x + p[2]) + p[3]
xdata = iso
ydata = extract(sim, sys.PKACI / sys.RItot)
p0 = [0.3, 0.01μM, 0.08]
lb = [0.0, 0.0, 0.0]
pkac1_fit = curve_fit(model, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pkac1_coef = coef(pkac1_fit)

#---
println("PKACI")
println("Basal activity: ", pkac1_coef[3])
println("Activated activity: ", pkac1_coef[1])
println("Michaelis constant: ", pkac1_coef[2], " μM")
println("RMSE: ", rmse(pkac1_fit))

#---
ypred = model.(xdata, Ref(pkac1_coef))
p1 = plot(xdata, [ydata ypred], lab=["Full model" "Fitted"], line=[:dash :dot], title="PKACI", legend=:topleft; xopts...)
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="PKACI error (%)", lab=false, xopts...)
plot(p1, p2, layout=(2, 1), size=(600, 600))

# ## Fitting active PKACII
xdata = iso
ydata = extract(sim, sys.PKACII / sys.RIItot)
p0 = [0.4, 0.01μM, 0.2]
lb = [0.0, 0.0, 0.0]
pkac2_fit = curve_fit(model, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pkac2_coef = coef(pkac2_fit)

#---
println("PKACII")
println("Basal activity: ", pkac2_coef[3])
println("Activated activity: ", pkac2_coef[1])
println("Michaelis constant: ", pkac2_coef[2], " μM")
println("RMSE: ", rmse(pkac2_fit))

#---
ypred = model.(xdata, Ref(pkac2_coef))
p1 = plot(xdata, [ydata ypred], lab=["Full model" "Fitted"], line=[:dash :dot], title="PKACII", legend=:topleft; xopts...)
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="PKACII error (%)", lab=false, xopts...)
plot(p1, p2, layout=(2, 1), size=(600, 600))

# ## Least-square fitting of PP1 activity
@. model_pp1(x, p) = p[1] * p[2] / (x + p[2]) + p[3]
xdata = iso
ydata = extract(sim, sys.PP1 / sys.PP1totBA)
p0 = [0.1, 3e-3μM, 0.8]
lb = [0.0, 0.0, 0.0]
pp1_fit = curve_fit(model_pp1, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pp1_coef = coef(pp1_fit)

#---
println("PP1")
println("Repressible activity: ", pp1_coef[1])
println("Minimal activity: ", pp1_coef[3])
println("Repressive Michaelis constant: ", pp1_coef[2], " μM")
println("RMSE: ", rmse(pp1_fit))

#---
ypred = model_pp1.(xdata, Ref(pp1_coef))
p1 = plot(xdata, [ydata ypred], lab=["Full model" "Fitted"], line=[:dash :dot], title="PP1", legend=:topright; xopts...)
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="PP1 error (%)", lab=false, xopts...)
plot(p1, p2, layout=(2, 1), size=(600, 600))

# ## Fitting PLBp
xdata = iso
ydata = extract(sim, sys.PLBp / sys.PLBtotBA)
plot(xdata, ydata, title="PLBp fraction", lab=false; xopts...)

# First try: Hill function
@. model_plb(x, p) = p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-2μM, 1.0, 0.1]
lb = [0.5, 1e-9μM, 1.0, 0.0]
fit = curve_fit(model_plb, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pestim = coef(fit)

#---
println("PLBp")
println("Basal activity: ", pestim[4])
println("Activated activity: ", pestim[1])
println("Michaelis constant: ", pestim[2], " μM")
println("Hill coefficient: ", pestim[3])
println("RMSE: ", rmse(fit))

#---
ypred = model_plb.(xdata, Ref(pestim))
p1 = plot(xdata, [ydata ypred], lab=["Full model" "Fitted"], line=[:dash :dot], title="PLBp", legend=:topleft; xopts...)
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="PLBp error (%)", lab=false, xopts...)
plot(p1, p2, layout=(2, 1), size=(600, 600))

# Second try: analytic PLB from the quadratic equation
function plbp_analytic(iso)
    PLBtot = 106μM
    PKACItot = 1.18μM
    PP1tot = 0.89μM
    k_PKA_PLB = 54Hz / μM
    Km_PKA_PLB = 21μM
    k_PP1_PLB = 8.5Hz / μM
    Km_PP1_PLB = 7.0μM
    PKACI_basal = 0.0734  ## basal activity
    PKACI_activated = 0.1995
    PKACI_KM = 0.0139μM
    PP1_basal = 0.8927
    PP1_activated = 0.0492
    PP1_KI = 0.00637μM

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

ypred = plbp_analytic.(xdata)
fit = (resid= (ypred .- ydata) ./ ydata,)

p1 = plot(xdata, [ydata ypred], lab=["Full model" "Quadratic"], line=[:dash :dot], title="PLBp", legend=:topleft; xopts...)
p2 = plot(xdata, fit.resid .* 100; title="PLBp error (%)", lab=false, xopts...)
plot(p1, p2, layout=(2, 1), size=(600, 600))

#---
println("RMSE: ", rmse(fit))

# ## Fitting PLMp
xdata = iso
ydata = extract(sim, sys.PLMp / sys.PLMtotBA)
plot(xdata, ydata, title="PLMp fraction", lab=false; xopts...)

#---
@. model_plm(x, p) = p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-2μM, 1.0, 0.1]
lb = [0.5, 1e-9μM, 1.0, 0.0]
fit = curve_fit(model_plm, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pestim = coef(fit)

#---
println("PLMp")
println("Basal activity: ", pestim[4])
println("Activated activity: ", pestim[1])
println("Michaelis constant: ", pestim[2], " μM")
println("Hill coefficient: ", pestim[3])
println("RMSE: ", rmse(fit))

#---
ypred = model_plm.(xdata, Ref(pestim))
p1 = plot(xdata, [ydata ypred], lab=["Full model" "Fitted"], line=[:dash :dot], title="PLMp", legend=:topleft; xopts...)
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="PLMp error (%)", lab=false, xopts...)
plot(p1, p2, layout=(2, 1), size=(600, 600))

# ## Fitting TnIp
xdata = iso
ydata = extract(sim, sys.TnIp / sys.TnItotBA)
plot(xdata, ydata, title="TnIp fraction", lab=false; xopts...)

#---
@. model_tni(x, p) = p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-2μM, 1.0, 0.1]
lb = [0.1, 1e-9μM, 1.0, 0.0]
fit = curve_fit(model_tni, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pestim = coef(fit)

#---
println("TnIp")
println("Basal activity: ", pestim[4])
println("Activated activity: ", pestim[1])
println("Michaelis constant: ", pestim[2], " μM")
println("Hill coefficient: ", pestim[3])
println("RMSE: ", rmse(fit))

#---
ypred = model_tni.(xdata, Ref(pestim))
p1 = plot(xdata, [ydata ypred], lab=["Full model" "Fitted"], line=[:dash :dot], title="TnIp", legend=:topleft; xopts...)
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="TnIp error (%)", lab=false, xopts...)
plot(p1, p2, layout=(2, 1), size=(600, 600))

# ## Fitting LCCap
xdata = iso
ydata = extract(sim, sys.LCCap / sys.LCCtotBA)
plot(xdata, ydata, title="LCCap fraction", lab=false; xopts...)

#---
@. model_lcc(x, p) = p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-2μM, 1.0, 0.1]
lb = [0.1, 1e-9μM, 0.5, 0.0]
fit = curve_fit(model_lcc, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pestim = coef(fit)

#---
println("LCCap")
println("Basal activity: ", pestim[4])
println("Activated activity: ", pestim[1])
println("Michaelis constant: ", pestim[2], " μM")
println("Hill coefficient: ", pestim[3])
println("RMSE: ", rmse(fit))

#---
ypred = model_lcc.(xdata, Ref(pestim))
p1 = plot(xdata, [ydata ypred], lab=["Full model" "Fitted"], line=[:dash :dot], title="LCCap", legend=:topleft; xopts...)
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="LCCap error (%)", lab=false, xopts...)
plot(p1, p2, layout=(2, 1), size=(600, 600))

# ## Fitting LCCbp
xdata = iso
ydata = extract(sim, sys.LCCbp / sys.LCCtotBA)
plot(xdata, ydata, title="LCCbp fraction", lab=false; xopts...)

#---
@. model_lcc(x, p) = p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-2μM, 1.0, 0.1]
lb = [0.1, 1e-9μM, 0.5, 0.0]
fit = curve_fit(model_lcc, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pestim = coef(fit)

#---
println("LCCbp")
println("Basal activity: ", pestim[4])
println("Activated activity: ", pestim[1])
println("Michaelis constant: ", pestim[2], " μM")
println("Hill coefficient: ", pestim[3])
println("RMSE: ", rmse(fit))

#---
ypred = model_lcc.(xdata, Ref(pestim))
p1 = plot(xdata, [ydata ypred], lab=["Full model" "Fitted"], line=[:dash :dot], title="LCCbp", legend=:topleft; xopts...)
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="LCCbp error (%)", lab=false, xopts...)
plot(p1, p2, layout=(2, 1), size=(600, 600))

# ## Fitting KURp
xdata = iso
ydata = extract(sim, sys.KURp / sys.IKurtotBA)
plot(xdata, ydata, title="LCCbp fraction", lab=false; xopts...)

#---
@. model_kur(x, p) = p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-2μM, 1.0, 0.1]
lb = [0.1, 1e-9μM, 0.5, 0.0]
fit = curve_fit(model_kur, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pestim = coef(fit)

#---
println("KURp")
println("Basal activity: ", pestim[4])
println("Activated activity: ", pestim[1])
println("Michaelis constant: ", pestim[2], " μM")
println("Hill coefficient: ", pestim[3])
println("RMSE: ", rmse(fit))

#---
ypred = model_lcc.(xdata, Ref(pestim))
p1 = plot(xdata, [ydata ypred], lab=["Full model" "Fitted"], line=[:dash :dot], title="KURp", legend=:topleft; xopts...)
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="KURp error (%)", lab=false, xopts...)
plot(p1, p2, layout=(2, 1), size=(600, 600))

# ## Fitting RyRp
xdata = iso
ydata = extract(sim, sys.RyR_PKAp)
plot(xdata, ydata, title="RyRp fraction", lab=false; xopts...)

#---
@. model(x, p) = p[1] * x / (x + p[2]) + p[3]
p0 = [0.3, 1e-2μM, 0.1]
lb = [0.0, 1e-9μM, 0.0]
fit = curve_fit(model, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pestim = coef(fit)

#---
println("RyRp")
println("Basal activity: ", pestim[3])
println("Activated activity: ", pestim[1])
println("Michaelis constant: ", pestim[2], " μM")
println("RMSE: ", rmse(fit))

#---
ypred = model.(xdata, Ref(pestim))
p1 = plot(xdata, [ydata ypred], lab=["Full model" "Fitted"], line=[:dash :dot], title="RyRp", legend=:topleft; xopts...)
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="RyRp error (%)", lab=false, xopts...)
plot(p1, p2, layout=(2, 1), size=(600, 600))
