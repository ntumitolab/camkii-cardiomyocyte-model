# # Simplifying b1AR system
# Fitting sensitivity to ISO.
using ModelingToolkit
using OrdinaryDiffEq, SteadyStateDiffEq, DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM, hil, Hz, hilr, second
Plots.default(lw=1.5)

# ## Setup b1AR system
@parameters ATP = 5000μM ISO = 0μM
@time "Build system" sys = get_bar_sys(ATP, ISO; simplify=true)
@time "Build problem" prob = SteadyStateProblem(sys, [])

# Log scale for ISO concentration.
alg = DynamicSS(Rodas5P())
iso = exp10.(range(log10(1e-4μM), log10(1μM), length=1001))
prob_func = (prob, i, repeat) -> remake(prob, p=[ISO => iso[i]])
trajectories = length(iso)
sol = solve(prob, alg; abstol=1e-10, reltol=1e-10) ## warmup
@time "Solve problem" sim = solve(EnsembleProblem(prob; prob_func, safetycopy=false), alg; trajectories, abstol=1e-10, reltol=1e-10)

#---
"""Extract values from ensemble simulations by a symbol"""
extract(sim, k) = map(s -> s[k], sim)
"""Calculate Root Mean Square Error (RMSE)"""
rmse(fit) = sqrt(sum(abs2, fit.resid) / length(fit.resid))

#---
xopts = (xlims=(iso[begin], iso[end]), minorgrid=true, xscale=:log10, xlabel="ISO (μM)",)
plot(iso, extract(sim, sys.cAMP); lab="cAMP", ylabel="Conc. (μM)", legend=:topleft, xopts...)

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
@time pkac1_fit = curve_fit(model, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
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

#---
savefig("pkaci_fit.png")
savefig("pkaci_fit.pdf")

#---
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="PKACI error (%)", lab=false, xopts...)

# ## Fitting active PKACII
xdata = iso
ydata = extract(sim, sys.PKACII / sys.RIItot)
p0 = [0.4, 0.01μM, 0.2]
lb = [0.0, 0.0, 0.0]
@time pkac2_fit = curve_fit(model, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
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

#---
savefig("pkacii_fit.png")
savefig("pkacii_fit.pdf")

#---
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="PKACII error (%)", lab=false, xopts...)

# ## Fitting PP1 activity
@. model_pp1(x, p) = p[1] * p[2] / (x + p[2]) + p[3]
xdata = iso
ydata = extract(sim, sys.PP1 / sys.PP1totBA)
p0 = [0.1, 3e-3μM, 0.8]
lb = [0.0, 0.0, 0.0]
@time pp1_fit = curve_fit(model_pp1, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
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

#---
savefig("pp1_fit.png")
savefig("pp1_fit.pdf")

#---
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="PP1 error (%)", lab=false, xopts...)

# ## Fitting PLBp
xdata = iso
ydata = extract(sim, sys.PLBp / sys.PLBtotBA)
plot(xdata, ydata, title="PLBp fraction", lab=false; xopts...)

#---
@. model_plb(x, p) = p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-2μM, 1.0, 0.1]
lb = [0.5, 1e-9μM, 1.0, 0.0]
@time fit = curve_fit(model_plb, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
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

#---
savefig("plbp_fit.png")
savefig("plbp_fit.pdf")

#---
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="PLBp error (%)", lab=false, xopts...)

# ## Fitting PLMp
xdata = iso
ydata = extract(sim, sys.PLMp / sys.PLMtotBA)
plot(xdata, ydata, title="PLMp fraction", lab=false; xopts...)

#---
@. model_plm(x, p) = p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-2μM, 1.0, 0.1]
lb = [0.5, 1e-9μM, 1.0, 0.0]
@time fit = curve_fit(model_plm, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
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

#---
savefig("plmp_fit.png")
savefig("plmp_fit.pdf")

#---
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="PLMp error (%)", lab=false, xopts...)

# ## Fitting TnIp
xdata = iso
ydata = extract(sim, sys.TnIp / sys.TnItotBA)
plot(xdata, ydata, title="TnIp fraction", lab=false; xopts...)

#---
@. model_tni(x, p) = p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-2μM, 1.0, 0.1]
lb = [0.1, 1e-9μM, 1.0, 0.0]
@time fit = curve_fit(model_tni, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
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

#---
savefig("tni_fit.png")
savefig("tni_fit.pdf")

#---
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="TnIp error (%)", lab=false, xopts...)

# ## Fitting LCCap
xdata = iso
ydata = extract(sim, sys.LCCap / sys.LCCtotBA)
plot(xdata, ydata, title="LCCap fraction", lab=false; xopts...)

#---
@. model_lcc(x, p) = p[1] * hil(x, p[2]) + p[3]
p0 = [0.8, 1e-2μM, 0.1]
lb = [0.1, 1e-9μM, 0.0]
@time fit = curve_fit(model_lcc, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pestim = coef(fit)

#---
println("LCCap")
println("Basal activity: ", pestim[3])
println("Activated activity: ", pestim[1])
println("Michaelis constant: ", pestim[2], " μM")
println("RMSE: ", rmse(fit))

#---
ypred = model_lcc.(xdata, Ref(pestim))
p1 = plot(xdata, [ydata ypred], lab=["Full model" "Fitted"], line=[:dash :dot], title="LCCap", legend=:topleft; xopts...)

#---
savefig("lcca_fit.png")
savefig("lcca_fit.pdf")

#---
plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="LCCap error (%)", lab=false, xopts...)

# ## Fitting LCCbp
xdata = iso
ydata = extract(sim, sys.LCCbp / sys.LCCtotBA)
plot(xdata, ydata, title="LCCbp fraction", lab=false; xopts...)

#---
p0 = [0.8, 1e-2μM, 0.1]
lb = [0.1, 1e-9μM, 0.0]
@time fit = curve_fit(model_lcc, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pestim = coef(fit)

#---
println("LCCbp")
println("Basal activity: ", pestim[3])
println("Activated activity: ", pestim[1])
println("Michaelis constant: ", pestim[2], " μM")
println("RMSE: ", rmse(fit))

#---
ypred = model_lcc.(xdata, Ref(pestim))
p1 = plot(xdata, [ydata ypred], lab=["Full model" "Fitted"], line=[:dash :dot], title="LCCbp", legend=:topleft; xopts...)

#---
savefig("lccbp_fit.png")
savefig("lccbp_fit.pdf")

#---
plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="LCCbp error (%)", lab=false, xopts...)

# ## Fitting KURp
xdata = iso
ydata = extract(sim, sys.KURp / sys.IKurtotBA)
plot(xdata, ydata, title="KURp fraction", lab=false; xopts...)

#---
@. model_kur(x, p) = p[1] * hil(x, p[2]) + p[3]
p0 = [0.8, 1e-2μM, 0.1]
lb = [0.1, 1e-9μM, 0.0]
@time fit = curve_fit(model_kur, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
pestim = coef(fit)

#---
println("KURp")
println("Basal activity: ", pestim[3])
println("Activated activity: ", pestim[1])
println("Michaelis constant: ", pestim[2], " μM")
println("RMSE: ", rmse(fit))

#---
ypred = model_kur.(xdata, Ref(pestim))
p1 = plot(xdata, [ydata ypred], lab=["Full model" "Fitted"], line=[:dash :dot], title="KURp", legend=:topleft; xopts...)

#---
savefig("kurp_fit.png")
savefig("kurp_fit.pdf")

#---
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="KURp error (%)", lab=false, xopts...)

# ## Fitting RyRp
xdata = iso
ydata = extract(sim, sys.RyR_PKAp)
plot(xdata, ydata, title="RyRp fraction", lab=false; xopts...)

#---
@. model(x, p) = p[1] * x / (x + p[2]) + p[3]
p0 = [0.3, 1e-2μM, 0.1]
lb = [0.0, 1e-9μM, 0.0]
@time fit = curve_fit(model, xdata, ydata, p0; lower=lb, autodiff=:forwarddiff)
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

#---
savefig("ryrp_fit.png")
savefig("ryrp_fit.pdf")

#---
p2 = plot(xdata, (ypred .- ydata) ./ ydata .* 100; title="RyRp error (%)", lab=false, xopts...)
