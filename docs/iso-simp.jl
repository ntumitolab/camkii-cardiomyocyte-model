#===
# b1AR system simplification

Fitting sensitivity to ISO.
===#
using ComponentArrays
using OrdinaryDiffEq, SteadyStateDiffEq
using Plots
using CurveFit
using CaMKIIModel
using CaMKIIModel: μM, hil, Hz, hilr, second

u0 = CaMKIIModel.get_bar_u0()
ps = CaMKIIModel.get_bar_ps()
prob = SteadyStateProblem(CaMKIIModel.bar_model!, u0, ps)

#---
alg = DynamicSS(KenCarp47())
iso = exp10.(range(log10(1e-4μM), log10(1μM), length=1001))
prob_func = (prob, i, repeat) -> (prob.p.ISO = iso[i]; prob)
sol = solve(prob, alg; abstol=1e-10, reltol=1e-10) ## warmup

@time "Solve problem" sim = solve(EnsembleProblem(prob; prob_func), alg; trajectories=length(iso), abstol=1e-10, reltol=1e-10);

#---
extract(sim, k) = map(s -> s.u[k], sim)

#---
xopts = (xlims=(iso[begin], iso[end]), minorgrid=true, xscale=:log10, xlabel="ISO (μM)",)

plot(iso, extract(sim, :cAMP); lab="cAMP", ylabel="Conc. (μM)", legend=:topleft, xopts...)

#---
pkaci = extract(sim, :PKACI) ./ ps.RItot
pkacii = extract(sim, :PKACII) ./ ps.RIItot
pp1 = map(s -> CaMKIIModel.bar_model(s.u, ps, nothing).PP1 / ps.PP1totBA, sim)

plot(iso, [pkaci pkacii pp1], lab=["PKACI" "PKACII" "PP1"], ylabel="Activation fraction", legend=:topleft; xopts...)

# ## PKACI activity
model(p, x) = @. p[1] * hil(x, p[2]) + p[3]
xdata = iso
ydata = extract(sim, :PKACI) ./ ps.RItot
p0 = [0.3, 0.01μM, 0.08]
prob = NonlinearCurveFitProblem(model, p0, xdata, ydata)
@time "Fit PKACI" fit_pkac1 = solve(prob)

#---
println("PKACI")
println("Basal activity: ", fit_pkac1.u[3])
println("Activated activity: ", fit_pkac1.u[1])
println("Michaelis constant: ", fit_pkac1.u[2], " μM")
println("RMSE: ", mse(fit_pkac1) |> sqrt)

#---
plot(xdata, [ydata predict(fit_pkac1)], lab=["Full model" "Fitted"], line=[:dash :dot], title="PKACI", legend=:topleft; xopts...)

#---
savefig("pkaci_fit.pdf")

#---
plot(xdata, residuals(fit_pkac1) ./ ydata .* 100; title="PKACI error (%)", lab=false, xopts...)

# ## PKACII activity
xdata = iso
ydata = extract(sim, :PKACII) ./ ps.RIItot
p0 = [0.4, 0.01μM, 0.2]
prob = NonlinearCurveFitProblem(model, p0, xdata, ydata)
@time fit_pkac2 = solve(prob)

#---
println("PKACII")
println("Basal activity: ", fit_pkac2.u[3])
println("Activated activity: ", fit_pkac2.u[1])
println("Michaelis constant: ", fit_pkac2.u[2], " μM")
println("RMSE: ", mse(fit_pkac2) |> sqrt)

#---
plot(xdata, [ydata predict(fit_pkac2)], lab=["Full model" "Fitted"], line=[:dash :dot], title="PKACII", legend=:topleft; xopts...)

#---
savefig("pkacii_fit.pdf")

#---
plot(xdata, residuals(fit_pkac2) ./ ydata .* 100; title="PKACII error (%)", lab=false, xopts...)

# ## PP1 activity
model_pp1(p, x) = @. p[1] * hil(p[2], x) + p[3]
xdata = iso
ydata = map(s -> CaMKIIModel.bar_model(s.u, ps, nothing).PP1 / ps.PP1totBA, sim)
p0 = [0.1, 3e-3μM, 0.8]
prob = NonlinearCurveFitProblem(model_pp1, p0, xdata, ydata)
@time fit_pp1 = solve(prob)
#---
println("PP1")
println("Repressible activity: ", fit_pp1.u[1])
println("Minimal activity: ", fit_pp1.u[3])
println("Repressive Michaelis constant: ", fit_pp1.u[2], " μM")
println("RMSE: ", mse(fit_pp1) |> sqrt)

#---
plot(xdata, [ydata predict(fit_pp1)], lab=["Full model" "Fitted"], line=[:dash :dot], title="PP1", legend=:topright; xopts...)

#---
savefig("pp1_fit.pdf")

#---
plot(xdata, residuals(fit_pp1) ./ ydata .* 100; title="PP1 error (%)", lab=false, xopts...)

# ## PLBp
xdata = iso
ydata = map(s -> CaMKIIModel.bar_model(s.u, ps, nothing).fracPLBp, sim)
model_plb(p, x) = @. p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-2μM, 1.0, 0.1]
prob = NonlinearCurveFitProblem(model_plb, p0, xdata, ydata)
@time fit_plb = solve(prob)

#---
println("PLBp")
println("Basal activity: ", fit_plb.u[4])
println("Activated activity: ", fit_plb.u[1])
println("Michaelis constant: ", fit_plb.u[2], " μM")
println("Hill coefficient: ", fit_plb.u[3])
println("RMSE: ", mse(fit_plb) |> sqrt)

#---
plot(xdata, [ydata predict(fit_plb)], lab=["Full model" "Fitted"], line=[:dash :dot], title="PLBp", legend=:topleft; xopts...)

#---
savefig("plbp_fit.pdf")

#---
plot(xdata, residuals(fit_plb) ./ ydata .* 100; title="PLBp error (%)", lab=false, xopts...)

# ## PLMp
xdata = iso
ydata = map(s -> CaMKIIModel.bar_model(s.u, ps, nothing).fracPLMp, sim)
model_plm(p, x) = @. p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-2μM, 1.0, 0.1]
prob = NonlinearCurveFitProblem(model_plm, p0, xdata, ydata)
@time fit_plm = solve(prob)

#---
println("PLMp")
println("Basal activity: ", fit_plm.u[4])
println("Activated activity: ", fit_plm.u[1])
println("Michaelis constant: ", fit_plm.u[2], " μM")
println("Hill coefficient: ", fit_plm.u[3])
println("RMSE: ", mse(fit_plm) |> sqrt)

#---
plot(xdata, [ydata predict(fit_plm)], lab=["Full model" "Fitted"], line=[:dash :dot], title="PLMp", legend=:topleft; xopts...)

#---
savefig("plmp_fit.pdf")

#---
plot(xdata, residuals(fit_plm) ./ ydata .* 100; title="PLMp error (%)", lab=false, xopts...)

# ## TnIp
xdata = iso
ydata = map(s -> CaMKIIModel.bar_model(s.u, ps, nothing).TnI_PKAp, sim)
model_tnip(p, x) = @. p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-2μM, 1.0, 0.1]
prob = NonlinearCurveFitProblem(model_tnip, p0, xdata, ydata)
@time fit_tnip = solve(prob)

#---
println("TnIp")
println("Basal activity: ", fit_tnip.u[4])
println("Activated activity: ", fit_tnip.u[1])
println("Michaelis constant: ", fit_tnip.u[2], " μM")
println("Hill coefficient: ", fit_tnip.u[3])
println("RMSE: ", mse(fit_tnip) |> sqrt)

#---
plot(xdata, [ydata predict(fit_tnip)], lab=["Full model" "Fitted"], line=[:dash :dot], title="TnIp", legend=:topleft; xopts...)

#---
savefig("tni_fit.pdf")

#---
plot(xdata, residuals(fit_tnip) ./ ydata .* 100; title="TnIp error (%)", lab=false, xopts...)

# ## LCCap
xdata = iso
ydata = map(s -> CaMKIIModel.bar_model(s.u, ps, nothing).LCCa_PKAp, sim)
model_lcc(p, x) = @. p[1] * hil(x, p[2]) + p[3]
p0 = [0.8, 1e-2μM, 0.1]
prob = NonlinearCurveFitProblem(model_lcc, p0, xdata, ydata)
@time fit_lcca = solve(prob)

#---
println("LCCap")
println("Basal activity: ", fit_lcca.u[3])
println("Activated activity: ", fit_lcca.u[1])
println("Michaelis constant: ", fit_lcca.u[2], " μM")
println("RMSE: ", mse(fit_lcca) |> sqrt)

#---
plot(xdata, [ydata predict(fit_lcca)], lab=["Full model" "Fitted"], line=[:dash :dot], title="LCCap", legend=:topleft; xopts...)

#---
savefig("lcca_fit.pdf")

#---
plot(xdata, residuals(fit_lcca) ./ ydata .* 100; title="LCCap error (%)", lab=false, xopts...)

# ## LCCbp
xdata = iso
ydata = map(s -> CaMKIIModel.bar_model(s.u, ps, nothing).LCCb_PKAp, sim)
p0 = [0.8, 1e-2μM, 0.1]
lb = [0.1, 1e-9μM, 0.0]
prob = NonlinearCurveFitProblem(model_lcc, p0, xdata, ydata)
@time fit_lccb = solve(prob)

#---
println("LCCbp")
println("Basal activity: ", fit_lccb.u[3])
println("Activated activity: ", fit_lccb.u[1])
println("Michaelis constant: ", fit_lccb.u[2], " μM")
println("RMSE: ", mse(fit_lccb) |> sqrt)

#---
plot(xdata, [ydata predict(fit_lccb)], lab=["Full model" "Fitted"], line=[:dash :dot], title="LCCbp", legend=:topleft; xopts...)

#---
savefig("lccbp_fit.pdf")

#---
plot(xdata, residuals(fit_lccb) ./ ydata .* 100; title="LCCbp error (%)", lab=false, xopts...)

# ## KURp
xdata = iso
ydata = map(s -> CaMKIIModel.bar_model(s.u, ps, nothing).IKUR_PKAp, sim)
model_kur(p, x) = @. p[1] * hil(x, p[2]) + p[3]
p0 = [0.8, 1e-2μM, 0.1]
lb = [0.1, 1e-9μM, 0.0]
prob = NonlinearCurveFitProblem(model_kur, p0, xdata, ydata)
@time fit_kur = solve(prob)

#---
println("KURp")
println("Basal activity: ", fit_kur.u[3])
println("Activated activity: ", fit_kur.u[1])
println("Michaelis constant: ", fit_kur.u[2], " μM")
println("RMSE: ", mse(fit_kur) |> sqrt)

#---
plot(xdata, [ydata predict(fit_kur)], lab=["Full model" "Fitted"], line=[:dash :dot], title="KURp", legend=:topleft; xopts...)
#---
savefig("kurp_fit.pdf")

#---
plot(xdata, residuals(fit_kur) ./ ydata .* 100; title="KURp error (%)", lab=false, xopts...)

# ## RyRp
xdata = iso
ydata = map(s -> CaMKIIModel.bar_model(s.u, ps, nothing).RyR_PKAp, sim)
model(p, x) = @. p[1] * x / (x + p[2]) + p[3]
p0 = [0.3, 1e-2μM, 0.1]
lb = [0.0, 1e-9μM, 0.0]
prob = NonlinearCurveFitProblem(model, p0, xdata, ydata)
@time fit_ryr = solve(prob)
#---
println("RyRp")
println("Basal activity: ", fit_ryr.u[3])
println("Activated activity: ", fit_ryr.u[1])
println("Michaelis constant: ", fit_ryr.u[2], " μM")
println("RMSE: ", mse(fit_ryr) |> sqrt)

#---
plot(xdata, [ydata predict(fit_ryr)], lab=["Full model" "Fitted"], line=[:dash :dot], title="RyRp", legend=:topleft; xopts...)

#---
savefig("ryrp_fit.pdf")

#---
plot(xdata, residuals(fit_ryr) ./ ydata .* 100; title="RyRp error (%)", lab=false, xopts...)
