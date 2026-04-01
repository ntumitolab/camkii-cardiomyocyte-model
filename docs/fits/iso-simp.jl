# # b1AR system simplification
# Fitting sensitivity to ISO.
using ModelingToolkit
using OrdinaryDiffEq, SteadyStateDiffEq, DiffEqCallbacks
using Plots
using CurveFit
using Model
using Model: μM, hil, Hz, hilr, second
Plots.default(lw=1.5)

# ## Setup b1AR system
@parameters ATP = 5000μM ISO = 0μM
@time "Build system" @mtkcompile sys = Model.get_bar_sys(ATP, ISO)
@time "Build problem" prob = SteadyStateProblem(sys, [])

alg = DynamicSS(KenCarp47())
@time "Solve problem" sol = solve(prob, alg; abstol=1e-10, reltol=1e-10) ## warmup

#---
@unpack cAMP, fracPKACI, fracPKACII, fracPP1, fracPLBp, fracPLMp, TnI_PKAp, LCCa_PKAp, LCCb_PKAp, IKUR_PKAp, RyR_PKAp = sys
tracked_states = [cAMP, fracPKACI, fracPKACII, fracPP1, fracPLBp, fracPLMp, TnI_PKAp, LCCa_PKAp, LCCb_PKAp, IKUR_PKAp, RyR_PKAp]
stimap = Dict(k=>i for (i, k) in enumerate(tracked_states))

@time "Extract tracked states" sol[tracked_states]
# ## Ensemble simulations
# Log scale for ISO concentration.
iso = exp10.(range(log10(1e-4μM), log10(1μM), length=501))

@time "Solve problems" sim = map(iso) do i
    newprob = remake(prob, p=[ISO => i])
    solve(newprob, alg; abstol=1e-10, reltol=1e-10)[tracked_states]
end

yymat = reduce(hcat, sim)'

#---
xopts = (xlims=(iso[begin], iso[end]), minorgrid=true, xscale=:log10, xlabel="ISO (μM)",)
plot(iso, yymat[:, stimap[cAMP]]; lab="cAMP", ylabel="Conc. (μM)", legend=:topleft, xopts...)

#---
plot(iso, yymat[:, stimap[fracPKACI]]; lab="PKACI", ylabel="Activation fraction")
plot!(iso, yymat[:, stimap[fracPKACII]], lab="PKACII")
plot!(iso, yymat[:, stimap[fracPP1]], lab="PP1", legend=:topleft; xopts...)

# ## PKACI activity
pkaci_model(p, x) = @. p[1] * hil(x, p[2]) + p[3]
xdata = iso
ydata = yymat[:, stimap[fracPKACI]]
p0 = [0.3, 0.01μM, 0.08]
prob = NonlinearCurveFitProblem(pkaci_model, p0, xdata, ydata)
@time "Fit PKACI" fit_pkac1 = solve(prob)

#---
println("PKACI")
println("Basal activity: ", fit_pkac1.u[3])
println("Activated activity: ", fit_pkac1.u[1])
println("Michaelis constant: ", fit_pkac1.u[2], " μM")
println("RMSE: ", mse(fit_pkac1) |> sqrt)

#---
p1 = plot(xdata, [ydata predict(fit_pkac1)], lab=["Full model" "Fitted"], line=[:dash :dot], title="PKACI", legend=:topleft; xopts...)

#---
savefig("pkaci_fit.pdf")

#---
p2 = plot(xdata, residuals(fit_pkac1) ./ ydata .* 100; title="PKACI error (%)", lab=false, xopts...)

# ## PKACII activity
pkacii_model(p, x) = @. p[1] * hil(x, p[2]) + p[3]
xdata = iso
ydata = yymat[:, stimap[fracPKACII]]
p0 = [0.4, 0.01μM, 0.2]
prob = NonlinearCurveFitProblem(pkacii_model, p0, xdata, ydata)
@time fit_pkac2 = solve(prob)

#---
println("PKACII")
println("Basal activity: ", fit_pkac2.u[3])
println("Activated activity: ", fit_pkac2.u[1])
println("Michaelis constant: ", fit_pkac2.u[2], " μM")
println("RMSE: ", mse(fit_pkac2) |> sqrt)

#---
p1 = plot(xdata, [ydata predict(fit_pkac2)], lab=["Full model" "Fitted"], line=[:dash :dot], title="PKACII", legend=:topleft; xopts...)

#---
savefig("pkacii_fit.pdf")

#---
p2 = plot(xdata, residuals(fit_pkac2) ./ ydata .* 100; title="PKACII error (%)", lab=false, xopts...)

# ## PP1 activity
pp1_model(p, x) = @. p[1] * hil(p[2], x) + p[3]
xdata = iso
ydata = yymat[:, stimap[fracPP1]]
p0 = [0.1, 3e-3μM, 0.8]
prob = NonlinearCurveFitProblem(pp1_model, p0, xdata, ydata)
@time fit_pp1 = solve(prob)

#---
println("PP1")
println("Repressible activity: ", fit_pp1.u[1])
println("Minimal activity: ", fit_pp1.u[3])
println("Repressive Michaelis constant: ", fit_pp1.u[2], " μM")
println("RMSE: ", mse(fit_pp1) |> sqrt)

#---
p1 = plot(xdata, [ydata predict(fit_pp1)], lab=["Full model" "Fitted"], line=[:dash :dot], title="PP1", legend=:topright; xopts...)

#---
savefig("pp1_fit.pdf")

#---
p2 = plot(xdata, residuals(fit_pp1) ./ ydata .* 100; title="PP1 error (%)", lab=false, xopts...)

# ## PLBp
xdata = iso
ydata = yymat[:, stimap[fracPLBp]]
plot(xdata, ydata, title="PLBp fraction", lab=false; xopts...)

#---
plb_model(p, x) = @. p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-2μM, 1.0, 0.1]
prob = NonlinearCurveFitProblem(plb_model, p0, xdata, ydata)
@time fit_plb = solve(prob)

#---
println("PLBp")
println("Basal activity: ", fit_plb.u[4])
println("Activated activity: ", fit_plb.u[1])
println("Michaelis constant: ", fit_plb.u[2], " μM")
println("Hill coefficient: ", fit_plb.u[3])
println("RMSE: ", mse(fit_plb) |> sqrt)

#---
p1 = plot(xdata, [ydata predict(fit_plb)], lab=["Full model" "Fitted"], line=[:dash :dot], title="PLBp", legend=:topleft; xopts...)

#---
savefig("plbp_fit.pdf")

#---
p2 = plot(xdata, residuals(fit_plb) ./ ydata .* 100; title="PLBp error (%)", lab=false, xopts...)

# ## PLMp
xdata = iso
ydata = yymat[:, stimap[fracPLMp]]
plot(xdata, ydata, title="PLMp fraction", lab=false; xopts...)

#---
plm_model(p, x) = @. p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-2μM, 1.0, 0.1]
prob = NonlinearCurveFitProblem(plm_model, p0, xdata, ydata)
@time fit_plm = solve(prob)

#---
println("PLMp")
println("Basal activity: ", fit_plm.u[4])
println("Activated activity: ", fit_plm.u[1])
println("Michaelis constant: ", fit_plm.u[2], " μM")
println("Hill coefficient: ", fit_plm.u[3])
println("RMSE: ", mse(fit_plm) |> sqrt)

#---
p1 = plot(xdata, [ydata predict(fit_plm)], lab=["Full model" "Fitted"], line=[:dash :dot], title="PLMp", legend=:topleft; xopts...)

#---
savefig("plmp_fit.pdf")

#---
p2 = plot(xdata, residuals(fit_plm) ./ ydata .* 100; title="PLMp error (%)", lab=false, xopts...)

## TnIp
xdata = iso
ydata = yymat[:, stimap[TnI_PKAp]]
plot(xdata, ydata, title="TnIp fraction", lab=false; xopts...)

#---
tni_model(p, x) = @. p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-2μM, 1.0, 0.1]
lb = [0.1, 1e-9μM, 1.0, 0.0]
prob = NonlinearCurveFitProblem(tni_model, p0, xdata, ydata)
@time fit_tni = solve(prob)
#---
println("TnIp")
println("Basal activity: ", fit_tni.u[4])
println("Activated activity: ", fit_tni.u[1])
println("Michaelis constant: ", fit_tni.u[2], " μM")
println("Hill coefficient: ", fit_tni.u[3])
println("RMSE: ", mse(fit_tni) |> sqrt)

#---
p1 = plot(xdata, [ydata predict(fit_tni)], lab=["Full model" "Fitted"], line=[:dash :dot], title="TnIp", legend=:topleft; xopts...)

#---
savefig("tni_fit.pdf")

#---
p2 = plot(xdata, residuals(fit_tni) ./ ydata .* 100; title="TnIp error (%)", lab=false, xopts...)

# ## LCCap
xdata = iso
ydata = yymat[:, stimap[LCCa_PKAp]]
plot(xdata, ydata, title="LCCap fraction", lab=false; xopts...)

#---
lcc_model(p, x) = @. p[1] * hil(x, p[2]) + p[3]
p0 = [0.8, 1e-2μM, 0.1]
prob = NonlinearCurveFitProblem(lcc_model, p0, xdata, ydata)
@time fit_lcca = solve(prob)

#---
println("LCCap")
println("Basal activity: ", fit_lcca.u[3])
println("Activated activity: ", fit_lcca.u[1])
println("Michaelis constant: ", fit_lcca.u[2], " μM")
println("RMSE: ", mse(fit_lcca) |> sqrt)

#---
p1 = plot(xdata, [ydata predict(fit_lcca)], lab=["Full model" "Fitted"], line=[:dash :dot], title="LCCap", legend=:topleft; xopts...)

#---
savefig("lcca_fit.pdf")

#---
plot(xdata, residuals(fit_lcca) ./ ydata .* 100; title="LCCap error (%)", lab=false, xopts...)

# ## LCCbp
xdata = iso
ydata = yymat[:, stimap[LCCb_PKAp]]
plot(xdata, ydata, title="LCCbp fraction", lab=false; xopts...)

#---
p0 = [0.8, 1e-2μM, 0.1]
lb = [0.1, 1e-9μM, 0.0]
prob = NonlinearCurveFitProblem(lcc_model, p0, xdata, ydata)
@time fit_lccb = solve(prob)

#---
println("LCCbp")
println("Basal activity: ", fit_lccb.u[3])
println("Activated activity: ", fit_lccb.u[1])
println("Michaelis constant: ", fit_lccb.u[2], " μM")
println("RMSE: ", mse(fit_lccb) |> sqrt)

#---
p1 = plot(xdata, [ydata predict(fit_lccb)], lab=["Full model" "Fitted"], line=[:dash :dot], title="LCCbp", legend=:topleft; xopts...)

#---
savefig("lccbp_fit.pdf")

#---
plot(xdata, residuals(fit_lccb) ./ ydata .* 100; title="LCCbp error (%)", lab=false, xopts...)

# ## KURp
xdata = iso
ydata = yymat[:, stimap[IKUR_PKAp]]
plot(xdata, ydata, title="KURp fraction", lab=false; xopts...)

#---
kur_model(p, x) = @. p[1] * hil(x, p[2]) + p[3]
p0 = [0.8, 1e-2μM, 0.1]
lb = [0.1, 1e-9μM, 0.0]
prob = NonlinearCurveFitProblem(kur_model, p0, xdata, ydata)
@time fit_kur = solve(prob)
#---
println("KURp")
println("Basal activity: ", fit_kur.u[3])
println("Activated activity: ", fit_kur.u[1])
println("Michaelis constant: ", fit_kur.u[2], " μM")
println("RMSE: ", mse(fit_kur) |> sqrt)

#---
p1 = plot(xdata, [ydata predict(fit_kur)], lab=["Full model" "Fitted"], line=[:dash :dot], title="KURp", legend=:topleft; xopts...)

#---
savefig("kurp_fit.pdf")

#---
p2 = plot(xdata, residuals(fit_kur) ./ ydata .* 100; title="KURp error (%)", lab=false, xopts...)

# ## RyRp
xdata = iso
ydata = yymat[:, stimap[RyR_PKAp]]
plot(xdata, ydata, title="RyRp fraction", lab=false; xopts...)

#---
ryr_model(p, x) = @. p[1] * x / (x + p[2]) + p[3]
p0 = [0.3, 1e-2μM, 0.1]
lb = [0.0, 1e-9μM, 0.0]
prob = NonlinearCurveFitProblem(ryr_model, p0, xdata, ydata)
@time fit_ryr = solve(prob)
#---
println("RyRp")
println("Basal activity: ", fit_ryr.u[3])
println("Activated activity: ", fit_ryr.u[1])
println("Michaelis constant: ", fit_ryr.u[2], " μM")
println("RMSE: ", mse(fit_ryr) |> sqrt)

#---
p1 = plot(xdata, [ydata predict(fit_ryr)], lab=["Full model" "Fitted"], line=[:dash :dot], title="RyRp", legend=:topleft; xopts...)

#---
savefig("ryrp_fit.pdf")

#---
p2 = plot(xdata, residuals(fit_ryr) ./ ydata .* 100; title="RyRp error (%)", lab=false, xopts...)
