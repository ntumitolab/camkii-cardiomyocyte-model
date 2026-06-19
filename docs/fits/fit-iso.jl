using CurveFit
using DiffEqCallbacks
using DifferentialEquations
using ModelingToolkit
using OrdinaryDiffEqSDIRK
using Plots
using Measures
using SteadyStateDiffEq
using Model
using Model: Hz, hil, hilr, second, μM
Plots.default(lw=1.5)

# ## Setup b1AR system
@parameters ATP = 5000μM ISO = 0μM
@time "Build system" sys = Model.get_bar_sys(ATP, ISO) |> mtkcompile
@time "Build problem" prob = SteadyStateProblem(sys, [])
alg = DynamicSS(KenCarp47())
@time "Solve problem" sol = solve(prob, alg; abstol=1e-10, reltol=1e-10)
# ## Ensemble simulations

"""Extract values from ensemble simulations by a symbol"""
extract(sim, k) = map(s -> s[k], sim)

# Log scale for ISO concentration.
iso = logrange(1e-4μM, 1μM, 1001)
@time "Solve problems" sim = map(iso) do i
    newprob = remake(prob, p=[ISO => i])
    solve(newprob, alg; abstol=1e-10, reltol=1e-10)
end

#---
@unpack cAMP, fracPKACI, fracPKACII, fracPP1, fracPLBp, fracPLMp, TnI_PKAp, LCCa_PKAp, LCCb_PKAp, IKUR_PKAp, RyR_PKAp = sys
xopts = (xlims=(iso[begin], iso[end]), minorgrid=true, xscale=:log10, xlabel="ISO (μM)",)

plot(iso, extract(sim, cAMP), lab="cAMP", ylabel="Conc. (μM)"; xopts...)

#---
plot(iso, extract(sim, fracPKACI), lab="PKACI", ylabel="Active fraction")
plot!(iso, extract(sim, fracPKACII), lab="PKACII")
plot!(iso, extract(sim, fracPP1), lab="PP1", legend=:topleft; xopts...)

# ## PKACI activity
pkaci_model(p, x) = @. p[1] * hil(x, p[2]) + p[3]
xdata = iso
ydata = extract(sim, fracPKACI)
p0 = [0.3, 0.01μM, 0.08]
@time fit_pkac1 = solve(NonlinearCurveFitProblem(pkaci_model, p0, xdata, ydata))

println("PKACI")
println("Basal activity: ", fit_pkac1.u[3])
println("Activated activity: ", fit_pkac1.u[1])
println("Michaelis constant: ", fit_pkac1.u[2], " μM")
println("RMSE: ", mse(fit_pkac1) |> sqrt)

figs2A = plot(xdata, [ydata predict(fit_pkac1)], lab=["PKACI" "Fitted"], line=[:dash :dot], title="A", legend=:topleft, titlelocation=:left; xopts...)

#---
plot(xdata, residuals(fit_pkac1) ./ ydata .* 100; title="PKACI fitting error (%)", lab=false, xopts...)

# ## PKACII activity
pkacii_model(p, x) = @. p[1] * hil(x, p[2]) + p[3]
xdata = iso
ydata = extract(sim, fracPKACII)
p0 = [0.4, 0.01μM, 0.2]
@time fit_pkac2 = solve(NonlinearCurveFitProblem(pkacii_model, p0, xdata, ydata))

println("PKACII")
println("Basal activity: ", fit_pkac2.u[3])
println("Activated activity: ", fit_pkac2.u[1])
println("Michaelis constant: ", fit_pkac2.u[2], " μM")
println("RMSE: ", mse(fit_pkac2) |> sqrt)

figs2B= plot(xdata, [ydata predict(fit_pkac2)], lab=["PKACII" "Fitted"], line=[:dash :dot], title="B", titlelocation=:left, legend=:topleft; xopts...)

#---
plot(xdata, residuals(fit_pkac2) ./ ydata .* 100; title="PKACII fitting error (%)", lab=false, xopts...)

# ## PP1 activity
pp1_model(p, x) = @. p[1] * hil(p[2], x) + p[3]
xdata = iso
ydata = extract(sim, fracPP1)
p0 = [0.1, 3e-3μM, 0.8]
@time fit_pp1 = solve(NonlinearCurveFitProblem(pp1_model, p0, xdata, ydata))

println("PP1")
println("Repressible activity: ", fit_pp1.u[1])
println("Minimal activity: ", fit_pp1.u[3])
println("Repressive Michaelis constant: ", fit_pp1.u[2], " μM")
println("RMSE: ", mse(fit_pp1) |> sqrt)

figs2C = plot(xdata, [ydata predict(fit_pp1)], lab=["PP1" "Fitted"], line=[:dash :dot], title="C", titlelocation=:left, legend=:topright; xopts...)

#---
plot(xdata, residuals(fit_pp1) ./ ydata .* 100; title="PP1 fitting error (%)", lab=false, xopts...)

# ## PLBp
plb_model(p, x) = @. p[1] * hil(x, p[2], p[3]) + p[4]
xdata = iso
ydata = extract(sim, fracPLBp)
p0 = [0.8, 1e-2μM, 1.0, 0.1]
@time fit_plb = solve(NonlinearCurveFitProblem(plb_model, p0, xdata, ydata))

println("PLBp")
println("Basal activity: ", fit_plb.u[4])
println("Activated activity: ", fit_plb.u[1])
println("Michaelis constant: ", fit_plb.u[2], " μM")
println("Hill coefficient: ", fit_plb.u[3])
println("RMSE: ", mse(fit_plb) |> sqrt)

figs2D = plot(xdata, [ydata predict(fit_plb)], lab=["PLBp" "Fitted"], line=[:dash :dot], title="D", titlelocation=:left, legend=:topleft; xopts...)

#---
plot(xdata, residuals(fit_plb) ./ ydata .* 100; title="PLBp fitting error (%)", lab=false, xopts...)

# ## PLMp
plm_model(p, x) = @. p[1] * hil(x, p[2], p[3]) + p[4]
xdata = iso
ydata = extract(sim, fracPLMp)
p0 = [0.8, 1e-2μM, 1.0, 0.1]
@time fit_plm = solve(NonlinearCurveFitProblem(plm_model, p0, xdata, ydata))

println("PLMp")
println("Basal activity: ", fit_plm.u[4])
println("Activated activity: ", fit_plm.u[1])
println("Michaelis constant: ", fit_plm.u[2], " μM")
println("Hill coefficient: ", fit_plm.u[3])
println("RMSE: ", mse(fit_plm) |> sqrt)

figs2E = plot(xdata, [ydata predict(fit_plm)], lab=["PLMp" "Fitted"], line=[:dash :dot], title="E", titlelocation=:left, legend=:topleft; xopts...)

#---
plot(xdata, residuals(fit_plm) ./ ydata .* 100; title="PLMp fitting error (%)", lab=false, xopts...)

# ## LCCap
lcc_model(p, x) = @. p[1] * hil(x, p[2]) + p[3]
xdata = iso
ydata = extract(sim, LCCa_PKAp)
p0 = [0.8, 1e-2μM, 0.1]
@time fit_lcca = solve(NonlinearCurveFitProblem(lcc_model, p0, xdata, ydata))

println("LCCap")
println("Basal activity: ", fit_lcca.u[3])
println("Activated activity: ", fit_lcca.u[1])
println("Michaelis constant: ", fit_lcca.u[2], " μM")
println("RMSE: ", mse(fit_lcca) |> sqrt)

figs2F = plot(xdata, [ydata predict(fit_lcca)], lab=["LCCap" "Fitted"], line=[:dash :dot], title="F", titlelocation=:left, legend=:topleft; xopts...)

#---
plot(xdata, residuals(fit_lcca) ./ ydata .* 100; title="LCCap fitting error (%)", lab=false, xopts...)

# ## LCCbp
xdata = iso
ydata = extract(sim, LCCb_PKAp)
p0 = [0.8, 1e-2μM, 0.1]
@time fit_lccb = solve(NonlinearCurveFitProblem(lcc_model, p0, xdata, ydata))

println("LCCbp")
println("Basal activity: ", fit_lccb.u[3])
println("Activated activity: ", fit_lccb.u[1])
println("Michaelis constant: ", fit_lccb.u[2], " μM")
println("RMSE: ", mse(fit_lccb) |> sqrt)

figs2G = plot(xdata, [ydata predict(fit_lccb)], lab=["LCCbp" "Fitted"], line=[:dash :dot], title="G", titlelocation=:left, legend=:topleft; xopts...)

#---
plot(xdata, residuals(fit_lccb) ./ ydata .* 100; title="LCCbp fitting error (%)", lab=false, xopts...)

# ## TnIp
xdata = iso
ydata = extract(sim, TnI_PKAp)

tni_model(p, x) = @. p[1] * hil(x, p[2], p[3]) + p[4]
p0 = [0.8, 1e-2μM, 1.0, 0.1]
@time fit_tni = solve(NonlinearCurveFitProblem(tni_model, p0, xdata, ydata))

println("TnIp")
println("Basal activity: ", fit_tni.u[4])
println("Activated activity: ", fit_tni.u[1])
println("Michaelis constant: ", fit_tni.u[2], " μM")
println("Hill coefficient: ", fit_tni.u[3])
println("RMSE: ", mse(fit_tni) |> sqrt)

figs2H = plot(xdata, [ydata predict(fit_tni)], lab=["TnIp" "Fitted"], line=[:dash :dot], title="H", titlelocation=:left, legend=:topleft; xopts...)

#---
plot(xdata, residuals(fit_tni) ./ ydata .* 100; title="TnIp fitting error (%)", lab=false, xopts...)

# ## RyRp
xdata = iso
ydata = extract(sim, RyR_PKAp)
ryr_model(p, x) = @. p[1] * x / (x + p[2]) + p[3]
p0 = [0.3, 1e-2μM, 0.1]
@time fit_ryr = solve(NonlinearCurveFitProblem(ryr_model, p0, xdata, ydata))

println("RyRp")
println("Basal activity: ", fit_ryr.u[3])
println("Activated activity: ", fit_ryr.u[1])
println("Michaelis constant: ", fit_ryr.u[2], " μM")
println("RMSE: ", mse(fit_ryr) |> sqrt)

figs2I = plot(xdata, [ydata predict(fit_ryr)], lab=["RyRp" "Fitted"], line=[:dash :dot], title="I", titlelocation=:left, legend=:topleft; xopts...)

#---
plot(xdata, residuals(fit_ryr) ./ ydata .* 100; title="RyRp fitting error (%)", lab=false, xopts...)

# ## KURp
xdata = iso
ydata = extract(sim, IKUR_PKAp)
kur_model(p, x) = @. p[1] * hil(x, p[2]) + p[3]
p0 = [0.8, 1e-2μM, 0.1]
@time fit_kur = solve(NonlinearCurveFitProblem(kur_model, p0, xdata, ydata))

println("KURp")
println("Basal activity: ", fit_kur.u[3])
println("Activated activity: ", fit_kur.u[1])
println("Michaelis constant: ", fit_kur.u[2], " μM")
println("RMSE: ", mse(fit_kur) |> sqrt)

figs2J = plot(xdata, [ydata predict(fit_kur)], lab=["KURp" "Fitted"], line=[:dash :dot], title="J", titlelocation=:left, legend=:topleft; xopts...)

#---
plot(xdata, residuals(fit_kur) ./ ydata .* 100; title="KURp fitting error (%)", lab=false, xopts...)

# ## Save figure
plot(figs2A, figs2B, figs2C, figs2D, figs2E, figs2F, figs2G, figs2H, figs2I, figs2J, layout=(5, 2), size=(800, 1200), left_margin = 10mm)

#---
savefig("figS2.png")
savefig("figS2.pdf")
