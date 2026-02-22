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
