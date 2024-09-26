# # Sensitivity to ISO
using ModelingToolkit
using DifferentialEquations
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM
Plots.default(lw=2)

# ## Setup b1AR system
@parameters ATP=5000μM ISO=0μM
sys = get_bar_sys(ATP, ISO; simplify=true)

#---
tend=10000.0
prob = ODEProblem(sys, [], tend)
alg = Rodas5()
callback = TerminateSteadyState()
# Log scale for ISO concentration
iso = exp10.(range(log10(1e-10μM), log10(1μM), length=1001))

prob_func = (prob, i, repeat) -> begin
    remake(prob, p=[ISO => iso[i]])
end
trajectories = length(iso)
sol = solve(prob, alg; callback) ## warmup
sim = solve(EnsembleProblem(prob; prob_func, safetycopy=false), alg; trajectories, callback)

# ## time series
# It took over 1000 seconds to reach equilibrium.
plot(sol, idxs=sys.PKACI/sys.RItot)

#---
"""Extract values from ensemble simulations by a symbol"""
extract(sim, k) = map(s->s[k][end], sim)

#---
plot(iso .* 1000, extract(sim, sys.cAMP), lab="cAMP", ylabel="Conc. (mM)", xlabel="ISO (μM)")

#---
plot(iso .* 1000, extract(sim, sys.PKACI/sys.RItot), lab="PKACI", xlabel="ISO (μM)", ylabel="Activation fraction")
plot!(iso .* 1000, extract(sim, sys.PKACII/sys.RIItot), lab="PKACII")

# ## Least-square fitting of PKACI activity
@. model(x, p) = p[1] * x / (x + p[2]) + p[3]
xdata = iso
ydata = extract(sim, sys.PKACI/sys.RItot)
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
plot(xdata, ydata, lab="Data", line=:dash, title="PKACI", xscale=:log10)
plot!(x-> model(x, pkac1_coef), xdata, lab="Fitted", line=:dot, legend=:topleft)

# ## Least-square fitting of PKACII activity
xdata = iso
ydata = extract(sim, sys.PKACII/sys.RIItot)
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
plot(xdata, ydata, lab="Data", line=:dash, title="PKACII", xscale=:log10)
plot!(x-> model(x, pkac2_coef), xdata, lab="Fitted", line=:dot, legend=:topleft)

# ## Least-square fitting of PP1 activity
@. model_pp1(x, p) = p[1] * p[2] / (x + p[2]) + p[3]
xdata = iso
ydata = extract(sim, sys.PP1/sys.PP1tot)
p0 = [0.1, 8e-6, 0.8]
lb = [0.0, 0.0, 0.0]
pp1_fit = curve_fit(model_pp1, xdata, ydata, p0; lower=lb)
pp1_coef = coef(pp1_fit)

println("PP1")
println("Basal activity: ", pp1_coef[1] + pp1_coef[3])
println("Minimal activity: ",  pp1_coef[3])
println("Repressive Michaelis constant: ", pp1_coef[2] * 1000, " μM")
confidence_inter = confint(pp1_fit; level=0.95)

#---
plot(xdata, ydata, lab="Data", line=:dash, title="PP1", xscale=:log10)
plot!(x-> model_pp1(x, pp1_coef), xdata, lab="Fitted", line=:dot)

# ## Least-square fitting of PLBp
# Signed Hill function
hil_s(x, k, n) = sign(x*k) * (abs(x)^n) / (abs(x)^n + abs(k)^n)
xdata = iso
ydata = extract(sim, sys.PLBp/sys.PLBtot)
@. model_plb(x, p) = p[1] * hil_s(x, p[2], p[3]) + p[4]
p0 = [0.8, 2e-6, 1.0, 0.1]
lb = [0.0, 0.0, 0.1, -Inf]
fit = curve_fit(model_plb, xdata, ydata, p0; lower=lb)
pestim = coef(fit)

#---
confidence_inter = confint(fit; level=0.95)

# It does not fit as good as the previous ones because it's in fact a cubic equation.
plot(xdata, ydata, lab="Data", line=:dash, title="PLBp", xscale=:log10)
plot!(x-> model_plb(x, pestim), xdata, lab="Fitted", line=:dot, legend=:topleft)
#---
plot(xdata, (model_plb.(xdata, Ref(pestim)) .- ydata)./ydata, lab="Error", xscale=:log10)
