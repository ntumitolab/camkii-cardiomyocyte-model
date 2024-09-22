# # Sensitivity to ISO
using ModelingToolkit
using DifferentialEquations
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM, hil
Plots.default(lw=2)

#---
@parameters ATP=5000μM ISO=0μM
sys = get_bar_sys(ATP, ISO; simplify=true)

#---
tend=100000.0
prob = ODEProblem(sys, [], tend)
alg = Rodas5()
callback = TerminateSteadyState()
iso = range(0.0, 0.5μM, length=1001)
prob_func = (prob, i, repeat) -> begin
    remake(prob, p=[ISO => iso[i]])
end
trajectories = length(iso)
sim = solve(EnsembleProblem(prob; prob_func, safetycopy=false), alg; trajectories, callback)

# time series
# it took over 1000 seconds to reach equilibrium
plot(sim[1], idxs=sys.PKACI/sys.RItot)

#---
"""Extract values from ensemble simulations by a symbol"""
extract(sim, k) = map(s->s[k][end], sim)

#---
plot(iso .* 1000, extract(sim, sys.cAMP), lab="cAMP", ylabel="Conc. (mM)", xlabel="ISO (μM)")

#---
plot(iso .* 1000, extract(sim, sys.PKACI/sys.RItot), lab="PKACI", xlabel="ISO (μM)", ylabel="Activation fraction")
plot!(iso .* 1000, extract(sim, sys.PKACII/sys.RIItot), lab="PKACII")

# Least-square fitting of PKACI activity
@. model(x, p) = p[1] * x / (x + p[2]) + p[3]
xdata = iso
ydata = extract(sim, sys.PKACI/sys.RItot)
p0 = [0.3, 0.04, 0.08]
lb = [0.0, 0.0, 0.0]
fit = curve_fit(model, xdata, ydata, p0; lower=lb)
pestim = coef(fit)

println("PKACI")
println("Basal activity: ", pestim[3])
println("Maximal activity: ", pestim[1] + pestim[3])
println("Michaelis constant: ", pestim[2] * 1000, " μM")
confidence_inter = confint(fit; level=0.95)

#---
plot(xdata, ydata, lab="Data", line=:dash, title="PKACI")
plot!(x-> model(x, pestim), xdata, lab="Fitted", line=:dot)

# Least-square fitting of PKACII activity
xdata = iso
ydata = extract(sim, sys.PKACII/sys.RIItot)
p0 = [0.4, 0.04, 0.2]
fit = curve_fit(model, xdata, ydata, p0; lower=lb)
pestim = coef(fit)

println("PKACII")
println("Basal activity: ", pestim[3])
println("Maximal activity: ", pestim[1] + pestim[3])
println("Michaelis constant: ", pestim[2] * 1000, " μM")
confidence_inter = confint(fit; level=0.95)
#---

plot(xdata, ydata, lab="Data", line=:dash, title="PKACII")
plot!(x-> model(x, pestim), xdata, lab="Fitted", line=:dot)

# PLB
plot(iso .* 1000, extract(sim, sys.PLBp/sys.PLBtot), lab="PLBp", xlabel="ISO (μM)")

# Signed Hill function
hil_s(x, k, n) = sign(x*k) * (abs(x)^n) / (abs(x)^n + abs(k)^n)
@. model(x, p) = p[1] * hil_s(x, p[2], p[3]) + 0.089
xdata = iso
ydata = extract(sim, sys.PLBp/sys.PLBtot)
p0 = [0.9, 2e-6, 1.0]
lb = [0.0, 0.0, 0.1]
fit = curve_fit(model, xdata, ydata, p0; lower=lb)
pestim = coef(fit)

#---
confidence_inter = confint(fit; level=0.95)
#---
plot(xdata, ydata, lab="Data", line=:dash, title="PLBp")
plot!(x-> model(x, pestim), xdata, lab="Fitted", line=:dot)
