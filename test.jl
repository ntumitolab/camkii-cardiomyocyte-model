using ProgressLogging
using ModelingToolkit
using DifferentialEquations
using Plots
using CaMKIIModel
using CaMKIIModel: μM
Plots.default(lw=2)

@parameters ATP=5000μM ISO=0μM
sys = get_bar_sys(ATP, ISO; simplify=true)
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

"""Extract values from ensemble simulations by a symbol"""
extract(sim, k) = map(s->s[k][end], sim)

plot(iso, extract(sim, sys.cAMP), lab="cAMP")
plot(iso, extract(sim, sys.PKACI/sys.RItot), lab="PKACI")
plot!(iso, extract(sim, sys.PKACII/sys.RIItot), lab="PKACII")

@time sol = solve(prob, alg; callback, progress=true, abstol=1e-6, reltol=1e-6)
sol[end]

# ECC
sys = build_neonatal_ecc_sys(simplify=true)
tend = 1000.0
prob = ODEProblem(sys, [], tend)

@unpack Istim = sys
callback = build_stim_callbacks(Istim, tend; period=1)
alg = FBDF()
@time sol = solve(prob, alg; callback, progress=true, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

# Add ISO
prob2 = remake(prob, p=[sys.ISO => 1E-5])
@time sol2 = solve(prob2, alg; callback, progress=true, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

plot(sol, idxs=sys.Cai_mean, tspan=(999, 1000), lab="ISO=0", title="Ca transient")
plot!(sol2, idxs=sys.Cai_mean, tspan=(999, 1000), lab="ISO=0.1")

plot(sol, idxs=sys.PKACI, lab="ISO=0", title="PKACI")
plot!(sol2, idxs=sys.PKACI, lab="ISO=0.1")

plot(sol, idxs=sys.cAMP, lab="ISO=0", title="cAMP")
plot!(sol2, idxs=sys.cAMP, lab="ISO=0.1")

@unpack LRG, RG, Gs, GsaGDP, GsaGTP, AC_GsaGTP = sys
plot(sol2, idxs = [LRG, RG, Gs], tspan=(0, 1))

@unpack b1AR, LR, LRG, RG, b1AR_S464, b1AR_S301 = sys
plot(sol2, idxs = [b1AR, LR, LRG, RG, b1AR_S464, b1AR_S301], tspan=(0, 10), size=(800, 800))
