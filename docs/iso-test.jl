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
