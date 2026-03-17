# # Parameter fitting
# Fitting CaMKII parameters to experimental decay rates
using ModelingToolkit
using SymbolicIndexingInterface
using OrdinaryDiffEq
using Optimization
using OptimizationLBFGSB
using ADTypes
using ForwardDiff
using CaMKIIModel
using CaMKIIModel: second

# Setup simulation without CaMKII A2 state
tend = 210second
@time "Build system" @mtkcompile sys = build_neonatal_ecc_sys()
@unpack r_CaMK, kb_CaMKP, kphos_CaMK, kdeph_CaMK, k_P1_P2, k_P2_P1, Istim = sys
@time "Build problem" prob = ODEProblem(sys, [k_P1_P2 => 0], tend)
sol = solve(prob, FBDF(), saveat=[30second], save_everystep=false) ## Baseline simulation
stimstart=30second
# Hyperparameters
ps = (;
    stimstart=stimstart,
    experiment_tau=[16.48, 16.73, 17.65, 18.08],
    callbacks=[build_stim_callbacks(Istim, stimstart + duration * second; period=1second, starttime=stimstart) for duration in (15.0, 30.0, 60.0, 90.0)],
    kphos_dephos=prob.ps[kphos_CaMK] / prob.ps[kdeph_CaMK],
    alg=FBDF(),
    a0=sol(stimstart, idxs=sys.CaMKAct),
)

# Loss function
function loss(x, ps)
    _p = remake(prob, p=[kb_CaMKP => exp10(x[1]), kdeph_CaMK => exp10(x[2]), kphos_CaMK => ps.kphos_dephos * exp10(x[2])])
    errors = 0.0
    for i in 1:4
        duration = (15.0, 30.0, 60.0, 90.0)[i]
        stimend = ps.stimstart + duration * second
        sample_time = stimend + ps.experiment_tau[i] * second
        saveat = [ps.stimstart + duration * second, sample_time]
        sol = solve(_p, ps.alg; callback=ps.callbacks[i], saveat=saveat, maxiters=100000, verbose=false)
        SciMLBase.successful_retcode(sol) || return Inf
        ## Activity should be 1/e (36.8%) at the sample time (after one time constant in exponential decay)
        yexpected = ps.a0 + (sol(stimend, idxs=sys.CaMKAct) - ps.a0) * exp(-1)
        ysim = sol(sample_time, idxs=sys.CaMKAct)
        errors += (ysim - yexpected)^2
    end
    return errors
end

x0 = [log10(prob.ps[kb_CaMKP]), log10(prob.ps[kdeph_CaMK])]
@time loss(x0, ps)

#---
optf = OptimizationFunction(loss, ADTypes.AutoForwardDiff())
optprob = OptimizationProblem(optf, x0, ps, lb=[-0.5 + x0[1], -2 + x0[2]], ub=[1 + x0[1], 2 + x0[2]])
@time sol = solve(optprob, OptimizationLBFGSB.LBFGSB(); verbose=false)
sol.objective
#---
fit = Dict(kb_CaMKP => exp10(sol.u[1]), kdeph_CaMK => exp10(sol.u[2]), kphos_CaMK => ps.kphos_dephos * exp10(sol.u[2]))

println("Fitted parameters:")
println("Dissociation time of CaMKP --> CaMKA:" , 1e-3 / fit[kb_CaMKP], " seconds.")
println("Dephosphorylation time of CaMKA:" , 1e-3 / fit[kdeph_CaMK], " seconds.")
println("Phosphorylation rate of CaMKB:" , fit[kphos_CaMK] * 1000, " Hz")
