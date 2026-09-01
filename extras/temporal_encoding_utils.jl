#===
# Temporal encoding utilities

Shared helper functions for the CaMKII "temporal encoding" analyses (`fig11-temporal-map.jl`
and `fig12-temporal-pattern.jl`):

- `build_schedule_callback`: assemble multi-block pacing schedules (needed for
  pulse-pause-pulse protocols) out of the existing `build_stim_callbacks`.
- `build_condition_prob`: inject a genotype (WT / R275H-like mutant) and a
  neuromodulatory condition (ISO, ROS) into an `ODEProblem`.
- `summarize_camkii` / `get_decay_tau`: reduce a solution to the three readouts
  requested for the temporal response map (peak activation, area-under-curve,
  post-pacing decay time constant).

No new package dependencies are introduced beyond what is already listed in the
repo's `Project.toml` (CSV, DataFrames, CurveFit, DiffEqCallbacks, ModelingToolkit,
OrdinaryDiffEq(SDIRK), Plots, StatsPlots).

Include this file (via `include(joinpath(@__DIR__, "temporal_encoding_utils.jl"))`)
before using it in a figure script, exactly like the other `docs/fig*.jl` scripts
`using Model`, etc.
===#
using Model
using Model: second, Hz, μM
using CurveFit
using DiffEqCallbacks
using ModelingToolkit
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK

"""
R275H-like mutant parameter set, taken directly from `docs/fig10-mutation.jl` /
Table 3 of the manuscript: forward CaM-Ca2/Ca4 association is doubled relative to
WT and refit against the reduced one-step CaMKII-CaM binding model (basal binding
kfb_CaMK -> 0).
"""
const R275H_PARAMS = (
    kfa2_CaMK=0.3134,
    kmCa2_CaMK=0.4732μM,
    kfa4_CaMK=0.1155,
    kmCa4_CaMK=0.3722μM,
    kfb_CaMK=0.0,
)

"""
    build_condition_prob(prob0, sys; mutant=false, iso=0μM, ros=0μM, ica_scale_iso=2.5, tend=nothing)

Remake `prob0` with a genotype (WT / R275H-like mutant, via `mutant::Bool`) and a
neuromodulatory condition (isoproterenol concentration `iso`, ROS/H2O2 concentration
`ros`). `ica_scale_iso=2.5` matches the ICaL scaling used elsewhere in the manuscript
for 0.1 uM ISO (Figs 6 and 10). Optionally override the simulated `tend` (in the
model's internal time units, i.e. already multiplied by `second`) so that short
protocols in a sweep don't pay for the full sweep's maximum simulation length.
"""
function build_condition_prob(prob0, sys; mutant=false, iso=0.0μM, ros=0.0μM, ica_scale_iso=2.5, tend=nothing)
    ps = Pair[sys.ISO => iso, sys.ROS => ros]
    iso > 0.0μM && push!(ps, sys.ICa_scale_ISO => ica_scale_iso)
    if mutant
        append!(ps, [
            sys.kfa2_CaMK => R275H_PARAMS.kfa2_CaMK,
            sys.kmCa2_CaMK => R275H_PARAMS.kmCa2_CaMK,
            sys.kfa4_CaMK => R275H_PARAMS.kfa4_CaMK,
            sys.kmCa4_CaMK => R275H_PARAMS.kmCa4_CaMK,
            sys.kfb_CaMK => R275H_PARAMS.kfb_CaMK,
        ])
    end
    return tend === nothing ? remake(prob0, p=ps) : remake(prob0, p=ps, tspan=(0.0, tend))
end

"""
    build_schedule_callback(Istim, blocks; kwargs...)

Build a `CallbackSet` implementing a pacing schedule made of one or more
constant-frequency blocks. `blocks` is a vector of `(start, stop, period)` tuples,
all in the model's internal time units (i.e. already scaled by `second`). This lets
pulse-pause-pulse (or any other non-continuous) protocol be assembled from the
existing single-block `build_stim_callbacks` without modifying it.

Example: 30 s of 1 Hz pacing, a 30 s pause, then another 30 s of 1 Hz pacing,
starting at t = 30 s:

    blocks = [(30.0second, 60.0second, 1.0second), (90.0second, 120.0second, 1.0second)]
    cb = build_schedule_callback(Istim, blocks)
"""
function build_schedule_callback(Istim, blocks; kwargs...)
    cbs = [build_stim_callbacks(Istim, stop; period=period, starttime=start, kwargs...) for (start, stop, period) in blocks]
    return CallbackSet(cbs...)
end

"Trapezoidal integration of y(x) sampled at (possibly irregular) points x."
function trapz(x, y)
    s = 0.0
    @inbounds for i in 2:length(x)
        s += (y[i] + y[i-1]) * (x[i] - x[i-1]) / 2
    end
    return s
end

"""
    summarize_camkii(sol, sys, stimstart, stimend, poststim_window=60.0; sample_dt=1.0)

Sample `sys.CaMKAct` from `stimstart` to `stimend + poststim_window` (`stimstart`,
`stimend` in internal time units; `poststim_window`, `sample_dt` in real seconds)
and return the peak activation and the area-under-curve (activity x seconds) over
that window, plus the raw sampled trace for plotting.
"""
function summarize_camkii(sol, sys, stimstart, stimend, poststim_window=60.0; sample_dt=1.0)
    tend_s = (stimend - stimstart) / second + poststim_window
    ts_s = collect(0.0:sample_dt:tend_s)
    ys = sol(stimstart .+ ts_s .* second, idxs=sys.CaMKAct).u
    return (; peak=maximum(ys), auc=trapz(ts_s, ys), ts_s, ys)
end

"""
    get_decay_tau(sol, sys, stimend; window=50.0, step=5.0)

Fit a single exponential decay to `sys.CaMKAct` from `stimend` to `stimend + window`
(real seconds) and return the decay time constant tau in seconds. Mirrors the
approach already used for Figs 4, 5, 7, 10 and `decay-sens.jl`.
"""
function get_decay_tau(sol, sys, stimend; window=50.0, step=5.0)
    ts_s = collect(0.0:step:window)
    ys = sol(stimend .+ ts_s .* second, idxs=sys.CaMKAct).u
    fit = solve(CurveFitProblem(ts_s, ys), ExpSumFitAlgorithm(n=1, withconst=true))
    return inv(-fit.u.λ[])
end
