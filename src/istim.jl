function build_stim_callbacks(sym, endtime; period=1second, duty=5e-4second, starttime=zero(endtime), strength=-80μAμF, baseline=0μAμF, proposeddt=0.01second)
    rise! = (integrator) -> begin
        integrator.ps[sym] = strength
        set_proposed_dt!(integrator, proposeddt)
    end

    fall! = (integrator) -> begin
        integrator.ps[sym] = baseline
        set_proposed_dt!(integrator, proposeddt)
    end

    ts = starttime:period:endtime

    riseevents = PresetTimeCallback(ts, rise!)
    fallevents = PresetTimeCallback(ts .+ duty, fall!)
    return CallbackSet(riseevents, fallevents)
end

function build_istim_events(starttime, endtime, period, duty, istim, strength, baseline=0)
    ts = starttime:period:endtime
    rise = [collect(ts) => [istim ~ strength]]
    fall = [collect(ts .+ duty) => [istim ~ baseline]]
    return [rise, fall]
end
