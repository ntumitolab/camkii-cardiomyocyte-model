using DiffEqCallbacks

function build_stim_callbacks(sym, endtime; period = 1second, duty = 5e-4second,starttime=zero(endtime), strength=-80μAμF, baseline=0μAμF, proposeddt=0.01second)
    rise! = (integrator) -> begin
        integrator[sym] = strength
        set_proposed_dt!(integrator, proposeddt)
    end

    fall! = (integrator) -> begin
        integrator[sym] = baseline
        set_proposed_dt!(integrator, proposeddt)
    end

    riseevents = PresetTimeCallback(starttime:period:endtime, rise!)
    fallevents = PresetTimeCallback(starttime+duty:period:endtime, fall!)
    return CallbackSet(riseevents, fallevents)
end

function build_istim_events(starttime, endtime, period, duty, istim, strength, baseline=0)
    rise = [collect(starttime:period:endtime) => [istim ~ strength]]
    fall = [collect(starttime+duty:period:endtime) => [istim ~ baseline]]
    return [rise, fall]
end
