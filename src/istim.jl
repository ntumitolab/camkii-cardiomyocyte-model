using DifferentialEquations

function build_stim_callbacks(starttime, endtime, period, duty, strength, sym; baseline=0, proposeddt=0.01)
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
