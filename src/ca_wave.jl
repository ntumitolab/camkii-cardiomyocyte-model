# Calcium waves
using ModelingToolkit

"Periodic assymetric calcium pulses"
function ca_wave(t;
    sharpness=0.38, assymetry=18, period=1 / 3second,
    ca_r=100nM, ca_rise=550nM, tstart=200.0second, tend=300.0second)
    tau = t / period
    x = assymetry * (tau - floor(tau))
    return ca_r + (ca_rise * (x * exp(1 - x))^sharpness) * (t >= tstart) * (t <= tend)
end

"Exponential decay calcium model"
function ca_decay_sys(;
    carest=50nM,
    decay_calcium=10.0Hz,
)
    @parameters begin
        CaResting = carest
        dCa = decay_calcium
        dCaRev = 1
    end
    @variables t Ca(t)
    D = Differential(t)
    eq = [D(Ca) ~ -dCa * dCaRev * (Ca - CaResting * dCaRev)]
    @named osys = ODESystem(eq, t)
    return osys
end
