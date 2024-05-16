#===
# Smooth calcium wave
===#
using Catalyst
using DifferentialEquations
using ModelingToolkit
using Plots
using CaMKIIModel: get_camkii_rn, μM, nM, second

# Reaction network
@parameters ROS=0μM period=1/3 ca_r=100nM ca_rise=550nM tstart=200.0second tend=300.0second
@variables t Ca(t)
rn = get_camkii_rn(Ca, ROS)
rnsys = convert(ODESystem, rn, remove_conserved=true)

# Periodic assymetric calcium pulses
function ca_wave(t;
    sharpness=0.38, assymetry=18, period=1 / 3second,
    ca_r=100nM, ca_rise=550nM, tstart=200.0second, tend=300.0second)
    tau = t / period
    x = assymetry * (tau - floor(tau))
    return ca_r + (ca_rise * (x * exp(1 - x))^sharpness) * (t >= tstart) * (t <= tend)
end

@register_symbolic ca_wave(t)

caeqs = [Ca ~ ca_wave(t; period, ca_r, ca_rise, tstart, tend)]
@named casys = ODESystem(caeqs, t)

odesys = extend(casys, rnsys)
sys = structural_simplify(odesys)


# ## First sumulation
tspan = (0.0, 400.0)
oprob = ODEProblem(sys, [], tspan)

alg = TRBDF2()  ## Rodas5 yields unstable results, so TRBDF2 is used
sol = solve(oprob, alg, tstops=200:1/3:300, abstol=1e-9, reltol=1e-9)  ## Lowered tolerance to make the reaction network respond to calcium changes

# Calcium waveform
@unpack Ca = sys
plot(sol, idxs=Ca, tspan=(200, 205), label="Calcium")

# CaMKII-CaM system
@unpack CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2 = sys

plot(
    sol,
    idxs=[CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2],
    size=(800, 600),
)

# Active CaMKII
plot(sol, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), label="Act. CaMKII", title="3Hz")
plot!(sol, idxs=Ca * 10, plotdensity=10_000, label="Ca * 10")

# Change frequency to 2Hz
oprob2 = remake(oprob, p=[period=>1/2second])
sol2 = solve(oprob2, alg, tstops=200:1/2:300, abstol=1e-9, reltol=1e-9)

# Calcium
plot(sol2, idxs=Ca, tspan=(200, 205))

# Components
plot(
    sol2,
    idxs=[CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2],
    size=(800, 600),
)

# Active CaMKII
plot(sol2, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), label="Act. CaMKII", title="2Hz")
plot!(sol2, idxs=Ca * 10, plotdensity=100_000, label="Ca * 10")

# Change frequency to 1Hz
oprob3 = remake(oprob, p=[period=>1second])
sol3 = solve(oprob3, alg, tstops=200:1:300, abstol=1e-9, reltol=1e-9)

# Calcium
plot(sol3, idxs=Ca, tspan=(200, 205))

# Components
plot(
    sol3,
    idxs=[CaM0, Ca2CaM_C, Ca2CaM_N, Ca4CaM, CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMK, CaMKP, CaMKP2],
    size=(800, 600),
)

# Active CaMKII
plot(sol3, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), label="Act. CaMKII", title="1Hz")
plot!(sol3, idxs=Ca * 10, plotdensity=100_000, label="Ca * 10")

# Frequency-dependent response
plot(sol, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), title="Act. CaMKII", label="3Hz")
plot!(sol2, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), label="2Hz")
plot!(sol3, idxs=sum([CaM0_CaMK, Ca2CaM_C_CaMK, Ca2CaM_N_CaMK, Ca4CaM_CaMK, CaM0_CaMKP, Ca2CaM_C_CaMKP, Ca2CaM_N_CaMKP, Ca4CaM_CaMKP, CaMKP, CaMKP2]), label="1Hz")
