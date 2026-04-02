# I_Na: Colancy-Rudy Markov Model
using Catalyst
using ModelingToolkit

function ina_markov_sys(CKIIOE=0)

    _α123(v, n1, d1, d2, b, k1, k2) = n1 / (d1 * exp(-(v + b) / k1) + d2 * exp(-(v + b) / k2))
    _β123(v, n, b, k) = n * exp(-(v - b) / k)

    @variables t v(t) (α(t))[1:8] (β(t))[1:8]
    @parameters P1a1 = 3.802
    @parameters P2a1 = 0.1027
    @parameters P3a1 = 2.5
    @parameters P4a1 = 17
    @parameters P5a1 = 0.20
    @parameters P6a1 = 150
    @parameters P4a2 = 15
    @parameters P5a2 = 0.23
    @parameters P4a3 = 12
    @parameters P5a3 = 0.25
    @parameters P1b1 = 0.1917
    @parameters P2b1 = 20.3
    @parameters P1b2 = 0.2
    @parameters P2b2 = 2.5
    @parameters P1b3 = 0.22
    @parameters P2b3 = 7.5
    @parameters P1a4 = 0.188495
    @parameters P2a4 = 16.6
    @parameters P3a4 = 0.393956
    @parameters P4a4 = 7
    @parameters P1a5 = 7e-7
    @parameters P2a5 = ifelse(CKIIOE == 1, 7.8, 7.2)
    @parameters P1b5 = ifelse(CKIIOE == 1, 0.00467, 0.0044)
    @parameters P2b5 = 2e-5
    @parameters P1a6 = 100
    @parameters P1b6 = ifelse(CKIIOE == 1, 6.5449e-7, 8.9554e-7)
    @parameters P2b6 = 11.3944
    @parameters P1a7 = ifelse(CKIIOE == 1, 0.3377e-3, 0.487e-4)
    @parameters P2a7 = 23.2696
    @parameters P1b7 = ifelse(CKIIOE == 1, 1.868e-4, 0.2868e-3)
    @parameters P2b7 = 35.9898
    @parameters P1a8 = ifelse(CKIIOE == 1, 6.5e-6, 0.1e-7)
    @parameters P1b8 = ifelse(CKIIOE == 1, 3.8e-3, 9.8e-3)

    rate_eqs = [
        α[1] ~ _α123(v, P1a1, P2a1, P5a1, P3a1, P4a1, P6a1),
        α[2] ~ _α123(v, P1a1, P2a1, P5a2, P3a1, P4a2, P6a1),
        α[3] ~ _α123(v, P1a1, P2a1, P5a3, P3a1, P4a3, P6a1),
        α[4] ~ 1 / (P1a4 * exp(-(v + P4a4) / P2a4) + P3a4),
        α[5] ~ _β123(v, P1a5, -P4a4, P2a5),
        α[6] ~ α[4] / P1a6,
        α[7] ~ _β123(v, P1a7, 0, -P2a7),
        α[8] ~ P1a8,
        β[1] ~ _β123(v, P1b1, -P3a1, P2b1),
        β[2] ~ _β123(v, P1b2, P2b2, P2b1),
        β[3] ~ _β123(v, P1b3, P2b3, P2b1),
        β[4] ~ (α[3] * α[4] * α[5]) / (β[3] * β[5]),
        β[5] ~ P1b5 + P2b5 * exp(v + P4a4),
        β[6] ~ _β123(v, P1b6, 0, P2b6),
        β[7] ~ _β123(v, P1b7, 0, P2b7),
        β[8] ~ P1b8
    ]

    rn = @reaction_network begin
        ($(α[1]), $(β[1])), (ICNa3, CNa3, LCNa3) <--> (ICNa2, CNa2, LCNa2)
        ($(α[2]), $(β[2])), (ICNa2, CNa2, LCNa2) <--> (IFNa, CNa1, LCNa1)
        ($(α[3]), $(β[3])), (CNa1, LCNa1) <--> (ONa, LONa)
        ($(α[4]), $(β[4])), ONa <--> IFNa
        ($(α[5]), $(β[5])), IFNa <--> CNa1
        ($(α[6]), $(β[6])), IFNa <--> I1Na
        ($(α[7]), $(β[7])), I1Na <--> I2Na
        ($(α[8]), $(β[8])), (CNa3, CNa2, CNa1, ONa) <--> (LCNa3, LCNa2, LCNa1, LONa)
    end

    osys = convert(System, rn, remove_conserved=true)
    extsys = extend(osys, System(rate_eqs, name=:rates))
    return extsys
end
