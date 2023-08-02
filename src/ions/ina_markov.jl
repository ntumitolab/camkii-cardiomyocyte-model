# I_Na: Colancy-Rudy Markov Model
using Catalyst

_α123(v, n1, d1, d2, b, k1, k2) = n1 / (d1 * exp(-(v + b)/k1) + d2 * exp(-(v + b)/k2))
_α1(v, P1a1, P2a1, P3a1, P4a1, P5a1, P6a1) = _α123(v, P1a1, P2a1, P5a1, P3a1, P4a1, P6a1)
_α2(v, P1a1, P2a1, P3a1, P4a2, P5a2, P6a1) = _α123(v, P1a1, P2a1, P5a2, P3a1, P4a2, P6a1)
_α3(v, P1a1, P2a1, P3a1, P4a3, P5a3, P6a1) = _α123(v, P1a1, P2a1, P5a3, P3a1, P4a3, P6a1)
_β123(v, n, b, k) = n * exp(-(v - b)/k)
_β1(v, P1b1, P3a1, P2b1) = _β123(v, P1b1, -P3a1, P2b1)
_β2(v, P1b2, P2b2, P2b1) = _β123(v, P1b2, P2b2, P2b1)
_β3(v, P1b3, P2b3, P2b1) = _β123(v, P1b3, P2b3, P2b1)
_α4(v, P1a4, P3a4, P4a4, P2a4) = inv(P1a4 * exp(-(v + P4a4)/P2a4) + P3a4)
_α5(v, P1a5, P4a4, P2a5) = _β123(v, P1a5, -P4a4, P2a5)
_β5(v, P1b5, P2b5, P4a4) = P1b5 + P2b5 * exp(v + P4a4)
_β6(v, P1b6, P2b6) = _β123(v, P1b6, 0, P2b6)
_α7(v, P1a7, P2a7) = _β123(v, P1a7, 0, -P2a7)
_β7(v, P1b7, P2b7) = _β123(v, P1b7, 0, P2b7)
_β4(α3, α4, α5, β3, β5) = (α3 * α4 * α5) / (β3 * β5)
_α6(v, P1a4, P3a4, P4a4, P2a4, P1a6) = _α4(v, P1a4, P3a4, P4a4, P2a4) / P1a6

function make_ina_markov_rn()
    return @reaction_network begin
        @variables v(t)
        @parameters P1a1 = 3.802
        @parameters P2a1 = 0.1027
        @parameters P3a1 = 2.5
        @parameters P4a1 = 17
        @parameters P5a1 = 0.20
        @parameters P6a1 = 150
        @parameters P4a2 = 15
        @parameters P5a2 = 0.23
        @parameters
        _α1(v, P1a1, P2a1, P3a1, P4a1, P5a1, P6a1), (ICNa3, CNa3, LCNa3) --> (ICNa2, CNa2, LCNa2)
        _β1(v, P1b1, P3a1, P2b1), (ICNa3, CNa3, LCNa3) <-- (ICNa2, CNa2, LCNa2)
        _α2(v, P1a1, P2a1, P3a1, P4a2, P5a2, P6a1), (ICNa2, CNa2, LCNa2) --> (IFNa, CNa1, LCNa1)
        _β2(v, P1b2, P2b2, P2b1), (ICNa2, CNa2, LCNa2) <-- (IFNa, CNa1, LCNa1)
        _α3(v, P1a1, P2a1, P3a1, P4a3, P5a3, P6a1), (CNa1, LCNa1) --> (ONa, LONa)
        _β3(v, P1b3, P2b3, P2b1), (CNa1, LCNa1) <-- (ONa, LONa)
        _α4(v, P1a4, P3a4, P4a4, P2a4), ONa --> IFNa
        _β4(_α3(v, P1a1, P2a1, P3a1, P4a3, P5a3, P6a1), _α4(v, P1a4, P3a4, P4a4, P2a4), _α5(v, P1a5, P4a4, P2a5), _β3(v, P1b3, P2b3, P2b1), _β5(v, P1b5, P2b5, P4a4)), ONa <-- IFNa
        (_α5(v, P1a5, P4a4, P2a5), _β5(v, P1b5, P2b5, P4a4)), IFNa <--> CNa1
        _α6(v, P1a4, P3a4, P4a4, P2a4, P1a6), IFNa --> I1Na
        _β6(v, P1b6, P2b6), IFNa <-- I1Na
        (_α7(v, P1a7, P2a7), _β7(v, P1b7, P2b7)), I1Na <--> I2Na
        (P1a8, P1b8), (CNa3, CNa2, CNa1, ONa) <--> (LCNa3, LCNa2, LCNa1, LONa)
    end
end

function set_ina_markov_rn!(rn; CKIIOE=0)
    CKIIflag = Bool(CKIIOE)
    setdefaults!(rn,[
        :P1a1=>3.802,
        :P2a1=>0.1027,
        :P3a1=>2.5,
        :P4a1=>17,
        :P5a1=>0.20,
        :P6a1=>150,
        :P4a2=>15,
        :P5a2=>0.23,
        :P4a3=>12,
        :P5a3=>0.25,
        :P1b1=>0.1917,
        :P2b1=>20.3,
        :P1b2=>0.2,
        :P2b2=>2.5,
        :P1b3=>0.22,
        :P2b3=>7.5,
        :P1a4=>0.188495,
        :P2a4=>16.6,
        :P3a4=>0.393956,
        :P4a4=>7,
        :P1a5=>7e-7,
        :P2a5=>ifelse(CKIIflag, 7.8, 7.2),
        :P1b5=>ifelse(CKIIflag, 0.00467, 0.0044),
        :P2b5=>2e-5,
        :P1a6=>100,
        :P1b6=>ifelse(CKIIflag, 6.5449e-7, 8.9554e-7),
        :P2b6=>11.3944,
        :P1a7=>ifelse(CKIIflag, 0.3377e-3, 0.487e-4),
        :P2a7=>23.2696,
        :P1b7=>ifelse(CKIIflag, 1.868e-4, 0.2868e-3),
        :P2b7=>35.9898,
        :P1a8=>ifelse(CKIIflag, 6.5e-6, 0.1e-7),
        :P1b8=>ifelse(CKIIflag, 3.8e-3, 9.8e-3),
    ])

    return rn
end
