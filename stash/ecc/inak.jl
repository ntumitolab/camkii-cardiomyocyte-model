# I_nak: Na/K Pump Current

using ModelingToolkit

function get_inak_eqs(vm, Fjunc_nak, Fsl_nak, Najunc, Nasl, Nao, Ko, CKIIflag=0, DigitalisFlag=0)
    sigma = (exp(Nao/67.3)-1)/7
    fnak = 1/(1+0.1245*exp(-0.1*vm*FoRT)+0.0365*sigma*exp(-vm*FoRT))

    NaGainFlag = CKIIflag

    @parameters IbarNaK = 5 * (1 - 0.1 * NaGainFlag) * (1 - 0.5 * DigitalisFlag)   # [uA/uF] changed from rabbit (1.90719)

    @parameters KmKo = 1.5      # [mM]
    @parameters Q10NaK = 1.63
    @parameters Q10KmNai = 1.39
    @parameters KmNaip = let
        KmNaip = 19 # [mM] changed from rabbit (11)
        fracPKA_PLMo = 0.116738 # Derived quantity (PLM_PKAp(baseline)/PLMtot)
        fracPKA_PLMiso = 0.859251 # Derived quantity (PLM_PKAp(ISO)/PLMtot)
        kPKA_PLM=KmNaip*(1-0.7019)/(fracPKA_PLMiso/fracPKA_PLMo-1) # PLM_PKAp ISO
        KmNaip_PKA=-kPKA_PLM+kPKA_PLM*(PLM_PKAp/fracPKA_PLMo)
        KmNaip-KmNaip_PKA
    end

    @variables t I_nak_junc(t) I_nak_sl(t) I_nak(t)

    return [
        I_nak_junc ~ Fjunc_nak * IbarNaK * fnak * hil(Ko, KmKo) * hil(Najunc, KmNaip, 4)
        I_nak_sl ~ Fsl_nak * IbarNaK * fnak * hil(Ko, KmKo) * hil(Nasl, KmNaip, 4)
        I_nak ~ I_nak_junc + I_nak_sl
    ]
end
