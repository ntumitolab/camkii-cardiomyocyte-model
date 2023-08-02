include("ions/ina_hh.jl")

function make_ecc_system(;
    freq=1,                    # [Hz] - CHANGE DEPENDING ON FREQUENCY
    ICa_MarkovFlag=1,          # Set ICa_MarkovFlag to 1 for Markov
    INa_MarkovFlag=0,          # Set INa_MarkovFlag to 1 for Markov INa, 0 otherwise
    NaClampFlag=0,             # Na clamp (set flag to 1, 0 otherwise)
    PLMkoFlag=0,               # PLM KO (set flag to 1, 0 otherwise)
    StrophFlag=0,              # Strophanthidin (dyastolic Na influx) (set flag to 1, 0 otherwise)
    CaffeineFlag=0,            # Caffeine-induced Ca transient (set flag to 1, 0 otherwise)
    DigitalisFlag=0,           # Digitalis (set flag to 1, 0 otherwise)
    CKIIflag=0,                # Set CKIIflag to 1 for overexpression, 0 otherwise
    ItoFlag=1,           # Set Ito to use either origl params (=0) or Grandi Params (=1)
    loop=0, # CaMKII-Na-Ca-CaMKII loop closed (set flag to 1, 0 otherwise)
    NaVsCaMKIIclamp=0 # if 1, CaMKII Clamp on NaV
)

    cycleLength = 1e3 / freq      # [ms]
    # Na loading parameters ON (set flag to 1, 0 otherwise)
    if CKIIflag == 1
        NaGainFlag = 1  # CaMKII-OE (default 1)
        inashift = -3.25
        alphaCKII = -0.18
    else
        NaGainFlag = 0   # WT (default 0)
        inashift = 0
        alphaCKII = 0
    end

    if CaffeineFlag == 1
        protocol = 2    # "vcRest"
    elseif StrophFlag == 1
        protocol = 3    # "none"
    else
        protocol = 1    # "pace"
    end

    @parameters begin
        R = 8314        # [J/kmol*K]
        Frdy = 96485    # [C/mol]
        Temp = 310      # [K] 310 K (37 C) for BT / 295 K (22 C) for RT
        FoRT = Frdy / R / Temp
        Qpow = (Temp - 310) / 10
        Acell = 20e3        # [um^2] (MOUSE)
        Cmem = Acell * 1e-14  # [F] 200 pF membrane capacitance (MOUSE)
        Fjunc = 3 / 16       # Fraction of junctional current
        Fsl = 1 - Fjunc       # Fraction of sarcoplasmic current
        Fjunc_nak = 9 / (17 * 1.6 + 31)
        Fsl_nak = 1 - Fjunc_nak
        Fjunc_ncx = Fjunc
        Fsl_ncx = Fsl
        Fjunc_CaL = 0.9
        Fsl_CaL = 1 - Fjunc_CaL

    end



    return nothing
end
