using ModelingToolkit

include("ecc/ina_hh.jl")
include("ecc/ica_markov.jl")


function make_ecc_sys(;
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
    loop=0,             # CaMKII-Na-Ca-CaMKII loop closed (set flag to 1, 0 otherwise)
    NaVsCaMKIIclamp=0   # if 1, CaMKII Clamp on NaV
)
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

    # Physical constants and cell physical properties
    @constants begin
        R = 8314        # [J/kmol*K]
        Frdy = 96485    # [C/mol]
        Temp = 310      # [K] 310 K (37 C) for BT / 295 K (22 C) for RT
        FoRT = Frdy / R / Temp
        Qpow = (Temp - 310) / 10
        Acell = 20e3          # [um^2] (MOUSE)
        Cmem = Acell * 1e-14  # [F] 200 pF membrane capacitance (MOUSE)
        Fjunc = 3 / 16        # Fraction of junctional current
        Fsl = 1 - Fjunc       # Fraction of sarcoplasmic current
        Fjunc_nak = 9 / (17 * 1.6 + 31)
        Fsl_nak = 1 - Fjunc_nak
        Fjunc_ncx = Fjunc
        Fsl_ncx = Fsl
        Fjunc_CaL = 0.9
        Fsl_CaL = 1 - Fjunc_CaL
        cellLength = 100                # cell length [um]
        cellRadius = 10.25              # cell radius [um]
        junctionLength = 15e-3          # junc length [um]
        junctionRadius = 160e-3         # junc radius [um]
        distSLcyto = 0.45               # dist. SL to cytosol [um]
        distJuncSL = 0.3                # dist. junc to SL [um] MOUSE
        DcaJuncSL = 1.64e-6             # Dca junc to SL [cm^2/sec]
        DcaSLcyto = 1.22e-6             # Dca SL to cyto [cm^2/sec]
        DnaJuncSL = 1.09e-5             # Dna junc to SL [cm^2/sec]
        DnaSLcyto = 1.79e-5             # Dna SL to cyto [cm^2/sec]
        Vcell = pi*cellRadius^2*cellLength*1e-15 # [L]
        Vmyo = 0.65*Vcell
        Vsr = 0.035*Vcell
        Vsl = 0.02*Vcell
        Vjunc = 0.0539*.01*Vcell
        SAsl = Fsl*Acell                                    # [um^2]  MOUSE
        Njunc = (Fjunc*Acell)/(pi*junctionRadius^2)         # [-]
        SAjunc = Njunc*pi*2*junctionLength*junctionRadius   # [um^2] MOUSE
        J_ca_juncsl = DcaJuncSL*SAjunc/distJuncSL*1e-10     # [L/msec] [m^2/sec] MOUSE
        J_ca_slmyo = DcaSLcyto*SAsl/distSLcyto*1e-10        # MOUSE
        J_na_juncsl = DnaJuncSL*SAjunc/distJuncSL*1e-10     # MOUSE
        J_na_slmyo = DnaSLcyto*SAsl/distSLcyto*1e-10        # MOUSE
    end

    @parameters begin
        cycleLength = 1e3 / freq    # [ms]
        LCCtotDyad = 31.4*.9        # [uM] - Total Dyadic [LCC] - (umol/l dyad)
        RyRtot = 382.6              # [uM] - Total RyR (in Dyad)
        PLBtot=106                  # [uM] - Total [PLB] in cytosolic units, MOUSE
        LCCtotBA = 0.025            # [uM] - [umol/L cytosol]
        PLMtotBA = 48               # [uM] - [umol/L cytosol] MOUSE
        RyRtotBA = 0.135            # [uM] - [umol/L cytosol]
        TnItotBA = 70               # [uM] - [umol/L cytosol]
        IKurtotBA = 0.025           # [uM] - [umol/L cytosol] MOUSE
        PLBtotBA = PLBtot           # [uM] - [umol/L cytosol]
    end



    return nothing
end
