module Module_ECC
    using DifferentialEquations
    using ModelingToolkit

    ## State variables
        # y(1) INa - gate m - HH model
        # y(2) INa - gate h - HH model
        # y(3) INa - gate j - HH model 
        # y(4-7) ICa - HH model (NOT USED)
        # y(8) Itos - gate x (NOT USED)
        # y(9) Itos - gate y (NOT USED)
        # y(10) Itof - gate x
        # y(11) Itof - gate y
        # y(12) IKr - gate x
        # y(13) IKs (NOT USED)
        # y(14) RyR (R)
        # y(15) RyR (O)
        # y(16) RyR (I)
        # y(17) Na Buffer cleft
        # y(18) Na Buffer sl
        # y(19) Ca Buffer myofilament
        # y(20) Ca Buffer TnCH
        # y(21) Mg Buffer TnCH
        # y(22) Ca Buffer CaM (NOT USED)
        # y(23) Ca Buffer Myosin
        # y(24) Mg Buffer Myosin
        # y(25) Ca Buffer SRB
        # y(26) Ca Buffer - low cleft
        # y(27) Ca Buffer - low sl
        # y(28) Ca Buffer - high cleft
        # y(29) Ca Buffer - high sl
        # y(30) Ca-Csqn
        # y(31) [Ca]sr
        # y(32) [Na]cleft
        # y(33) [Na]sl
        # y(34) [Na]i
        # y(35) [K]i (NOT USED)
        # y(36) [Ca]cleft
        # y(37) [Ca]sl
        # y(38) [Ca]i
        # y(39) Vm
        # y(40) Itos - gate r (NOT USED)
        # y(41-42) (NOT USED)
        # y(43-46) used for analysis of influx/efflux via LTCC, PMCA, NCX, ICa,bkg
        # y(47) INa late - gate h
        # y(48-59) INa Markov Model - Placeholder (NOT USED)
        # y(60-65) states ICa Markov Model (m1-cleft)           �
        # y(66-71) states ICa Markov Model (m2-cleft)           �
        # y(72-77) states ICa Markov Model (m1-sl)              �
        # y(78-83) states ICa Markov Model (m2-sl)              �
        # y(84) IKslow1,2 - gate x
        # y(85) IKslow1 - gate y
        # y(86) Iss - gate x
        # y(87) IKslow2 - gate y

    function morotti_eccODE!(du,u,p,t)
    
        (y1, y2, y3, y4, y5, y6, y7, y8, y9, y10,
         y11, y12, y13, y14, y15, y16, y17, y18, y19, y20,
         y21, y22, y23, y24, y25, y26, y27, y28, y29, y30,
         y31, y32, y33, y34, y35, y36, y37, y38, y39, y40,
         y41, y42, y43, y44, y45, y46, y47, y48, y49, y50,
         y51, y52, y53, y54, y55, y56, y57, y58, y59, y60,
         y61, y62, y63, y64, y65, y66, y67, y68, y69, y70,
         y71, y72, y73, y74, y75, y76, y77, y78, y79, y80,
         y81, y82, y83, y84, y85, y86, y87) = u

    ## Parameters
        # P2~P4              CaMKII phosphorylated targets (%)
        # P5~P11, P16, P17   PKA phosphorylated targets (%)
        # P12                Flag for CaMKII-ODE
        # P13, P14           Protocols
        # P15                PKA 
        (cycleLength, LCC_CKp, RyR_CKp, PLB_CKp, LCCa_PKAp,
         LCCb_PKAp, PLB_PKAn, RyR_PKAp, TnI_PKAp, IKs_PKAp, 
         ICFTR_PKAp, CKIIOE, rest, variablePar, Ligtot,
         IKur_PKAp, PLM_PKAp) = p

         CKIIflag = CKIIOE

    ## Flags

        # Set CKIIflag to 1 for CaMKII-OE, 0 otherwise
        # CKIIflag = CKIIOE
        # Set ICa_MarkovFlag to 1 for Markov ICa, 0 otherwise
        ICa_MarkovFlag = 1
        # Set INa_MarkovFlag to 1 for Markov INa, 0 otherwise
        INa_MarkovFlag = 0
        # Set Ito to use either origl params (=0) or Grandi Params (=1)
        ItoFlag = 1

        # Na clamp (set flag to 1, 0 otherwise)
        NaClampFlag = 0
        # PLM KO (set flag to 1, 0 otherwise)
        PLMkoFlag = 0
        # Strophanthidin (dyastolic Na influx) (set flag to 1, 0 otherwise)
        StrophFlag = 0
        # Caffeine-induced Ca transient (set flag to 1, 0 otherwise)
        CaffeineFlag = 0
        # Digitalis (set flag to 1, 0 otherwise)
        DigitalisFlag = 0

        ## Na loading and CaMKII-Na-Ca-CaMKII loop properties

        # Na loading parameters ON (set flag to 1, 0 otherwise)
        NaGainFlag =0 # WT (default 0)
        if CKIIflag == 1
            NaGainFlag = 1 # CaMKII-OE (default 1)
        end

        # CaMKII-Na-Ca-CaMKII loop closed (set flag to 1, 0 otherwise)
        loop = 0 # (default 0)

    ## Model Parameters

        # Constants
        R = 8314        # [J/kmol*K]  
        Frdy = 96485    # [C/mol]  
        Temp = 310      # [K] 310 K (37 C) for BT / 295 K (22 C) for RT
        FoRT = Frdy/R/Temp
        Qpow = (Temp-310)/10

        # Cell geometry
        Acell = 20e3        # [um^2] MOUSE
        #Cmem = 1.3810e-10  # [F] membrane capacitance RABBIT
        Cmem = Acell*1e-14  # [F] 200 pF membrane capacitance MOUSE

        # Fractional currents in compartments
        #Fjunc = 0.11 // Fsl = 1-Fjunc # RABBIT
        #Fjunc = 0.3 // Fsl = 1-Fjunc # MOUSE
        Fjunc = 17/(17+31)*7/17+31/(17+31)*2/31
        Fsl = 1-Fjunc                   # MOUSE Fjunc = 0.1875
        #Fjunc_nak = Fjunc // Fsl_nak = 1-Fjunc_nak # RABBIT
        Fjunc_nak = 1.6*17/(1.6*17+31)*7/17+31/(1.6*17+31)*2/31 
        Fsl_nak = 1-Fjunc_nak           # MOUSE Fjunc = 0.2268
        Fjunc_ncx = Fjunc 
        Fsl_ncx = 1-Fjunc_ncx
        #Fjunc_ncx = Fjunc_nak // Fsl_ncx = 1-Fjunc_ncx
        #Fjunc_ncx = 0.5 // Fsl_ncx = 1-Fjunc_ncx
        #Fjunc_ncx = 0.9 // Fsl_ncx = 1-Fjunc_ncx
        Fjunc_CaL = 0.9 
        Fsl_CaL = 1-Fjunc_CaL

        cellLength = 100                # cell length [um]
        cellRadius = 10.25              # cell radius [um]
        junctionLength = 15e-3          # junc length [um]
        junctionRadius = 160e-3         # junc radius [um]
        distSLcyto = 0.45               # dist. SL to cytosol [um]
        # distJuncSL = 0.5              # dist. junc to SL [um] RABBIT
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
        #SAjunc = 20150*pi*2*junctionLength*junctionRadius # [um^2] RABBIT
        #SAsl = pi*2*cellRadius*cellLength # [um^2] RABBIT
        SAsl = Fsl*Acell                                    # [um^2]  MOUSE
        Njunc = (Fjunc*Acell)/(pi*junctionRadius^2)         # [-]
        SAjunc = Njunc*pi*2*junctionLength*junctionRadius   # [um^2] MOUSE
        #spacing=sqrt(2*Acell/Njunc) # [um] not used -> reduction in distJuncSL

        #J_ca_juncsl = 1/1.2134e12          # [L/msec] [m^2/sec] RABBIT
        #J_ca_slmyo = 1/2.68510e11          # RABBIT
        #J_na_juncsl = 1/(1.6382e12/3*100)  # RABBIT 
        #J_na_slmyo = 1/(1.8308e10/3*100)   # RABBIT 
        J_ca_juncsl = DcaJuncSL*SAjunc/distJuncSL*1e-10     # [L/msec] [m^2/sec] MOUSE
        J_ca_slmyo = DcaSLcyto*SAsl/distSLcyto*1e-10        # MOUSE
        J_na_juncsl = DnaJuncSL*SAjunc/distJuncSL*1e-10     # MOUSE
        J_na_slmyo = DnaSLcyto*SAsl/distSLcyto*1e-10        # MOUSE

        # Fixed ion concentrations     
        Cli = 15    # Intracellular Cl  [mM]
        Clo = 150   # Extracellular Cl  [mM]
        Ko = 5.4    # Extracellular K   [mM]
        Nao = 140   # Extracellular Na  [mM]
        Cao = 1     # Extracellular Ca  [mM] MOUSE # 1.8 mM in RABBIT
        Mgi = 1     # Intracellular Mg  [mM]

    ## Na transport parameters

        GNa = 10 # [mS/uF] changed from rabbit (16)
        GNaB = 4.5 * 0.297e-3 # [mS/uF] changed from rabbit

        IbarNaK = 5 # [uA/uF] changed from rabbit (1.90719)
        if NaGainFlag == 1
            GNaB = GNaB * 4
            IbarNaK = IbarNaK*0.9
        end
        if DigitalisFlag == 1
            IbarNaK = IbarNaK*0.5 # 50# block
        end
        KmNaip = 19 # [mM] changed from rabbit (11)
        KmKo = 1.5 # [mM]
        Q10NaK = 1.63
        Q10KmNai = 1.39
        if PLMkoFlag == 1
            PLM_PKAp = 1
            GNaB = GNaB * 48 / 20
            IbarNaK = IbarNaK * 0.8
        end
        if StrophFlag == 1
            IbarNaK = 0
        end

        # INa Markov Model parameters
        GNa2 = 10.64  #23      # [mS/uF]

        P1a1 = 3.802
        P2a1 = 0.1027
        P3a1 = 2.5
        P4a1 = 17
        P5a1 = 0.20
        P6a1 = 150
        P4a2 = 15
        P5a2 = 0.23
        P4a3 = 12
        P5a3 = 0.25
        P1b1 = 0.1917
        P2b1 = 20.3
        P1b2 = 0.2
        P2b2 = 2.5
        P1b3 = 0.22
        P2b3 = 7.5
        P1a4 = 0.188495
        P2a4 = 16.6
        P3a4 = 0.393956
        P4a4 = 7
        P1a5 = 7e-7
        P2a5 = 7.2 # TG 7.8
        P1b5 = 0.0044 # TG 0.00467
        P2b5 = 2e-5
        P1a6 = 100
        P1b6 = 8.9554e-7 # TG 6.5449e-7
        P2b6 = 11.3944
        P1a7 = 0.487e-4 # TG 0.3377e-3
        P2a7 = 23.2696
        P1b7 = 0.2868e-3 # TG 1.868e-4
        P2b7 = 35.9898
        P1a8 = 0.1e-7 # TG 6.5e-6
        P1b8 = 9.8e-3 # TG 3.8e-3
        if CKIIflag == 1 # MOUSE - CaMKII-OE
            P2a5 = 7.8
            P1b5 = 0.00467
            P1b6 = 6.5449e-7
            P1a7 = 0.3377e-3
            P1b7 = 1.868e-4
            P1a8 = 6.5e-6
            P1b8 = 3.8e-3
        end

    ## K currents parameters

        Kcoeff = 1 # K current modulation
                
        pNaK = 0.01833
        GtoSlow = 0 # [mS/uF] changed from rabbit (0.06): NO ItoSlow in MOUSE
        GtoFast = 0.44 # [mS/uF] changed from rabbit (0.02)
        if CKIIflag == 1 # MOUSE
            GtoFast = GtoFast * 2/3 # chronic CaMKII-OE effect
        end
        #Gkur = 0.3  # [mS/uF] only in MOUSE
        Gkur1 = 1.1 * 0.16 # fast
        Gkur2 = 0.14 # slow
        Gss = 0.15 # [mS/uF] only in MOUSE
        gkr = 0.03 * sqrt(Ko/5.4)
        gkp = 0.001

        # Cl current parameters
        GClCa = 0.109625 # [mS/uF]
        GClB = 9e-3 # [mS/uF]
        KdClCa = 100e-3 # [mM]

    ## LTCC parameters

        K_Ica = 1.65 # MOUSE
        pNa = K_Ica * 1.5e-8 # [cm/sec]
        pCa = K_Ica * 5.4e-4 # [cm/sec] - Ca permeability
        pK = K_Ica * 2.7e-7 # [cm/sec]
        KmCa = 0.6e-3 # [mM]
        Q10CaL = 1.8

    ## Ca transport parameters

        IbarNCX = 1 # [uA/uF] changed from rabbit (9)
        if CKIIflag == 1
            IbarNCX = 1.5 * IbarNCX
        end
        KmCai = 3.59e-3         # [mM]
        KmCao = 1.3             # [mM]
        KmNai = 12.29           # [mM]
        KmNao = 87.5            # [mM]
        ksat = 0.27             # [none]  
        nu = 0.35               # [none]
        Kdact = 1/2 * 0.256e-3  # [mM] changed from rabbit
        Q10NCX = 1.57           # [none]
        IbarSLCaP = 0.0673      # [uA/uF]
        KmPCa = 0.5e-3          # [mM]
        GCaB = 3 * 2.513e-4     # [uA/uF] changed from rabbit (2.513e-4)
        Q10SLCaP = 2.35         # [none]

        # SR flux parameters
        Q10SRCaP = 2.6                  # [none]
        Vmax_SRCaP = 1.15*1.15*2.86e-4  # [mM/msec] (mmol/L cytosol/msec) changed
        Kmf = 0.3e-3                    # [mM] changed from rabbit (0.246e-3) # from Yang-Saucerman
        Kmr = 2.1                       # [mM]L cytosol changed from rabbit (1.7) # from Yang-Saucerman
        hillSRCaP = 1.787               # [mM]
        ks = 25                         # [1/ms]      
        koCa = 10                       # [mM^-2 1/ms]
        kom = 0.06                      # [1/ms]     
        kiCa = 0.5                      # [1/mM/ms]
        kim = 0.005                     # [1/ms]
        ec50SR = 0.45 + 0.05            # [mM] changed from rabbit (0.45)

        if CaffeineFlag == 1
            koCa = koCa * 7.5
            GCaB = 0
            Vmax_SRCaP = 0
        end

    ## Buffering parameters

        Bmax_Naj = 7.561        # [mM] 
        Bmax_Nasl = 1.65        # [mM]
        koff_na = 1e-3          # [1/ms]
        kon_na = 0.1e-3         # [1/mM/ms]
        Bmax_TnClow = 70e-3     # [mM]                      # TnC low affinity
        koff_tncl = 19.6e-3     # [1/ms] 
        kon_tncl = 32.7         # [1/mM/ms]
        Bmax_TnChigh = 140e-3   # [mM]                      # TnC high affinity 
        koff_tnchca = 0.032e-3  # [1/ms] 
        kon_tnchca = 2.37       # [1/mM/ms]
        koff_tnchmg = 3.33e-3   # [1/ms] 
        kon_tnchmg = 3e-3       # [1/mM/ms]
        Bmax_myosin = 140e-3    # [mM]                      # Myosin buffering
        koff_myoca = 0.46e-3    # [1/ms]
        kon_myoca = 13.8        # [1/mM/ms]
        koff_myomg = 0.057e-3   # [1/ms]
        kon_myomg = 0.0157      # [1/mM/ms]
        Bmax_SR = 19 * 0.9e-3   # [mM] 
        koff_sr = 60e-3         # [1/ms]
        kon_sr = 100            # [1/mM/ms]
        Bmax_SLlowsl = 37.38e-3 * Vmyo/Vsl          # [mM]    # SL buffering
        Bmax_SLlowj = 4.62e-3 * Vmyo/Vjunc * 0.1    # [mM]    
        koff_sll = 1300e-3      # [1/ms]
        kon_sll = 100           # [1/mM/ms]
        Bmax_SLhighsl = 13.35e-3 * Vmyo/Vsl         # [mM] 
        Bmax_SLhighj = 1.65e-3 * Vmyo/Vjunc * 0.1   # [mM] 
        koff_slh = 30e-3        # [1/ms]
        kon_slh = 100           # [1/mM/ms]
        Bmax_Csqn = 2.7         # 140e-3*Vmyo/Vsr  [mM] 
        koff_csqn = 65          # [1/ms] 
        kon_csqn = 100          # [1/mM/ms] 

        # PKA-dependent phosphoregulation of TnI (increases Kd of TnC)
        fracTnIpo = 0.062698  # Derived quantity (TnI_PKAp(baseline)/TnItot)
        #fPKA_TnI = (1.45-0.45*(1-TnI_PKAp)/(1-fracTnIpo))  # Max effect +45#
        fPKA_TnI = (1.61 - 0.61 * (1 - TnI_PKAp) / (1 - fracTnIpo)) # Max effect +61#
        koff_tncl = koff_tncl * fPKA_TnI


    ## Model Parameters
        # Nernst Potentials
        @variables ena_junc(t) ena_junc ~ (1/FoRT)*log(Nao/y32)      # [mV]
        @variables ena_sl(t) ena_sl ~ (1/FoRT)*log(Nao/y33)          # [mV]
        @variables ek(t) ek ~ (1/FoRT)*log(Ko/y35)                   # [mV]
        @variables eca_junc(t) eca_junc ~ (1/FoRT/2)*log(Cao/y36)    # [mV]
        @variables eca_sl(t) eca_sl ~ (1/FoRT/2)*log(Cao/y37)        # [mV]
        @variables ecl(t) ecl ~ (1/FoRT)*log(Cli/Clo)                # [mV]

    ## I_Na: Fast Na Current

        # Max INa alterations with CaMKII hyperactivity as in Hund & Rudy 2008
        if CKIIflag == 1 # acute effects
            inashift = -3.25
            alphaCKII = -0.18
            #deltGbarNal_CKII = 2  # RABBIT
            if NaGainFlag == 1
                deltGbarNal_CKII = 3  # MOUSE
            else
                deltGbarNal_CKII = 0  # no Na Gain in OE
            end
        else
            inashift = 0
            alphaCKII = 0
            deltGbarNal_CKII = 0
        end
        
        if loop == 1
            RyRp_WT_mean = 0.2101 
            RyRp_OE_mean = 0.7387 # Derived (1 Hz, no loop)
            RyRp_OEloop_min = 0.7033 # Derived (1 Hz, OE loop)
            delta_loop = (3/(RyRp_OE_mean-RyRp_WT_mean)) * RyR_CKp - (3/(RyRp_OE_mean-RyRp_WT_mean)) * RyRp_WT_mean
            NaVsCaMKIIclamp = 0 # if 1, CaMKII Clamp on NaV
            if NaVsCaMKIIclamp == 1
                delta_loop = (3/(RyRp_OE_mean-RyRp_WT_mean)) * RyRp_OEloop_min - (3/(RyRp_OE_mean-RyRp_WT_mean)) * RyRp_WT_mean
            end
            GNaB=(4.5) * 0.297e-3 * (1+delta_loop)
            if CKIIflag == 1 # OE
                if NaGainFlag == 1 # acute
                    deltGbarNal_CKII=delta_loop
                else
                    deltGbarNal_CKII = 0
                end
            else # WT
                deltGbarNal_CKII = 0
            end          
        end
        
        @variables am(t) am ~ 0.32 * (y39 + 47.13)/(1-exp(-0.1*(y39 + 47.13)))
        @variables bm(t) bm ~ 0.08 * exp(-(y39)/11)
        if (y39-inashift) >= -40
            ah = 0 
            aj = 0
            #bh = 1/(0.13*(1+exp(-((y(39)-inashift)+10.66)/11.1))) # RABBIT
            @variables bh(t) bh ~ 0.66 * 1/(0.13 * (1 + exp(-((y39-inashift) + 10.66)/11.1))) # MOUSE
            @variables bj(t) bj ~ 0.3 * exp(-2.535e-7 * (y39-inashift))/(1+exp(-0.1 * ((y39-inashift) + 32)))
        else
            @variables ah(t) ah ~ 0.135 * exp((80+(y39-inashift))/(-6.8))
            #bh = 3.56*exp(0.079*(y(39)-inashift))+3.1e5*exp(0.35*(y(39)-inashift)) # RABBIT
            @variables bh(t) bh ~ 1.1 * 3.56 * exp(0.079 * (y39-inashift-2)) + 3.1e5 * exp(0.35 * (y39-inashift-2)) # MOUSE
            # Including alteration to aj as in Hund and Rudy 2008
            @variables aj(t) aj ~ (1+alphaCKII)*((-1.2714e5*exp(0.2444*(y39-inashift))-3.474e-5*exp(-0.04391*(y39-inashift)))*((y39-inashift)+37.78)/(1+exp(0.311*((y39-inashift)+79.23))))
            @variables bj(t) bj ~ 0.1212 * exp(-0.01052*(y39-inashift))/(1+exp(-0.1378*((y39-inashift)+40.14)))
        end
        
        du[1] = am*(1-y1)-bm*y1
        du[2] = ah*(1-y2)-bh*y2
        du[3] = aj*(1-y3)-bj*y3

        @variables I_Na_junc1(t) I_Na_junc1 ~ Fjunc*GNa*(y1^3)*y2*y3*(y39-ena_junc)
        @variables I_Na_sl1(t) I_Na_sl1 ~ Fsl*GNa*(y1^3)*y2*y3*(y39-ena_sl)

    ## I_Na,L: Late INa current (as in Hund & Rudy 2008)

        GbarNal = 0.0065*(1+deltGbarNal_CKII)*2 # deltGbar assigned in 'Fast INa' section
        
        # h-gate (note: m-gate is same as INa m-gate -> using y(1) for this)
        @variables hlss(t) hlss ~ 1/(1+exp((y39+91)/6.1))
        tauhl = 600 # ms
        @variables I_Nalj(t) I_Nalj ~ Fjunc*GbarNal*(y1^3)*y47*(y39-ena_junc)
        @variables I_Nalsl(t) I_Nalsl ~ Fsl*GbarNal*(y1^3)*y47*(y39-ena_sl)
        #I_Nal = I_Nalj+I_Nalsl

    ## I_Na: alternative Markov Model - unused

        if INa_MarkovFlag == 1
            # State variables
            @variables CNa2(t) CNa2 ~ y48
            @variables CNa1(t) CNa1 ~ y49
            @variables ONa(t) ONa ~ y50
            @variables IFNa(t) IFNa ~ y51
            @variables I1Na(t) I1Na ~ y52
            @variables CNa3(t) CNa3 ~ y53
            @variables ICNa2(t) ICNa2 ~ y54
            @variables ICNa3(t) ICNa3 ~ y55
            @variables LONa(t) LONa ~ y56
            @variables LCNa1(t) LCNa1 ~ y57
            @variables LCNa2(t) LCNa2 ~ y58
            @variables LCNa3(t) LCNa3 ~ y59
            @variables I2Na(t) I2Na ~ (1-(ONa+CNa1+CNa2+CNa3+IFNa+I1Na+ICNa2+ICNa3+LONa+LCNa1+LCNa2+LCNa3))
            # Transition rates
            @variables alphaNa1(t) alphaNa1 ~ P1a1/(P2a1*exp(-(y39+P3a1)/P4a1)+P5a1*exp(-(y39+P3a1)/P6a1))
            @variables alphaNa2(t) alphaNa2 ~ P1a1/(P2a1*exp(-(y39+P3a1)/P4a2)+P5a2*exp(-(y39+P3a1)/P6a1))
            @variables alphaNa3(t) alphaNa3 ~ P1a1/(P2a1*exp(-(y39+P3a1)/P4a3)+P5a3*exp(-(y39+P3a1)/P6a1))
            @variables betaNa1(t) betaNa1 ~ P1b1*exp(-(y39+P3a1)/P2b1) # shift
            @variables betaNa2(t) betaNa2 ~ P1b2*exp(-(y39-P2b2)/P2b1)
            @variables betaNa3(t) betaNa3 ~ P1b3*exp(-(y39-P2b3)/P2b1)
            @variables alphaNa4(t) alphaNa4 ~ 1/(P1a4*exp(-(y39+P4a4)/P2a4)+P3a4)
            @variables alphaNa5(t) alphaNa5 ~ P1a5*exp(-(y39+P4a4)/P2a5)
            @variables betaNa5(t) betaNa5 ~ (P1b5+P2b5*(y39+P4a4))
            @variables betaNa6(t) betaNa6 ~ P1b6*exp(-y39/P2b6)
            @variables alphaNa7(t) alphaNa7 ~ P1a7*exp(y39/P2a7)
            @variables betaNa7(t) betaNa7 ~ P1b7*exp(-y39/P2b7)
            alphaNa8 = P1a8
            betaNa8 = P1b8
            @variables betaNa4(t) betaNa4 ~ (alphaNa3*alphaNa4*alphaNa5)/(betaNa3*betaNa5)
            @variables alphaNa6(t) alphaNa6 ~ alphaNa4/P1a6
            # ODEs
            dCNa3  = betaNa8*LCNa3+betaNa1*CNa2+alphaNa5*ICNa3-(alphaNa1+betaNa5+alphaNa8)*CNa3
            dCNa2  = betaNa8*LCNa2+alphaNa1*CNa3+betaNa2*CNa1+alphaNa5*ICNa2-(betaNa1+alphaNa2+betaNa5+alphaNa8)*CNa2
            dCNa1  = betaNa8*LCNa1+alphaNa2*CNa2+betaNa3*ONa+alphaNa5*IFNa-(betaNa2+alphaNa3+betaNa5+alphaNa8)*CNa1
            dONa   = betaNa8*LONa+alphaNa3*CNa1+betaNa4*IFNa-(betaNa3+alphaNa4+alphaNa8)*ONa
            dIFNa  = alphaNa4*ONa+betaNa5*CNa1+betaNa6*I1Na+alphaNa2*ICNa2-(betaNa4+alphaNa5+alphaNa6+betaNa2)*IFNa
            dI1Na  = alphaNa6*IFNa+betaNa7*I2Na-(betaNa6+alphaNa7)*I1Na
            dICNa2 = alphaNa1*ICNa3+betaNa2*IFNa+betaNa5*CNa2-(betaNa1+alphaNa2+alphaNa5)*ICNa2
            dICNa3 = betaNa1*ICNa2+betaNa5*CNa3-(alphaNa1+alphaNa5)*ICNa3
            dLONa  = alphaNa3*LCNa1+alphaNa8*ONa-(betaNa8+betaNa3)*LONa
            dLCNa1 = alphaNa8*CNa1+alphaNa2*LCNa2+betaNa3*LONa-(betaNa8+betaNa2+alphaNa3)*LCNa1
            dLCNa2 = betaNa2*LCNa1+alphaNa8*CNa2+alphaNa1*LCNa3-(betaNa8+betaNa1+alphaNa2)*LCNa2
            dLCNa3 = alphaNa8*CNa3+betaNa1*LCNa2-(betaNa8+alphaNa1)*LCNa3
            du[48] = dCNa2
            du[49] = dCNa1
            du[50] = dONa
            du[51] = dIFNa
            du[52] = dI1Na
            du[53] = dCNa3
            du[54] = dICNa2
            du[55] = dICNa3
            du[56] = dLONa
            du[57] = dLCNa1
            du[58] = dLCNa2
            du[59] = dLCNa3
        else # If not using INa Markov, set ODEs to zero to speed simulations
            du[48] = 0
            du[49] = 0
            du[50] = 0
            du[51] = 0
            du[52] = 0
            du[53] = 0
            du[54] = 0
            du[55] = 0
            du[56] = 0
            du[57] = 0
            du[58] = 0
            du[59] = 0
        
        end
        @variables I_Na_junc2(t) I_Na_junc2 ~ Fjunc*GNa2*(y50+y56)*(y39-ena_junc) # junc current
        @variables I_Na_sl2(t) I_Na_sl2 ~ Fsl*GNa2*(y50+y56)*(y39-ena_sl) # sl current

    ## I_Na: compute total current (fast and late components)

        if INa_MarkovFlag == 1
            @variables I_Na_junc(t) I_Na_junc ~ I_Na_junc2
            @variables I_Na_sl(t) I_Na_sl ~ I_Na_sl2
        else
            @variables I_Na_junc(t) I_Na_junc ~ I_Na_junc1 + I_Nalj
            @variables I_Na_sl(t) I_Na_sl ~ I_Na_sl1 + I_Nalsl
        end
        
        # I_Na = I_Na_junc + I_Na_sl
        # I_Na_store(tStep) = I_Na

    ## I_nabk: Na Background Current

        @variables I_nabk_junc(t) I_nabk_junc ~ Fjunc*GNaB*(y39-ena_junc)
        @variables I_nabk_sl(t) I_nabk_sl ~ Fsl*GNaB*(y39-ena_sl)
        @variables I_nabk(t) I_nabk ~ I_nabk_junc+I_nabk_sl
        
        # global I_Nabk_store
        # I_Nabk_store(tStep) = I_nabk

    ## I_nak: Na/K Pump Current

        sigma = (exp(Nao/67.3)-1)/7
        @variables fnak(t) fnak ~ 1/(1+0.1245*exp(-0.1*y39*FoRT)+0.0365*sigma*exp(-y39*FoRT))
        #KmNaip = 14 # PKA effect - 1 mM (Despa 2005)
        fracPKA_PLMo = 0.116738 # Derived quantity (PLM_PKAp(baseline)/PLMtot)
        fracPKA_PLMiso = 0.859251 # Derived quantity (PLM_PKAp(ISO)/PLMtot)
        #kPKA_PLM=(KmNaip-14)/(fracPKA_PLMiso/fracPKA_PLMo-1) # PLM_PKAp ISO
        kPKA_PLM=KmNaip*(1-0.7019)/(fracPKA_PLMiso/fracPKA_PLMo-1) # PLM_PKAp ISO
        KmNaip_PKA=-kPKA_PLM+kPKA_PLM*(PLM_PKAp/fracPKA_PLMo)
        KmNaip = KmNaip-KmNaip_PKA
        
        @variables I_nak_junc(t) I_nak_junc ~ Fjunc_nak*IbarNaK*fnak*Ko /(1+(KmNaip/y32)^4) /(Ko+KmKo)
        @variables I_nak_sl(t) I_nak_sl ~ Fsl_nak*IbarNaK*fnak*Ko /(1+(KmNaip/y33)^4) /(Ko+KmKo)
        @variables I_nak(t) I_nak ~ I_nak_junc+I_nak_sl

        # global I_NaK_store
        # I_NaK_store(tStep) = I_nak

    ## I_kur - IK,slow

        @variables xurss(t) xurss ~ 1/(1+exp(-(y39+15)/14))
        @variables yurss(t) yurss ~ 1/(1+exp((y39+48)/6.2))
        
        tauxur = 0.95+0.05*exp(-0.08*y39)
        du[84] = (xurss-y84)/tauxur         # IKslow1
        tauxur2 = 1+7/(1+exp(-(y(39)+45)/8))+20*exp(-((y(39)+35)/10)^2) #8+20*exp(-((y(39)+35)/10)^2)
        du[13] = (xurss-y13)/tauxur2        # IKslow2
        #tauyur = 400+900*exp(-((y(39)+55)/16)^2)+150/(1+exp(-(y(39)+60)/8))
        #ydot(85) = (yurss-y(85))/tauyur
        tauyur1 = 400+900*exp(-((y(39)+55)/16)^2)-250/(1+exp(-(y(39)+60)/8)) # fast
        du[85] = (yurss-y85)/tauyur1 
        tauyur2 = 400+900*exp(-((y(39)+55)/16)^2)+550/(1+exp(-(y(39)+60)/8)) # slow
        du[87] = (yurss-y87)/tauyur2
        
        # PKA-dependent phosphoregulation of Ik,slow1 (increases Gkur1)
        fracIKurp0 = 0.437635  # Derived quantity (IKur_PKAp(baseline)/IKurtot)
        fracIKurpISO = 0.718207 # Derived quantity (IKur_PKAp(ISO)/IKurtot)
        a_Kur = (1.20-1)/(fracIKurpISO/fracIKurp0-1)
        fracIKuravail = (1-a_Kur)+a_Kur*(IKur_PKAp/fracIKurp0) # +20# with 0.1 uM ISO
        
        #I_kur = Kcoeff*fracIKuravail*Gkur*y(84)*y(85)*(y(39)-ek)
        @variables I_kur1(t) I_kur1 ~ Kcoeff*fracIKuravail*Gkur1*y84*y85*(y39-ek) # IKslow1
        #I_kur2 = Kcoeff*fracIKuravail*Gkur2*y(84)*y(87)*(y(39)-ek) # IKslow2
        @variables I_kur2(t) I_kur2 ~ Kcoeff*Gkur2*y84*y87*(y39-ek) # IKslow2 # no PKA effect
        @variables I_kur(t) I_kur ~ I_kur1 + I_kur2
        
        # global I_kur1_store, I_kur2_store
        # I_kur1_store(tStep) = I_kur1
        # I_kur2_store(tStep) = I_kur2

    ## I_ss

        @variables xssss(t) xssss ~ xurss # = 1/(1+exp(-(y(39)+15)/14))
        #tauxss = 14+0.8*exp(-0.08*y(39))
        tauxss = 70*exp(-((y39+43)/30)^2)+14 # Iss
        du[86] = (xssss-y86)/tauxss
        @variables I_ss(t) I_ss ~ Kcoeff*Gss*y86*(y39-ek) #store
        
        # global I_ss_store
        # I_ss_store(tStep) = I_ss

    ## I_kr: Rapidly Activating K Current

        @variables xrss(t) xrss ~ 1/(1+exp(-(y39+50)/7.5))
        @variables tauxr(t) tauxr ~ 1/(1.38e-3*(y39+7)/(1-exp(-0.123*(y39+7)))+6.1e-4*(y39+10)/(exp(0.145*(y39+10))-1))
        du[12] = (xrss-y12)/tauxr
        @variables rkr(t) rkr ~ 1/(1+exp((y39+33)/22.4))
        @variables I_kr(t) I_kr ~ Kcoeff*gkr*y12*rkr*(y39-ek)
        
        # global I_kr_store
        # I_kr_store(tStep) = I_kr

    ## I_ks: Slowly Activating K Current

        # Phosphoregulation of IKs by PKA parameters
        # fracIKspo = 0.07344   # Derived quantity (IKs_PKAp(baseline)/IKstot)
        # fracIKsavail = (0.2*(IKs_PKAp/fracIKspo)+0.8)
        # Xs05 = 1.5*(2.0 - IKs_PKAp/fracIKspo)
        # 
        # pcaks_junc = -log10(y(36))+3.0
        # pcaks_sl = -log10(y(37))+3.0
        # gks_junc = fracIKsavail*0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_junc)/0.6))) # Now regulated by PKA
        # gks_sl = fracIKsavail*0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_sl)/0.6)))     # Now regulated by PKA
        # eks = (1/FoRT)*log((Ko+pNaK*Nao)/(y(35)+pNaK*y(34)))
        # # xsss = 1/(1+exp(-(y(39)-1.5)/16.7)) # Original version
        # xsss = 1/(1+exp(-(y(39)-Xs05)/16.7))  # Now regulated by PKA
        # tauxs = 1/(7.19e-5*(y(39)+30)/(1-exp(-0.148*(y(39)+30)))+1.31e-4*(y(39)+30)/(exp(0.0687*(y(39)+30))-1)) 
        #ydot(13) = 0 #(xsss-y(13))/tauxs
        # state variable y(13) is now used for IKslow2 activation
        I_ks_junc = 0   #Fjunc*gks_junc*y(13)^2*(y(39)-eks) # No IKs in mouse
        I_ks_sl = 0     #Fsl*gks_sl*y(13)^2*(y(39)-eks) # No IKs in mouse
        I_ks = I_ks_junc+I_ks_sl
        
        # global IKs_store
        # IKs_store(tStep) = I_ks

    ## I_kp: Plateau K current

        @variables kp_kp(t) kp_kp ~ 1/(1+exp(7.488-y39/5.98))
        @variables I_kp_junc(t) I_kp_junc ~ Kcoeff*Fjunc*gkp*kp_kp*(y39-ek)
        @variables I_kp_sl(t) I_kp_sl ~ Kcoeff*Fsl*gkp*kp_kp*(y39-ek)
        @variables I_kp(t) I_kp ~ I_kp_junc+I_kp_sl

    ## I_to: Transient Outward K Current (slow and fast components)

        # Itos (ABSENT IN MOUSE)
        @variables xtoss(t) xtoss = 1/(1+exp(-(y39+3.0)/13))
        @variables ytoss(t) ytoss = 1/(1+exp((y39+48)/5))
        @variables rtoss(t) rtoss = 1/(1+exp((y39+33.5)/10)) # Rto not used in MOUSE model
        
        @variables tauxtos(t) tauxtos ~ 0.08+0.7*exp(-((y39+25)/30)^2)
        if ItoFlag == 0 # (not used)
            # Shannon Versions
            #tauytos = 3e3/(1+exp((y(39)+60.0)/10))+30
            @variables taurtos(t) taurtos ~ 2.8e3/(1+exp((y39+60.0)/10))+220 # no Rto
            @variables tauytos(t) tauytos ~ 100+400/(1+exp((y39+25)/5))
        elseif ItoFlag == 1 && CKIIflag == 0 # WT
            # Grandi Versions
            Py = 182
            Pr1 = 8085 
            Pr2 = 313            # Normal
            #tauytos = Py/(1+exp((y(39)+33.5)/10))+1
            @variables taurtos(t) taurtos ~ Pr1/(1+exp((y39+33.5)/10))+Pr2 # no Rto
            @variables tauytos(t) tauytos ~ 100+400/(1+exp((y39+25)/5))
        elseif ItoFlag == 1 && CKIIflag == 1            # CaMKII-OE acute effect
            Py = 15
            Pr1 = 3600 
            Pr2 = 500
            #GtoSlow = GtoSlow*1.5  # Rabbit
            #tauytos = Py/(1+exp((y(39)+33.5)/10))+1
            @variables taurtos(t) taurtos ~ Pr1/(1+exp((y39+33.5)/10))+Pr2 # no Rto
            #tauytos = 35+400/(1+exp((y(39)+25)/5))    # MOUSE tau rec -73#
            @variables tauytos(t) tauytos ~ 100+35/(1+exp((y39+25)/5))     # MOUSE tau rec -73#
        end
        
        du[8] = (xtoss-y8)/tauxtos
        du[9] = (ytoss-y9)/tauytos
        du[40] = (rtoss-y40)/taurtos        # Rto not used in MOUSE model
        #I_tos = GtoSlow*y(8)*(y(9)+0.5*y(40))*(y(39)-ek) # [uA/uF]
        @variables I_tos(t) I_tos ~ 0*Kcoeff*GtoSlow*y8*y9*(y39-ek) # N0 Itos in MOUSE # [uA/uF]
        
        # Itof
        @variables xtofs(t) xtofs ~ 1/(1+exp(-(y39+3.0)/13)) # = xtoss
        @variables ytofs(t) ytofs ~ 1/(1+exp((y39+48)/5)) # = ytoss
        
        #tauxtof = 3.5*exp(-y(39)*y(39)/30/30)+1.5 # Original
        #tauxtof = 3.5*exp(-((y(39)+3)/30)^2)+1.5 # Version in Grandi Code (does not change AP shape)
        @variables tauxtof(t) tauxtof ~ 0.08+0.7*exp(-((y39+25)/30)^2) # = tauxtos 
        @variables tauytof(t) tauytof ~ 10+32*exp(-((y39+55)/16)^2)+8/(1+exp(-(y39+60)/8))
        if CKIIflag == 1 # MOUSE (CaMKII-OE acute effect) # tau rec -38#
            @variables tauytof(t) tauytof ~ 5+32*exp(-((y39+55)/16)^2)+12/(1+exp(-(y39+60)/8))
        end
        
        du[10] = (xtofs-y10)/tauxtof
        du[11] = (ytofs-y11)/tauytof
        @variables I_tof(t) I_tof ~ Kcoeff*GtoFast*y10*y11*(y39-ek)
        
        @variables I_to(t) I_to ~ I_tos + I_tof
        
        # global I_to_store
        # I_to_store(1,tStep) = I_to     # Total I_to
        # I_to_store(2,tStep) = I_tof    # Fast Component
        # I_to_store(3,tStep) = I_tos    # Slow component

    ## I_k1: Time-Independent K Current (I_ki)

        @variables aki(t) aki ~ 1.02/(1+exp(0.2385*(y39-ek-59.215)))
        @variables bki(t) bki ~ (0.49124*exp(0.08032*(y39+5.476-ek))+exp(0.06175*(y39-ek-594.31)))/(1 + exp(-0.5143*(y39-ek+4.753)))
        @variables kiss(t) kiss ~ aki/(aki+bki)
        #I_ki = 0.9*sqrt(Ko/5.4)*kiss*(y(39)-ek) # RABBIT
        if CKIIflag == 1 # MOUSE (chronic CaMKII-OE effect)
            @variables I_ki(t) I_ki ~ 1/2*0.3*sqrt(Ko/5.4)*kiss*(y39-ek)*Kcoeff
        else
            @variables I_ki(t) I_ki ~ 0.3*sqrt(Ko/5.4)*kiss*(y39-ek)*Kcoeff
        end
        
        # global I_K1_store
        # I_K1_store(tStep) = I_ki

    ## I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current

        @variables I_ClCa_junc(t) I_ClCa_junc ~ Fjunc*GClCa/(1+KdClCa/y36)*(y39-ecl)
        @variables I_ClCa_sl(t) I_ClCa_sl ~ Fsl*GClCa/(1+KdClCa/y37)*(y39-ecl)
        @variables I_ClCa(t) I_ClCa ~ I_ClCa_junc+I_ClCa_sl
        @variables I_Clbk(t) I_Clbk ~ GClB*(y39-ecl)

    ## Original H-H formulation for LCC - unused if ICa_MarkovFlag = 1

        @variables dss(t) dss ~ 1/(1+exp(-(y39+14.5)/6.0))
        @variables taud(t) taud ~ dss*(1-exp(-(y39+14.5)/6.0))/(0.035*(y39+14.5))
        @variables fss(t) fss ~ 1/(1+exp((y39+35.06)/3.6))+0.6/(1+exp((50-y39)/20))
        @variables tauf(t) tauf ~ 1/(0.0197*exp( -(0.0337*(y39+14.5))^2 )+0.02)
        
        # ydot(4) = (dss-y(4))/taud
        # ydot(5) = (fss-y(5))/(tauf)
        # ydot(6) = (1.7)*y(36)*(1-y(6))-11.9e-3*y(6) # fCa_junc  
        # ydot(7) = 1.7*y(37)*(1-y(7))-11.9e-3*y(7) # fCa_sl
        du[4] = 0
        du[5] = 0
        du[6] = 0
        du[7] = 0
        
        @variables ibarca_j(t) ibarca_j ~ pCa*4*(y39*Frdy*FoRT) * (0.341*y36*exp(2*y39*FoRT)-0.341*Cao) /(exp(2*y39*FoRT)-1)
        @variables ibarca_sl(t) ibarca_sl ~ pCa*4*(y39*Frdy*FoRT) * (0.341*y37*exp(2*y39*FoRT)-0.341*Cao) /(exp(2*y39*FoRT)-1)
        @variables ibark(t) ibark ~ pK*(y39*Frdy*FoRT)*(0.75*y35*exp(y39*FoRT)-0.75*Ko) /(exp(y39*FoRT)-1)
        @variables ibarna_j(t) ibarna_j ~ pNa*(y39*Frdy*FoRT) *(0.75*y32*exp(y39*FoRT)-0.75*Nao)  /(exp(y39*FoRT)-1)
        @variables ibarna_sl(t) ibarna_sl ~ pNa*(y39*Frdy*FoRT) *(0.75*y33*exp(y39*FoRT)-0.75*Nao)  /(exp(y39*FoRT)-1)
        
        @variables I_Ca_junc1(t) I_Ca_junc1 ~ (Fjunc_CaL*ibarca_j*y(4)*y(5)*(1-y(6))*Q10CaL^Qpow)*0.45
        @variables I_Ca_sl1(t) I_Ca_sl1 ~ (Fsl_CaL*ibarca_sl*y(4)*y(5)*(1-y(7))*Q10CaL^Qpow)*0.45
        @variables I_CaK1(t) I_CaK1 ~ (ibark*y(4)*y(5)*(Fjunc_CaL*(1-y(6))+Fsl_CaL*(1-y(7)))*Q10CaL^Qpow)*0.45
        @variables I_CaNa_junc1(t) I_CaNa_junc1 ~ (Fjunc_CaL*ibarna_j*y(4)*y(5)*(1-y(6))*Q10CaL^Qpow)*0.45
        @variables I_CaNa_sl1(t) I_CaNa_sl1 ~ (Fsl_CaL*ibarna_sl*y(4)*y(5)*(1-y(7))*Q10CaL^Qpow)*0.45

    ## LCC MARKOV MODEL - based on Mahajan et al. (2008)

        # This portion contains Markov state transitions for four channel types:
        # 'mode 1' channels in the junction and sl and 'mode 2' channels in the
        # same two compartments. Markov state transitions are computed for each
        # channel type independently - total currents are the sum of the two
        # channel types in each compartment (i.e. ICatot_junc = ICa_mode1_junc +
        # ICa_mode2_junc). Ca-dependent transition rates differ between juncitonal
        # and sl channels, whereas closing rate (r2) is adjusted to define mode1
        # vs. mode2 channels. Parameters determined through microscopic
        # reversibility are redefined to preserve constraint.
        
        # CaMKII shifts distribution of junctional and subsarcolemmal channels to 
        # either mode 2 at the expense of mode 1 channels (i.e. 10# mode 2 results 
        # in 90# mode 1).
        
        # PKA alters overall availability of channels (favail term that changes
        # overall scaling factor for currents) and also shifts distribution of
        # mode1/2 channels. PKA actions act on both junctional and sarcolemmal
        # channels.

        # To allow for CDI KO
        # cajLCC = y36
        # caslLCC = y37
        
        # LCC Current Fixed Parameters
        taupo = 1          # [ms] - Time constant of activation
        TBa = 450          # [ms] - Time constant
        s1o = 0.0221
        k1o = 0.03
        kop = 2.5e-3       # [mM]
        cpbar = 8e-3       # [mM]
        tca = 78.0312
        ICa_scale = 5.25
        recoveryReduc = 3
        
        # PKA PHOSPHOREGULATION OF LCC AVAILABLILITY (beta subunit phosph)
        fracLCCbp0 = 0.250657 # Derived quantity - (LCCbp(baseline)/LCCbtot)
        fracLCCbpISO = 0.525870 # Derived quantity - (LCCbp(ISO)/LCCbtot)
        #a_favail=(1.50-1)/(fracLCCbpISO/fracLCCbp0-1)  # fracLCCbp ISO
        a_favail=(1.56-1)/(fracLCCbpISO/fracLCCbp0-1) # fracLCCbp ISO (x1.56 o.1 ISO)
        favail = (1-a_favail)+a_favail*(LCCb_PKAp/fracLCCbp0) # Test (max x2.52 100# phosph)
        #favail = 1  # no PKA effect on LTCCb
        ICa_scale =  ICa_scale*favail
        
        SSAshift=0 
        SSIshift=0
        # Voltage- and Ca-dependent Parameters
        @variables poss(t) poss ~ 1/(1+exp(-(y39+SSAshift)/8))
        @variables fcaj(t) fcaj ~ 1/(1+(kop/y36)^3) # cajLCC = y36
        @variables Rv(t) Rv ~ 10 + 4954*exp(y39/15.6)
        @variables PrLCC(t) PrLCC ~ 1-1/(1+exp(-(y39+40)/4))
        @variables PsLCC(t) PsLCC ~ 1/(1+exp(-(y39+40+SSIshift)/11.32))
        @variables TCaj(t) TCaj ~ (tca + 0.1*(1+(y36/cpbar)^2))/(1+(y36/cpbar)^2)
        @variables tauCaj(t) tauCaj ~ (Rv-TCaj)*PrLCC + TCaj    
        @variables tauBa(t) tauBa ~ (Rv-TBa)*PrLCC + TBa
        
        # Tranisition Rates (20 rates)
        @variables alphaLCC(t) alphaLCC ~ poss/taupo
        @variables betaLCC(t) betaLCC ~ (1-poss)/taupo
        r1 = 0.3                            # [1/ms] - Opening rate
        r2 = 3                              # [1/ms] - closing rate
        @variables s1(t) s1 ~ s1o*fcaj 
        s1p = 0.00195                       # [ms] - Inactivation rate
        @variables k1(t) k1 ~ k1o*fcaj  
        k1p = 0.00413                       # [ms] - Inactivation rate
        k2 = 1e-4                           # [ms] - Inactivation rate
        k2p = 0.00224                       # [ms] - Inactivation rate
        @variables s2(t) s2 ~ s1*(k2/k1)*(r1/r2)
        s2p = s1p*(k2p/k1p)*(r1/r2)
        @variables k3(t) k3 ~ exp(-(y39+40)/3)/(3*(1+exp(-(y39+40)/3)))
        @variables k3p(t) k3p ~ k3
        @variables k6(t) k6 ~ (fcaj*PsLCC)/tauCaj
        # k5 = (1-PsLCC)/tauCaj
        # k5p = (1-PsLCC)/tauBa
        
        # Recovery terms
        # k5 = k5/recoveryReduc
        @variables k5(t) k5 ~ (1-PsLCC)/(tauCaj*recoveryReduc)
        # k5p = k5p/recoveryReduc
        @variables k5p(t) k5p ~ (1-PsLCC)/(tauBa*recoveryReduc)
        
        @variables k6p(t) k6p ~ PsLCC/tauBa
        @variables k4(t) k4 ~ k3*(alphaLCC/betaLCC)*(k1/k2)*(k5/k6)
        @variables k4p(t) k4p ~ k3p*(alphaLCC/betaLCC)*(k1p/k2p)*(k5p/k6p)
        
        # global gates
        # gates(1,tStep) = s1
        # gates(2,tStep) = k1
        
        # State transitions for MODE 1 junctional LCCs
        # O = no differential C2 = 60 C1 = 61 I1Ca = 62 I2Ca = 63
        # I1Ba = 64 I2Ba = 65
        @variables Po_LCCj_m1(t) Po_LCCj_m1 ~ 1.0-y60-y61-y62-y63-y64-y65                    # O_m1j
        du[60] = betaLCC*y61 + k5*y63 + k5p*y65 - (k6+k6p+alphaLCC)*y60                     # C2_m1j
        du[61] = alphaLCC*y60 + k2*y62 + k2p*y64 + r2*Po_LCCj_m1 - (r1+betaLCC+k1+k1p)*y61  # C1_m1j
        du[62] = k1*y61 + k4*y63 + s1*Po_LCCj_m1 - (k2+k3+s2)*y62                           # I1Ca_m1j
        du[63] = k3*y62 + k6*y60 - (k4+k5)*y63                                              # I2Ca_m1j
        du[64] = k1p*y61 + k4p*y65 + s1p*Po_LCCj_m1 - (k2p+k3p+s2p)*y64                     # I1Ba_m1j
        du[65] = k3p*y64 + k6p*y60 - (k5p+k4p)*y65                                          # I2Ba_m1j
        @variables ibarca_jm1(t) ibarca_jm1 ~ (4*pCa*y39*Frdy*FoRT)*(0.001*exp(2*y39*FoRT)-0.341*Cao)/(exp(2*y39*FoRT)-1)
        @variables I_Ca_junc_m1(t) I_Ca_junc_m1 ~ (Fjunc_CaL*ibarca_jm1*Po_LCCj_m1*Q10CaL^Qpow)*ICa_scale
        
        # Re-define all parameters as mode 2 specific parameters
        s1om2 = 0.0221
        k1om2 = 0.03
        kopm2 = 2.5e-3
        cpbarm2 = 8e-3
        tcam2 = 78.0312
        
        @variables possm2(t) possm2 ~ 1/(1+exp(-(y39+SSAshift)/8))
        @variables fcajm2(t) fcajm2 ~ 1/(1+(kopm2/y36)^3) # Depends on junctional Ca
        @variables Rvm2(t) Rvm2 ~ 10 + 4954*exp(y39/15.6)
        @variables PrLCCm2(t) PrLCCm2 ~ 1-1/(1+exp(-(y39+40)/4))
        @variables PsLCCm2(t) PsLCCm2 ~ 1/(1+exp(-(y39+40+SSIshift)/11.32))
        @variables TCajm2(t) TCajm2 ~ (tcam2 + 0.1*(1+(y36/cpbarm2)^2))/(1+(y36/cpbarm2)^2) # Caj dependent
        @variables tauCajm2(t) tauCajm2 ~ (Rvm2-TCajm2)*PrLCCm2 + TCajm2 # Caj dependence
        @variables tauBam2(t) tauBam2 ~ (Rvm2-TBa)*PrLCCm2 + TBa
        
        @variables alphaLCCm2(t) alphaLCCm2 ~ possm2/taupo
        @variables betaLCCm2(t) betaLCCm2 ~ (1-possm2)/taupo
        r1m2 = 0.3                               # [1/ms] - Opening rate
        r2m2 = 3/8 # [1/ms] - closing rate,  changed from rabbit (3/10) - MOUSE
        @variables s1m2(t) s1m2 ~ s1om2*fcajm2 
        s1pm2 = 0.00195                           # [ms] - Inactivation rate
        @variables k1m2(t) k1m2 ~ k1om2*fcajm2
        k1pm2 = 0.00413                           # [ms] - Inactivation rate
        k2m2 = 1e-4                              # [ms] - Inactivation rate
        k2pm2 = 0.00224                           # [ms] - Inactivation rate
        @variables s2m2(t) s2m2 ~ s1m2*(k2m2/k1m2)*(r1m2/r2m2)
        s2pm2 = s1pm2*(k2pm2/k1pm2)*(r1m2/r2m2)
        @variables k3m2(t) k3m2 ~ exp(-(y(39)+40)/3)/(3*(1+exp(-(y(39)+40)/3)))
        @variables k3pm2(t) k3pm2 ~ k3m2
        # k5m2 = (1-PsLCCm2)/tauCajm2
        @variables k6m2(t) k6m2 ~ (fcajm2*PsLCCm2)/tauCajm2
        # k5pm2 = (1-PsLCCm2)/tauBam2
        # k5m2 = k5m2/recoveryReduc      # reduced for recovery
        # k5pm2 = k5pm2/recoveryReduc    # reduced for recovery
        @variables k5m2(t) k5m2 ~ (1-PsLCCm2)/(tauCajm2*recoveryReduc)
        @variables k5pm2(t) k5pm2 ~ (1-PsLCCm2)/(tauBam2*recoveryReduc)
        @variables k6pm2(t) k6pm2 ~ PsLCCm2/tauBam2
        @variables k4m2(t) k4m2 ~ k3m2*(alphaLCCm2/betaLCCm2)*(k1m2/k2m2)*(k5m2/k6m2)
        @variables k4pm2(t) k4pm2 ~ k3pm2*(alphaLCCm2/betaLCCm2)*(k1pm2/k2pm2)*(k5pm2/k6pm2)
        
        # State transitions for MODE 2 junctional LCCs
        # O = no differential C2 = 66 C1 = 67 I1Ca = 68 I2Ca = 69
        # I1Ba = 70 I2Ba = 71
        @variables Po_LCCj_m2(t) Po_LCCj_m2 ~ 1.0-y66-y67-y68-y69-y70-y71                                                # O_m2j
        du[66] = betaLCCm2*y67 + k5m2*y69 + k5pm2*y71 - (k6m2+k6pm2+alphaLCCm2)*y66                                     #C2_m2j
        du[67] = alphaLCCm2*y66 + k2m2*y68 + k2pm2*y70 + r2m2*Po_LCCj_m2 - (r1m2+betaLCCm2+k1m2+k1pm2)*y67              # C1_m2j
        du[68] = k1m2*y67 + k4m2*y69 + s1m2*Po_LCCj_m2 - (k2m2+k3m2+s2m2)*y68                                           # I1Ca_m2j
        du[69] = k3m2*y68 + k6m2*y66 - (k4m2+k5m2)*y69                                                                  # I2Ca_m2j
        du[70] = k1pm2*y67 + k4pm2*y71 + s1pm2*Po_LCCj_m2 - (k2pm2+k3pm2+s2pm2)*y70                                     # I1Ba_m2j
        du[71] = k3pm2*y70 + k6pm2*y66 - (k5pm2+k4pm2)*y71                                                              # I2Ba_m2j
        @variables ibarca_jm2(t) ibarca_jm2 ~ (4*pCa*y39*Frdy*FoRT)*(0.001*exp(2*y39*FoRT)-0.341*Cao)/(exp(2*y39*FoRT)-1)
        @variables I_Ca_junc_m2(t) I_Ca_junc_m2 ~ (Fjunc_CaL*ibarca_jm2*(Po_LCCj_m2)*Q10CaL^Qpow)*ICa_scale
        
        # CaMKII AND PKA-DEPENDENT SHIFTING OF DYADIC LCCS TO MODE 2
        #fpkam2 = 0.1543*LCCa_PKAp - .0043 # Assumes max phosphorylation results in 15# mode 2
        fracLCCap0 = 0.219577 # Derived
        frac_fpkam2 = (0.15*fracLCCap0)/(1-fracLCCap0)
        fpkam2 = (0.15+frac_fpkam2)*LCCa_PKAp - frac_fpkam2 # Assumes max (100#) phosphorylation results in 15# mode 2 channels
        #(fpkam2 = 0 with NO ISO)
        #fpkam2 = 0 # no PKA effect on LTCCa
            #fpkam2 = 0.0765*0.8
        fckiim2 = LCC_CKp*0.1 # Assumes max phosphorylation results in 10# mode 2 channels (max LCC_CKp = 1)
        # Sum up total fraction of CKII and PKA-shifted mode 2 channels
        junc_mode2 = fckiim2 + fpkam2
        # Total junctional ICa
        @variables I_Ca_junc2(t) I_Ca_junc2 ~ (1-junc_mode2)*I_Ca_junc_m1 + junc_mode2*I_Ca_junc_m2
        
        # SUB-SARCOLEMMAL LCCs
        # Re-assign necessary params to be Casl sensitive
        @variables fcasl(t) fcasl ~ 1/(1+(kop/y37)^3)    # Depends on sl Ca
        @variables TCasl(t) TCasl ~ (tca + 0.1*(1+(y37/cpbar))^2)/(1+(y37/cpbar)^2)
        @variables tauCasl(t) tauCasl ~ (Rv-TCasl)*PrLCC + TCasl
        
        # Re-assign necessary rates to be Casl sensitive
        @variables s1sl(t) s1sl ~ s1o*fcasl
        @variables k1sl(t) k1sl ~ k1o*fcasl
        @variables s2sl(t) s2sl ~ s1sl*(k2/k1sl)*(r1/r2)
        s2psl = s1p*(k2p/k1p)*(r1/r2)
        # k5sl = (1-PsLCC)/tauCasl
        # k5sl = k5sl/recoveryReduc  # Reduced for recovery
        @variables k5sl(t) k5sl ~ (1-PsLCC)/(tauCasl*recoveryReduc)
        @variables k6sl(t) k6sl ~ (fcasl*PsLCC)/tauCasl
        @variables k4sl(t) k4sl ~ k3*(alphaLCC/betaLCC)*(k1sl/k2)*(k5sl/k6sl)
        @variables k4psl(t) k4psl ~ k3p*(alphaLCC/betaLCC)*(k1p/k2p)*(k5p/k6p)
        
        # State transitions for 'mode 1' sarcolemmal LCCs
        # O = no differential C2 = 72 C1 = 73 I1Ca = 74 I2Ca = 75
        # I1Ba = 76 I2Ba = 77
        @variables Po_LCCsl_m1(t) Po_LCCsl_m1 ~ 1-y72-y73-y74-y75-y76-y77                        # O_m1sl
        du[72] = betaLCC*y73 + k5sl*y75 + k5p*y77 - (k6sl+k6p+alphaLCC)*y72                     # C2_m1sl
        du[73] = alphaLCC*y72 + k2*y74 + k2p*y76 + r2*Po_LCCsl_m1 - (r1+betaLCC+k1sl+k1p)*y73   # C1_m1sl
        du[74] = k1sl*y73 + k4sl*y75 + s1sl*Po_LCCsl_m1 - (k2+k3+s2sl)*y74                      # I1Ca_m1sl
        du[75] = k3*y74 + k6sl*y72 - (k4sl+k5sl)*y75                                            # I2Ca_m1sl
        du[76] = k1p*y73 + k4psl*y77 + s1p*Po_LCCsl_m1 - (k2p+k3p+s2psl)*y76                    # I1Ba_m1sl
        du[77] = k3p*y76 + k6p*y72 - (k5p+k4psl)*y77                                            # I2Ba_m1sl
        @variables ibarca_slm1(t) ibarca_slm1 ~ (4*pCa*y39*Frdy*FoRT)*(0.001*exp(2*y39*FoRT)-0.341*Cao)/(exp(2*y39*FoRT)-1)
        @variables I_Casl_m1(t) I_Casl_m1 ~ (Fsl_CaL*ibarca_slm1*Po_LCCsl_m1*Q10CaL^Qpow)*ICa_scale
        
        # Adjust closing rate for 'mode 2' sarcolemmal LCCs
        r2slm2 = r2m2
        @variables s2slm2(t) s2slm2 ~ s1sl*(k2/k1sl)*(r1/r2slm2)
        s2pslm2 = s1p*(k2p/k1p)*(r1/r2slm2)
        
        # State transitions for mode 2 sarcolemmal LCCs
        # O = no differential C2 = 78 C1 = 79 I1Ca = 80 I2Ca = 81 I1Ba = 82 I2Ba = 83
        @variables Po_LCCsl_m2(t) Po_LCCsl_m2 ~ 1-y78-y79-y80-y81-y82-y83                            # O_m2sl
        du[78] = betaLCC*y79 + k5sl*y81 + k5p*y83 - (k6sl+k6p+alphaLCC)*y78                         # C2_m2sl
        du[79] = alphaLCC*y78 + k2*y80 + k2p*y82 + r2slm2*Po_LCCsl_m2 - (r1+betaLCC+k1sl+k1p)*y79   # C1_m2sl
        du[80] = k1sl*y79 + k4sl*y81 + s1sl*Po_LCCsl_m2 - (k2+k3+s2slm2)*y80                        # I1Ca_m2sl
        du[81] = k3*y80 + k6sl*y78 - (k4sl+k5sl)*y81                                                # I2Ca_m2sl
        du[82] =  k1p*y79 + k4psl*y83 + s1p*Po_LCCsl_m2 - (k2p+k3p+s2pslm2)*y82                     # I1Ba_m2sl
        du[83] = k3p*y82 + k6p*y78 - (k5p+k4psl)*y83                                                # I2Ba_m2sl
        @variables ibarca_slm2(t) ibarca_slm2 ~ (4*pCa*y39*Frdy*FoRT)*(0.001*exp(2*y39*FoRT)-0.341*Cao)/(exp(2*y39*FoRT)-1)
        @variables I_Casl_m2(t) I_Casl_m2 ~ (Fsl_CaL*ibarca_slm2*Po_LCCsl_m2*Q10CaL^Qpow)*ICa_scale
        
        # Sum mode 1 and mode 2 sl channels for total sl current
        fckiim2_sl = 0 # Set to zero since SL LCCp by CaMKII is negligible
        sl_mode2 = fckiim2_sl + fpkam2
        @variables I_Ca_sl2(t) I_Ca_sl2 ~ (1-sl_mode2)*I_Casl_m1 + sl_mode2*I_Casl_m2
        
        # Na and K currents through LCC
        @variables I_CaKj2(t) I_CaKj2 ~ ibark*Fjunc_CaL*((1-junc_mode2)*Po_LCCj_m1 + junc_mode2*Po_LCCj_m2)*Q10CaL^Qpow*ICa_scale
        @variables I_CaKsl2(t) I_CaKsl2 ~ ibark*Fsl_CaL*((1-sl_mode2)*Po_LCCsl_m1 + sl_mode2*Po_LCCsl_m2)*Q10CaL^Qpow*ICa_scale
        @variables I_CaK2(t) I_CaK2 ~ I_CaKj2+I_CaKsl2
        @variables I_CaNa_junc2(t) I_CaNa_junc2 ~ (Fjunc_CaL*ibarna_j*((1-junc_mode2)*Po_LCCj_m1+junc_mode2*Po_LCCj_m2)*Q10CaL^Qpow)*ICa_scale
        @variables I_CaNa_sl2(t) I_CaNa_sl2 ~ Fsl_CaL*ibarna_sl*((1-sl_mode2)*Po_LCCsl_m1 + sl_mode2*Po_LCCsl_m2)*Q10CaL^Qpow*ICa_scale
        
        # These are now able to switch depending on whether or not the flag to
        # switch to Markov model of ICa is ON
        @variables I_Ca_junc(t) I_Ca_junc ~ (1-ICa_MarkovFlag)*I_Ca_junc1 + ICa_MarkovFlag*I_Ca_junc2
        @variables I_Ca_sl(t) I_Ca_sl ~ (1-ICa_MarkovFlag)*I_Ca_sl1 + ICa_MarkovFlag*I_Ca_sl2
        @variables I_Ca(t) I_Ca ~ I_Ca_junc+I_Ca_sl   # Total Ca curren throuhgh LCC
        @variables I_CaNa_junc(t) I_CaNa_junc ~ (1-ICa_MarkovFlag)*(I_CaNa_junc1) + (ICa_MarkovFlag)*(I_CaNa_junc2)
        @variables I_CaNa_sl(t) I_CaNa_sl ~ (1-ICa_MarkovFlag)*(I_CaNa_sl1) + (ICa_MarkovFlag)*(I_CaNa_sl2)
        @variables I_CaNa(t) I_CaNa ~ I_CaNa_junc + I_CaNa_sl   # Total Na current through LCC
        @variables I_CaK(t) I_CaK ~ (1-ICa_MarkovFlag)*(I_CaK1) + ICa_MarkovFlag*(I_CaK2)  # Total K current through LCC
        
        # Collect all currents through LCC
        @variables I_Catot(t) I_Catot ~ I_Ca+I_CaK+I_CaNa
        du[43] = -I_Ca*Cmem/(Vmyo*2*Frdy)*1e3
        
        # global I_Ca_store ibar_store
        # I_Ca_store(tStep) = I_Catot
        # ibar_store(tStep) = ibarca_j

    ## I_ncx: Na/Ca Exchanger flux

        @variables Ka_junc(t) Ka_junc ~ 1/(1+(Kdact/y36)^3)
        @variables Ka_sl(t) Ka_sl ~ 1/(1+(Kdact/y37)^3)
        @variables s1_junc(t) s1_junc ~ exp(nu*y39*FoRT)*y32^3*Cao
        @variables s1_sl(t) s1_sl ~ exp(nu*y39*FoRT)*y33^3*Cao
        @variables s2_junc(t) s2_junc ~ exp((nu-1)*y39*FoRT)*Nao^3*y36
        @variables s3_junc(t) s3_junc ~ (KmCai*Nao^3*(1+(y32/KmNai)^3)+KmNao^3*y36+ KmNai^3*Cao*(1+y36/KmCai)+KmCao*y32^3+y32^3*Cao+Nao^3*y36)*(1+ksat*exp((nu-1)*y39*FoRT))
        @variables s2_sl(t) s2_sl ~ exp((nu-1)*y39*FoRT)*Nao^3*y37
        @variables s3_sl(t) s3_sl ~ (KmCai*Nao^3*(1+(y33/KmNai)^3) + KmNao^3*y37+KmNai^3*Cao*(1+y37/KmCai)+KmCao*y33^3+y33^3*Cao+Nao^3*y37)*(1+ksat*exp((nu-1)*y39*FoRT))
        @variables I_ncx_junc(t) I_ncx_junc ~ Fjunc_ncx*IbarNCX*Q10NCX^Qpow*Ka_junc*(s1_junc-s2_junc)/s3_junc
        @variables I_ncx_sl(t) I_ncx_sl ~ Fsl_ncx*IbarNCX*Q10NCX^Qpow*Ka_sl*(s1_sl-s2_sl)/s3_sl
        @variables I_ncx(t) I_ncx ~ I_ncx_junc+I_ncx_sl
        du[45] = 2*I_ncx*Cmem/(Vmyo*2*Frdy)*1e3     # uM/ms
        
        #　global Incx  
        # Incx(tStep) = I_ncx

    ## I_pca: Sarcolemmal Ca Pump Current

        @variables I_pca_junc(t) I_pca_junc ~ Fjunc*Q10SLCaP^Qpow*IbarSLCaP*y(36)^1.6/(KmPCa^1.6+y36^1.6)
        @variables I_pca_sl(t) I_pca_sl ~ Fsl*Q10SLCaP^Qpow*IbarSLCaP*y(37)^1.6/(KmPCa^1.6+y37^1.6)
        @variables I_pca(t) I_pca ~ I_pca_junc+I_pca_sl
        du[44] = -I_pca*Cmem/(Vmyo*2*Frdy)*1e3
        
        # global Ipca_store
        # Ipca_store(tStep) = I_pca
    ## I_cabk: Ca Background Current
        
        @variables I_cabk_junc(t) I_cabk_junc ~ Fjunc*GCaB*(y39-eca_junc)
        @variables I_cabk_sl(t) I_cabk_sl ~ Fsl*GCaB*(y39-eca_sl)
        @variables I_cabk(t) I_cabk ~ I_cabk_junc+I_cabk_sl
        du[46] =-I_cabk*Cmem/(Vmyo*2*Frdy)*1e3

    ## I_CFTR or I_cl_(cAMP) - Cystic Fibrosis Transmembrane Conductance Reg.
        # This is an Em- and time-independent current that is activated by PKA
        # fact_pka_cftr = 1.1933*ICFTR_PKAp - 0.1933    # Derived?
        # gCFTR = fact_pka_cftr*4.9e-3                  # [A/F] - Max value as in Shannon et al. (2005)
        Icftr = 0       # gCFTR*(y(39) - ecl)           # NO Icftr in MOUSE
        
        # global ICFTR
        # ICFTR(tStep) = Icftr

    ## RyR model - SR release fluxes and leak

        # CaMKII and PKA-dependent phosphoregulation of RyR Po
        fCKII_ec50SR = 1.16 - 4/5*RyR_CKp
        ec50SR = fCKII_ec50SR*ec50SR # MOUSE - 60# 
        
        MaxSR = 15 
        MinSR = 1
        @variables kCaSR(t) kCaSR ~ MaxSR - (MaxSR-MinSR)/(1+(ec50SR/y31)^2.5)
        # koSRCa = koCa/kCaSR
        @variables kiSRCa(t) kiSRCa ~ kiCa*kCaSR
        kleak = 2*5.348e-6 # [1/ms] changed from rabbit (5.348e-6)
        
        #fCKII_RyR = (20*RyR_CKp/3 - 1/3)  # 1 at basal condition - RABBIT
        fCKII_RyR = (10*RyR_CKp - 1) # 1 at basal condition - MOUSE
        
        #fPKA_RyR = RyR_PKAp*1.025 + 0.9750  # 1 with NO ISO
        frac_RyRo = 0.204276 # Derived (RyR_PKAp(basal)/RyRtot)
        a_RyR = (2-1)/(1/frac_RyRo-1) # Max effect: fPKA_RyR=2
        fPKA_RyR = 1-a_RyR+a_RyR*(RyR_PKAp/frac_RyRo)
        @variables koSRCa(t) koSRCa ~ (fCKII_RyR + fPKA_RyR - 1)*(koCa/kCaSR) # *koSRCa = *(koCa/kCaSR)
        
        # ODEs for RyR states and SR release through open RyRs
        @variables RI(t) RI ~ 1-y14-y15-y16
        du[14] = (kim*RI-kiSRCa*y36*y14)-(koSRCa*y36^2*y14-kom*y15)   # R
        du[15] = (koSRCa*y36^2*y14-kom*y15)-(kiSRCa*y36*y15-kim*y16)  # O
        du[16] = (kiSRCa*y36*y15-kim*y16)-(kom*y16-koSRCa*y36^2*RI)   # I
        @variables J_SRCarel(t) J_SRCarel ~ ks*y15*(y31-y36)           # [mmol/L SR/ ms]
        
        # Passive RyR leak - includes CaMKII regulation of leak flux
        #kleak = (1/3 + 10*RyR_CKp/3)*kleak # RABBIT
        kleak = (1/2 + 5*RyR_CKp/2)*kleak # MOUSE (reduced CaMKII effect on leak)
        @variables J_SRleak(t) J_SRleak ~ kleak*(y31-y36) # [mmol/L cyt/ms]
        
        # global Jleak 
        # Jleak(tStep,1) = J_SRCarel*Vsr/Vmyo + J_SRleak    # Total leak [mmol/L cyt/ms]
        # Jleak(tStep,2) = J_SRleak                         # Passive leak only [mmol/L cyt/ms] 

    ## SERCA model - SR uptake fluxes

        # CaMKII and PKA-dependent phosphoregulation of PLB (changes to SERCA flux)
        fCKII_PLB = (1-0.5*PLB_CKp)                         # Max effect: fCKII_PLB=0.5
        fracPKA_PLBo = 1-0.079755                           # Derived quantity - (1 - (PLBp(baseline)/PLBtot))
        #fPKA_PLB = (PLB_PKAn/fracPKA_PLBo)*3/4 + 1/4       # Max effect: fPKA_PLB=0.25
        fPKA_PLB = (PLB_PKAn/fracPKA_PLBo)*(100-55.31)/100 + 55.31/100      # Max effect: fPKA_PLB=0.45
        
        # Select smaller value (resulting in max reduction of Kmf)
        if fCKII_PLB < fPKA_PLB
            Kmf = Kmf*fCKII_PLB     #fCKII_PLB
        elseif fPKA_PLB < fCKII_PLB
            Kmf = Kmf*fPKA_PLB      #fPKA_PLB
        end
        
        @variables J_serca(t)
        J_serca ~ Q10SRCaP^Qpow*Vmax_SRCaP*((y38/Kmf)^hillSRCaP-(y31/Kmr)^hillSRCaP)/(1+(y38/Kmf)^hillSRCaP+(y31/Kmr)^hillSRCaP) # [mM/msec] 
        
        # global Jserca
        # Jserca(tStep) = J_serca       # [mM/msec] or [mmol/L cyt msec]

    ## Na and Ca Buffering

        du[17] = kon_na*y32*(Bmax_Naj-y17)-koff_na*y17       # NaBj      [mM/ms]
        du[18] = kon_na*y33*(Bmax_Nasl-y18)-koff_na*y18      # NaBsl     [mM/ms]
        
        # Cytosolic Ca Buffers
        du[19] = kon_tncl*y38*(Bmax_TnClow-y19)-koff_tncl*y19           # TnCL      [mM/ms]
        du[20] = kon_tnchca*y38*(Bmax_TnChigh-y20-y21)-koff_tnchca*y20  # TnCHc     [mM/ms]
        du[21] = kon_tnchmg*Mgi*(Bmax_TnChigh-y20-y21)-koff_tnchmg*y21  # TnCHm     [mM/ms]
        du[22] = 0   # commented b/c buffering done by CaM module // kon_cam*y(38)*(Bmax_CaM-y(22))-koff_cam*y(22)  # CaM       [mM/ms]
        # kon_cam*y(38)*(Bmax_CaM-y(22))-koff_cam*y(22) # CaM       [mM/ms]
        du[23] = kon_myoca*y38*(Bmax_myosin-y23-y24)-koff_myoca*y23     # Myosin_ca [mM/ms]
        du[24] = kon_myomg*Mgi*(Bmax_myosin-y23-y24)-koff_myomg*y24     # Myosin_mg [mM/ms]
        du[25] = kon_sr*y38*(Bmax_SR-y25)-koff_sr*y25                   # SRB       [mM/ms]
        #J_CaB_cytosol = sum(ydot(19:25))   # wrong formulation
        J_CaB_cytosol = du[19] + du[20] + du[22] + du[23] + du[25]
        
        # Junctional and SL Ca Buffers
        du[26] = kon_sll*y36*(Bmax_SLlowj-y26)-koff_sll*y26             # SLLj      [mM/ms]
        du[27] = kon_sll*y37*(Bmax_SLlowsl-y27)-koff_sll*y27            # SLLsl     [mM/ms]
        du[28] = kon_slh*y36*(Bmax_SLhighj-y28)-koff_slh*y28            # SLHj      [mM/ms]
        du[29] = kon_slh*y37*(Bmax_SLhighsl-y29)-koff_slh*y29           # SLHsl     [mM/ms]
        J_CaB_junction = du[26] + du[28]
        J_CaB_sl = du[27] + du[29]

    ## Ion concentrations

        # SR Ca Concentrations
        du[30] = kon_csqn*y31*(Bmax_Csqn-y30)-koff_csqn*y30             # Csqn      [mM/ms]
        du[31] = J_serca*Vmyo/Vsr-(J_SRleak*Vmyo/Vsr+J_SRCarel)-du[30]  # Ca_sr     [mM/ms] # Ratio 3 leak current
        
        # Na Concentrations
        @variables I_Na_tot_junc(t) I_Na_tot_junc ~ I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc   # [uA/uF]
        @variables I_Na_tot_sl(t) I_Na_tot_sl ~ I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl   #[uA/uF]
        du[32] = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(y33-y32)-du[17]
        du[33] = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(y32-y33)+J_na_slmyo/Vsl*(y34-y33)-du[18]
        if NaClampFlag == 1
            du[34] = 0                              # Na clamp
        else
            du[34] = J_na_slmyo/Vmyo*(y33-y34)      # [mM/msec]
        end
        
        # K Concentration
        @variables I_K_tot(t) I_K_tot ~ I_to+I_kr+I_ks+I_ki-2*I_nak+I_CaK+I_kp+I_kur+I_ss    # [uA/uF]
        du[35] = 0          #-I_K_tot*Cmem/(Vmyo*Frdy)                                      # [mM/msec]
        
        # Ca Concentrations
        @variables I_Ca_tot_junc(t) I_Ca_tot_junc ~ I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc    # [uA/uF]
        @variables I_Ca_tot_sl(t) I_Ca_tot_sl ~ I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl                # [uA/uF]
        # Ca_j
        du[36] = -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(y37-y36)-J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc   
        # Ca_sl 
        du[37] = -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(y36-y37)+ J_ca_slmyo/Vsl*(y38-y37)-J_CaB_cytosol         
        # Cai                  
        du[38] = -J_serca-J_CaB_cytosol +J_ca_slmyo/Vmyo*(y37-y38)
        @variables junc_sl(t) junc_sl ~ J_ca_juncsl/Vsl*(y36-y37)
        @variables sl_junc(t) sl_junc ~ J_ca_juncsl/Vjunc*(y37-y36)
        @variables sl_myo(t) sl_myo ~ J_ca_slmyo/Vsl*(y38-y37)
        @variables myo_sl(t) myo_sl ~ J_ca_slmyo/Vmyo*(y37-y38)

    ## Simulation type

        protocol = 1        # "pace"
        
        if CaffeineFlag==1
            protocol = 2    # "vcRest"
        end
        if StrophFlag == 1
            protocol = 3    # "none"
        end
        
        # AP Waveform for AP clamp
        
        if protocol == 1
            if mod(t,cycleLength) <= 5 
                I_app = 9.5
            else
                I_app = 0.0
            end
        elseif protocol == 2 
            V_clamp = -83   #y(39)#-80
            R_clamp = 0.01
            I_app = (V_clamp-y39)/R_clamp
        else
            I_app = 0
        end          

    ## Membrane Potential

        @variables I_Na_tot(t) I_Na_tot ~ I_Na_tot_junc + I_Na_tot_sl                 # [uA/uF]
        @variables I_Cl_tot(t) I_Cl_tot ~ I_ClCa+I_Clbk+Icftr                        # [uA/uF]
        @variables I_Ca_tot(t) I_Ca_tot ~ I_Ca_tot_junc+I_Ca_tot_sl                   # [uA/uF]
        @variables I_tot(t) I_tot ~ I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot            # [uA/uF]
        du[39] = -(I_tot-I_app)

        # global dVm_store 
        # dVm_store(tStep) = ydot(39)


        return du
    end

    export morotti_eccODE
end