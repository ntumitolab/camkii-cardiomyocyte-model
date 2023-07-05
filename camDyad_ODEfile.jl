module Module_camDyad
    using DifferentialEquations
    using ModelingToolkit
    # This function describes the ODE's for CaM, CaMKII, and CaN.
    
    function morotti_camDyadODE!(du, u,  p, CaDyad,t)
    ## State variables
        ## Descriptions for state variables

        # Ca2CaM = y(1)  # 2 Ca bound to C terminal sites
        # Ca4CaM = y(2)  # 4 Ca bound
        # CaMB = y(3)   
        # Ca2CaMB = y(4)
        # Ca4CaMB = y(5)
        # Pb2 = y(6)     # probability of a Ca2CaM bound CaMKII subunit
        # Pb = y(7)      # probability of a Ca4CaM bound CaMKII subunit
        # Pt = y(8)      # probability of a Ca4CaM bound autophosphorylated CaMKII subunit
        # Pt2 = y(9)     # probability of a Ca2CaM bound autophosphorylated CaMKII subunit
        # Pa = y(10)     # probability of an autonomous autophosphorylated CaMKII subunit
        # CaMCaN = y(11)
        # Ca2CaMCaN = y(12)
        # Ca4CaMCaN = y(13)  # active calcineurin
        # description of intermediate variables
        # CaM- Ca free CaM
        (CaM,Ca2CaM,Ca4CaM,CaMB,Ca2CaMB,Ca4CaMB,Pb2,Pb,Pt,Pt2,Pa,
         Ca4CaN,CaMCa4CaN,Ca2CaMCa4CaN,Ca4CaMCa4CaN) = u

    ## Parameters    
        (K,Mg,CaMtot,Btot,CaMKIItot,CaNtot,PP1tot,Ca,cyclelength,compartment) = p

        
    ## Parameters
        # Ca/CaM parameters
        if Mg <= 1
            Kd02 = 0.0025*(1+K/0.94-Mg/0.012)*(1+K/8.1+Mg/0.022)  # [uM^2]
            Kd24 = 0.128*(1+K/0.64+Mg/0.0014)*(1+K/13.0-Mg/0.153) # [uM^2]
        else
            Kd02 = 0.0025*(1+K/0.94-1/0.012+(Mg-1)/0.060)*(1+K/8.1+1/0.022+(Mg-1)/0.068)   # [uM^2]
            Kd24 = 0.128*(1+K/0.64+1/0.0014+(Mg-1)/0.005)*(1+K/13.0-1/0.153+(Mg-1)/0.150)  # [uM^2]
        end
        k20 = 10               # [s^-1]      
        k02 = k20/Kd02         # [uM^-2 s^-1]
        k42 = 500              # [s^-1]      
        k24 = k42/Kd24         # [uM^-2 s^-1]
        
        # CaM buffering (B) parameters
        k0Boff = 0.0014        # [s^-1] 
        k0Bon = k0Boff/0.2   # [uM^-1 s^-1] kon = koff/Kd
        k2Boff = k0Boff/100    # [s^-1] 
        k2Bon = k0Bon          # [uM^-1 s^-1]
        k4Boff = k2Boff        # [s^-1]
        k4Bon = k0Bon          # [uM^-1 s^-1]

        # using thermodynamic constraints
        k20B = k20/100 # [s^-1] thermo constraint on loop 1
        k02B = k02     # [uM^-2 s^-1] 
        k42B = k42     # [s^-1] thermo constraint on loop 2
        k24B = k24     # [uM^-2 s^-1]
        
        # CaMKII parameters
        # Wi Wa Wt Wp
        kbi = 2.2      # [s^-1] (Ca4CaM dissocation from Wb)
        kib = kbi/33.5e-3 # [uM^-1 s^-1]
        kib2 = kib
        kb2i = kib2*5
        kb24 = k24
        kb42 = k42*33.5e-3/5
        kpp1 = 1.72     # [s^-1] (PP1-dep dephosphorylation rates)
        Kmpp1 = 11.5    # [uM]
        kta = kbi/1000  # [s^-1] (Ca4CaM dissociation from Wt)
        kat = kib       # [uM^-1 s^-1] (Ca4CaM reassociation with Wa)
        kt42 = k42*33.5e-6/5
        kt24 = k24
        kat2 = kib
        kt2a = kib*5
        
        # CaN parameters
        kcanCaoff = 1              # [s^-1] 
        kcanCaon = kcanCaoff/0.5   # [uM^-1 s^-1] 
        kcanCaM4on = 46            # [uM^-1 s^-1]
        kcanCaM4off = 1.3e-3       # [s^-1]
        kcanCaM2on = kcanCaM4on
        kcanCaM2off = 2508*kcanCaM4off
        kcanCaM0on = kcanCaM4on
        kcanCaM0off = 165*kcanCaM2off
        k02can = k02
        k20can = k20/165
        k24can = k24
        k42can = k20/2508
        
    ## Fluxes
        # CaM Reaction fluxes
        @variables rcn02(t) rcn02 ~ k02*Ca^2*CaM - k20*Ca2CaM
        @variables rcn24(t) rcn24 ~ k24*Ca^2*Ca2CaM - k42*Ca4CaM

        # CaM buffer fluxes
        @variables B(t) B ~ Btot - CaMB - Ca2CaMB - Ca4CaMB
        @variables rcn02B(t) rcn02B ~ k02B*Ca^2*CaMB - k20B*Ca2CaMB
        @variables rcn24B(t) rcn24B ~ k24B*Ca^2*Ca2CaMB - k42B*Ca4CaMB
        @variables rcn0B(t) rcn0B ~ k0Bon*CaM*B - k0Boff*CaMB
        @variables rcn2B(t) rcn2B ~ k2Bon*Ca2CaM*B - k2Boff*Ca2CaMB
        @variables rcn4B(t) rcn4B ~ k4Bon*Ca4CaM*B - k4Boff*Ca4CaMB

        # CaN reaction fluxes 
        @variables Ca2CaN(t) Ca2CaN ~ CaNtot - Ca4CaN - CaMCa4CaN - Ca2CaMCa4CaN - Ca4CaMCa4CaN
        @variables rcnCa4CaN(t) rcnCa4CaN ~ kcanCaon*Ca^2*Ca2CaN - kcanCaoff*Ca4CaN
        @variables rcn02CaN(t) rcn02CaN ~ k02can*Ca^2*CaMCa4CaN - k20can*Ca2CaMCa4CaN
        @variables rcn24CaN(t) rcn24CaN ~ k24can*Ca^2*Ca2CaMCa4CaN - k42can*Ca4CaMCa4CaN
        @variables rcn0CaN(t) rcn0CaN ~ kcanCaM0on*CaM*Ca4CaN - kcanCaM0off*CaMCa4CaN
        @variables rcn2CaN(t) rcn2CaN ~ kcanCaM2on*Ca2CaM*Ca4CaN - kcanCaM2off*Ca2CaMCa4CaN
        @variables rcn4CaN(t) rcn4Ca ~ kcanCaM4on*Ca4CaM*Ca4CaN - kcanCaM4off*Ca4CaMCa4CaN

        # CaMKII reaction fluxes
        @variables Pi(t) Pi = 1 - Pb2 - Pb - Pt - Pt2 - Pa
        @variables rcnCKib2(t) rcnCKib2 = kib2*Ca2CaM*Pi - kb2i*Pb2
        @variables rcnCKb2b(t) rcnCKb2b = kb24*Ca^2*Pb2 - kb42*Pb
        @variables rcnCKib(t) rcnCKib = kib*Ca4CaM*Pi - kbi*Pb
        @variables T(t) T = Pb + Pt + Pt2 + Pa
        @variables kbt(t) kbt = 0.055*T + 0.0074*T^2 + 0.015*T^3
        @variables rcnCKbt(t) rcnCKbt = kbt*Pb - kpp1*PP1tot*Pt/(Kmpp1+CaMKIItot*Pt)
        @variables rcnCKtt2(t) rcnCKtt = kt42*Pt - kt24*Ca^2*Pt2
        @variables rcnCKta(t) rcnCKta = kta*Pt - kat*Ca4CaM*Pa
        @variables rcnCKt2a(t) rcnCKt2a = kt2a*Pt2 - kat2*Ca2CaM*Pa
        @variables rcnCKt2b2(t) rcnCKt2b2 = kpp1*PP1tot*Pt2/(Kmpp1+CaMKIItot*Pt2)
        @variables rcnCKai(t) rcnCKai = kpp1*PP1tot*Pa/(Kmpp1+CaMKIItot*Pa)

    ## Equations

        # CaM equations
        du[1] = 1e-3*(-rcn02 - rcn0B - rcn0CaN)                                             # CaM
        du[2] = 1e-3*(rcn02 - rcn24 - rcn2B - rcn2CaN + CaMKIItot.*(-rcnCKib2 + rcnCKt2a) ) # Ca2CaM
        du[3] = 1e-3*(rcn24 - rcn4B - rcn4CaN + CaMKIItot.*(-rcnCKib+rcnCKta) )             # Ca4CaM
        du[4] = 1e-3*(rcn0B-rcn02B)                                                         # CaMB
        du[5] = 1e-3*(rcn02B + rcn2B - rcn24B)                                              # Ca2CaMB
        du[6] = 1e-3*(rcn24B + rcn4B)                                                       # Ca4CaMB
        
        # CaMKII equations
        du[7] = 1e-3*(rcnCKib2 - rcnCKb2b + rcnCKt2b2)  # Pb2
        du[8] = 1e-3*(rcnCKib + rcnCKb2b - rcnCKbt)     # Pb
        du[9] = 1e-3*(rcnCKbt-rcnCKta-rcnCKtt2)         # Pt
        du[10] = 1e-3*(rcnCKtt2-rcnCKt2a-rcnCKt2b2)     # Pt2
        du[11] = 1e-3*(rcnCKta+rcnCKt2a-rcnCKai)        # Pa
        
        # CaN equations
        du[12] = 1e-3*(rcnCa4CaN - rcn0CaN - rcn2CaN - rcn4CaN) # Ca4CaN
        du[13] = 1e-3*(rcn0CaN - rcn02CaN)                      # CaMCa4CaN
        du[14] = 1e-3*(rcn2CaN+rcn02CaN-rcn24CaN)               # Ca2CaMCa4CaN
        du[15] = 1e-3*(rcn4CaN+rcn24CaN)                        # Ca4CaMCa4CaN

    ## write to global variables for adjusting Ca buffering in EC coupling model
        @variables CaDyad(t) CaDyad ~ 1e-3*(2*CaMKIItot*(rcnCKtt2-rcnCKb2b) - 2*(rcn02+rcn24+rcn02B+rcn24B+rcnCa4CaN+rcn02CaN+rcn24CaN))   # [uM/msec]
    
        return du
    end

    export morotti_camDyadODE
    export CaDyad

end