include("ECC_ODEfile.jl")
include("camDyad_ODEfile.jl")
include("camSL_ODEfile.jl")
include("camCyt_ODEfile.jl")
include("CAMKII_ODEfile.jl")
include("BAR_ODEfile.jl")
using DifferentialEquations
using XLSX
using Plots
using .Module_BAR
using .Module_camDyad
using .Module_camSL
using .Module_camCyt
using .Module_CAMKII
using .Module_ECC
using ModelingToolkit

## This file merges morotti_et_al_mouse_masterODEfile.m and morotti_et_al_mouse_masterCompute.m
## This file calls the ODE function files for EC coupling, CaM reactions, CaMKII and PKA phosphorylation modules.

K = 135     # [mM]
Mg = 1      # [mM]

## Parameters for external modules

# ECC and CaM modules
freq = 1                    # [Hz] - CHANGE DEPENDING ON FREQUENCY
# cycleLength = 1e3/freq      # [ms]
# CaMtotDyad = 418            # [uM]
# BtotDyad = 1.54/8.293e-4    # [uM]
# CaMKIItotDyad = 120         # [uM]
# CaNtotDyad = 3e-3/8.293e-4  # [uM]
# PP1totDyad = 96.5           # [uM]
# CaMtotSL = 5.65             # [uM]
# BtotSL = 24.2               # [uM]
# CaMKIItotSL = 120*8.293e-4  # [uM]
# CaNtotSL = 3e-3             # [uM]
# PP1totSL = 0.57             # [uM]
# CaMtotCyt = 5.65            # [uM]
# BtotCyt = 24.2              # [uM]
# CaMKIItotCyt = 120*8.293e-4 # [uM]
# CaNtotCyt = 3e-3            # [uM] 
# PP1totCyt = 0.57            # [uM]

# ADJUST CaMKII ACTIVITY LEVELS (expression = "WT", "OE", or "KO")
expression = 1  # expression = 1 ("WT"), =2 ("OE"), or =3 ("KO")

if expression == 1          # "WT"
    tempCKIIOE = 0
elseif expression == 2     # "OE"
    tempCKIIOE = 1 # Flag for CaMKKII-OE (0=WT, 1=OE)
    n_OE=6
    CaMKIItotDyad = 120*n_OE            # [uM] 
    CaMKIItotSL = 120*8.293e-4*n_OE     # [uM]
    CaMKIItotCyt = 120*8.293e-4*n_OE    # [uM]
else                         # "KO"
    tempCKIIOE = 0
    CaMKIItotDyad = 0                   # [uM] 
    CaMKIItotSL = 0                     # [uM]
    CaMKIItotCyt = 0                    # [uM]
end

#=
CKIIOE = 0 # Should be zero during 'WT' and 'KO' runs
if strcmp(expression,"OE") # OE
    CKIIOE = 1 # Flag for CaMKKII-OE (0=WT, 1=OE)
    n_OE=6
    CaMKIItotDyad = 120*n_OE         # [uM] 
    CaMKIItotSL = 120*8.293e-4*n_OE  # [uM]
    CaMKIItotCyt = 120*8.293e-4*n_OE # [uM]
elseif strcmp(expression,"KO")
    CaMKIItotDyad = 0          # [uM] 
    CaMKIItotSL = 0            # [uM]
    CaMKIItotCyt = 0           # [uM]
end
=#

#plb_val=38 # RABBIT
plb_val=106 # MOUSE

# Parameters for CaMKII module
# LCCtotDyad = 31.4*0.9      # [uM] - Total Dyadic [LCC] - (umol/l dyad)
# LCCtotSL = 0.0846          # [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
# RyRtot = 382.6             # [uM] - Total RyR (in Dyad)
# PP1_dyad = 95.7            # [uM] - Total dyadic [PP1]
# PP1_SL = 0.57              # [uM] - Total Subsarcolemmal [PP1]
# PP2A_dyad = 95.76          # [uM] - Total dyadic PP2A
# OA = 0                     # [uM] - PP1/PP2A inhibitor Okadaic Acid
# PLBtot = plb_val           # [uM] - Total [PLB] in cytosolic units

# Parameters for BAR module
# Ligtot = 0.1               # [uM] - SET LIGAND CONCENTRATION (0 or 0.1)
# LCCtotBA = 0.025           # [uM] - [umol/L cytosol]
# RyRtotBA = 0.135           # [uM] - [umol/L cytosol]
# PLBtotBA = plb_val         # [uM] - [umol/L cytosol]
# TnItotBA = 70              # [uM] - [umol/L cytosol]
# IKstotBA = 0.025           # [uM] - [umol/L cytosol]
# ICFTRtotBA = 0.025         # [uM] - [umol/L cytosol]
# PP1_PLBtot = 0.89          # [uM] - [umol/L cytosol]
# IKurtotBA = 0.025          # [uM] - [umol/L cytosol] MOUSE
# PLMtotBA = 48              # [uM] - [umol/L cytosol] MOUSE

# For Recovery from inactivation of LTCC
# recoveryTime = 10 # initialize to smallest value

# Parameter varied in protocol simulation
# variablePar = 20 # initilization



## This part is combining all ODEs scattered in different modules
function combined_ode!(du, u, p, t)
    # Create separate instances of Module_CAM
    # Module_CAM_instance1 = Module_CAM.initialize()
    # Module_CAM_instance2 = Module_CAM.initialize()
    # Module_CAM_instance3 = Module_CAM.initialize()

## Collect ODEs
    # ODEs from Module_ECC // y_ecc
    Module_ECC.morotti_eccODE(du[1:87], u[1:87], pECC, LCC_CKdyadp, RyR_CKp, PLB_CKp, t)
    # ODEs from Module_camDyad 
    Module_camDyad.morotti_camODE(du[87+1:87+15], u[87+1:87+15], pCaMDyad, CaDyad, t)
    # ODEs from Module_camSL
    Module_camSL.morotti_camODE(du[87+15+1:87+15+15], u[87+15+1:87+15+15], pCaMSL, CaSL, t)
    # ODEs from Module_camCyt
    Module_camCyt.morotti_camODE(du[87+30+1:87+30+15], u[87+30+1:87+30+15], pCaMCyt, CaCyt, t)
    # ODEs from Module_CAMKII // y_camkii
    Module_CAM.morotti_camkiiODE(du[87+45+1:87+45+6], u[87+45+1:87+45+6], pCaMKII, t)
    # ODEs from Module_BAR // y_BAR
    Module_BAR.morotti_barODE(du[87+45+6+1:87+45+6+36], u[87+45+6+1:87+45+6+36], pBAR, t)

    ## Define Parameters
    pECC = [cycleLength, LCCa_PKAp, LCCb_PKAp, PLB_PKAn, RyR_PKAp, TnI_PKAp, IKs_PKAp, ICFTR_PKAp, 
            CKIIOE, recoveryTime, variablePar, Ligtot, IKur_PKAp, PLM_PKAp]
    pCaMDyad = [K, Mg, CaMtotDyad, 0, CaMKIItotDyad, CaNtotDyad, PP1totDyad, CaDyad, cycleLength, compart_dyad]
    pCaMSL = [K, Mg, CaMtotSL, BtotSL, CaMKIItotSL, CaNtotSL, PP1totSL, CaSL, cycleLength,compartSL] 
    pCaMCyt = [K, Mg, CaMtotCyt, BtotCyt, CaMKIItotCyt, CaNtotCyt, PP1totCyt, CaCyt, cycleLength,compartCyt]
    pCaMKII = [CaMKIIact_Dyad,LCCtotDyad,RyRtot,PP1_dyad,PP2A_dyad,OA,PLBtot,CaMKIIact_SL,LCCtotSL,PP1_SL,PP1_PLB_avail]
    pBAR = [Ligtot, LCCtotBA, RyRtotBA, PLBtotBA, TnItotBA, IKstotBA, ICFTRtotBA, PP1_PLBtot, IKurtotBA, PLMtotBA]
    


## Calculate the value of parameters based on variables of other module 

    # CaM module
    @variables CaDyad(t) CaDyad ~ u[36]*1e3      # from ECC model, *** Converting from [mM] to [uM] ***
    @variables CaSL(t) CaSL ~ u[37]*1e3          # from ECC model, *** Converting from [mM] to [uM] ***
    @variables CaCyt(t) CaCyt ~ u[38]*1e3        # from ECC model, *** Converting from [mM] to [uM] ***

    # CaMKII phosphorylation module 
    @variables CaMKIIact_Dyad(t) CaMKIIact_Dyad ~ CaMKIItotDyad.*sum(u[87+8:87+11]) # Multiply total by fraction
    @variables CaMKIIact_SL(t) CaMKIIact_SL ~ CaMKIItotSL.*sum(u[87+15+8:87+15+11])
    #PP1_PLB_avail = y(87+45+6+22)./PP1_PLBtot + .0091  # Active PP1 near PLB / total PP1 conc + basal value
    @variables PP1_PLB_avail(t) PP1_PLB_avail ~ 1 - u[87+45+6+24]/PP1_PLBtot + 0.081698  # NEW ODEs
    # PP1_PLB_avail is 100# (1) with NO ISO (Derived)
    @variables LCC_CKdyadp(t) LCC_CKdyadp ~ u[87+45+2] ./ LCCtotDyad     # fractional CaMKII-dependent LCC dyad phosphorylation
    @variables RyR_CKp(t) RyR_CKp ~ u[87+45+4] ./ RyRtot                 # fractional CaMKII-dependent RyR phosphorylation
    @variables PLB_CKp(t) PLB_CKp ~ u[87+45+5] ./ PLBtot                 # fractional CaMKII-dependent PLB phosphorylation

    # BAR (PKA phosphorylation) module
    @variables LCCa_PKAp(t) LCCa_PKAp ~ u[87+45+6+28]./LCCtotBA
    @variables LCCb_PKAp(t) LCCb_PKAp ~ u[87+45+6+29]./LCCtotBA
    @variables PLB_PKAn(t) PLB_PKAn ~ (PLBtotBA - u[87+45+6+26])./PLBtotBA
    @variables RyR_PKAp(t) RyR_PKAp ~ u[87+45+6+30]./RyRtotBA
    @variables TnI_PKAp(t) TnI_PKAp ~ u[87+45+6+31]./TnItotBA
    @variables IKs_PKAp(t) IKs_PKAp ~ u[87+45+6+34]./IKstotBA
    @variables ICFTR_PKAp(t) ICFTR_PKAp ~ u[87+45+6+35]./ICFTRtotBA
    @variables IKur_PKAp(t) IKur_PKAp ~ u[87+45+6+36]./IKurtotBA
    # PLM_PKAn = (PLMtotBA - y(87+45+6+32))./PLMtotBA
    @variables PLM_PKAp(t) PLM_PKAp ~ u[87+45+6+27]./PLMtotBA

    # incorporate Ca buffering from CaM, convert JCaCyt from uM/msec to mM/msec
    du[35] = du[36] + 1e-3*JCaDyad           #  incorporate Ca buffering from CaM, convert JCaDyad from uM/msec to mM/msec
    du[37] = du[37] + 1e-3*JCaSL             #  incorporate Ca buffering from CaM, convert JCaSL from uM/msec to mM/msec
    du[38] = du[38] + 1e-3*JCaCyt            #  incorporate Ca buffering from CaM, convert JCaCyt from uM/msec to mM/msec

    # incorporate CaM diffusion between compartments
    Vmyo = 2.1454e-11       # [L]
    Vdyad = 1.7790e-014     # [L]
    VSL = 6.6013e-013       # [L]
    # kDyadSL = 3.6363e-16	# [L/msec]
    kSLmyo = 8.587e-15      # [L/msec]
    k0Boff = 0.0014         # [s^-1] 
    k0Bon = k0Boff/0.2      # [uM^-1 s^-1] kon = koff/Kd
    k2Boff = k0Boff/100     # [s^-1] 
    k2Bon = k0Bon           # [uM^-1 s^-1]
    # k4Boff = k2Boff       # [s^-1]
    k4Bon = k0Bon           # [uM^-1 s^-1]
    @variables CaMtotDyad(t) CaMtotDyad ~ sum(u[87+1:87+6])+CaMKIItotDyad*sum(u[87+7:87+10])+sum(u[87+13:87+15])
    @variables Bdyad(t) Bdyad ~ BtotDyad - CaMtotDyad # [uM dyad]
    @variables J_cam_dyadSL(t) J_cam_dyadSL ~ 1e-3*(k0Boff*u[87+1] - k0Bon*Bdyad*u[87+15+1]) # [uM/msec dyad]
    @variables J_ca2cam_dyadSL(t) J_ca2cam_dyadSL ~ 1e-3*(k2Boff*u[87+2] - k2Bon*Bdyad*u[87+15+2]) # [uM/msec dyad]
    @variables J_ca4cam_dyadSL(t) J_ca4cam_dyadSL ~ 1e-3*(k2Boff*u[87+3] - k4Bon*Bdyad*u[87+15+3]) # [uM/msec dyad]
    @variables J_cam_SLmyo(t) J_cam_SLmyo ~ kSLmyo*(u[87+15+1]-u[87+30+1]) # [umol/msec]
    @variables J_ca2cam_SLmyo(t) J_ca2cam_SLmyo ~ kSLmyo*(u[87+15+2]-u[87+30+2]) # [umol/msec]
    @variables J_ca4cam_SLmyo(t) J_ca4cam_SLmyo  ~ kSLmyo*(u[87+15+3]-u[87+30+3]) # [umol/msec]
    du[87+1] = du[87+1] - J_cam_dyadSL
    du[87+2] = du[87+2] - J_ca2cam_dyadSL
    du[87+3] = du[87+3] - J_ca4cam_dyadSL
    du[87+15+1] = du[87+15+1] + J_cam_dyadSL*Vdyad/VSL - J_cam_SLmyo/VSL
    du[87+15+2] = du[87+15+2] + J_ca2cam_dyadSL*Vdyad/VSL - J_ca2cam_SLmyo/VSL
    du[87+15+3] = du[87+15+3] + J_ca4cam_dyadSL*Vdyad/VSL - J_ca4cam_SLmyo/VSL
    du[87+30+1] = du[87+30+1] + J_cam_SLmyo/Vmyo
    du[87+30+2] = du[87+30+2] + J_ca2cam_SLmyo/Vmyo
    du[87+30+3] = du[87+30+3] + J_ca4cam_SLmyo/Vmyo

    return du
end


## Collect all parameters and define mass matrix for BAR module

p = [cycleLength => 1e3/freq, recoveryTime => 10, variablePar => 20, CaMtotDyad => 418, BtotDyad => 1.54/8.293e-4, CaMKIItotDyad => 120,
     CaNtotDyad => 3e-3/8.293e-4, PP1totDyad => 96.5, CaMtotSL => 5.65, BtotSL => 24.2, CaMKIItotSL => 120*8.293e-4, CaNtotSL => 3e-3,
     PP1totSL => 0.57, CaMtotCyt => 5.65, BtotCyt => 24.2, CaMKIItotCyt => 120*8.293e-4, CaNtotCyt => 3e-3, PP1totCyt => 0.57,
     LCCtotDyad => 31.4*0.9, RyRtot => 382.6, PP1_dyad => 95.7, PP2A_dyad => 95.76, OA => 0, PLBtot => plb_val, LCCtotSL => 0.0846,
     PP1_SL => 0.57, Ligtot => 0.1, LCCtotBA => 0.025, RyRtotBA => 0.135, PLBtotBA => plb_val, TnItotBA => 70, IKstotBA => 0.025, 
     ICFTRtotBA => 0.025, PP1_PLBtot => 0.89, IKurtotBA => 0.025, PLMtotBA => 48, CKIIOE => tempCKIIOE]

u0 = readxlsheet("Morotti_WT_1Hz.xls", "Sheet1")

tspan = (0, 3e3)


# Solve the ODE
prob = ODEProblem(combined_ode, u0, tspan, p)
sol = solve(prob, Tsit5())

plot(sol.t, (@. sol[38,:]), title="Morotti's Model", xlabel="Time(s)", label="[Ca2+]i")

