# the ODE's for CaM, CaMKII, and CaN.

using Catalyst
using ModelingToolkit

function make_cam_rn(;
    Mg=1.0,  # Magnesium (mM)
    K=135.0, # Potassium (mM)
    CaMKII=120.0, # Dyad: 120; SL:120*8.293e-4; Cyt:120*8.293e-4
    PP1tot=96.5, # Dyad: 96.5; SL:0.57; Cyt:0.57
)
    function kbt(Pb, Pt, Pt2, Pa)
        T = Pb + Pt + Pt2 + Pa
        return 0.055*T + 0.0074*T^2 + 0.015*T^3
    end

    rn = @reaction_network begin
        @variables Ca(t)

        # CaM: Calmodulin
        # CaN: Calcineurin
        # B: buffer
        (k02 * Ca^2, k20), CaM <--> Ca2CaM
        (k24 * Ca^2, k42), Ca2CaM <--> Ca4CaM
        (k02B * Ca^2, k20B), CaMB <--> Ca2CaMB
        (k24B * Ca^2, k42B), Ca2CaMB <--> Ca4CaMB
        (k0Bon, k0Boff), CaM + B <--> CaMB
        (k2Bon, k2Boff), Ca2CaM + B <--> Ca2CaMB
        (k4Bon, k4Boff), Ca4CaM + B <--> Ca4CaMB
        (kcanCaon * Ca^2, kcanCaoff), Ca2CaN <--> Ca4CaN
        (k02can * Ca^2, k20can), CaMCa4CaN <--> Ca2CaMCa4CaN
        (k24can * Ca^2, k42can), Ca2CaMCa4CaN <--> Ca4CaMCa4CaN
        (kcanCaM0on, kcanCaM0off), CaM + Ca4CaN <--> CaMCa4CaN
        (kcanCaM2on, kcanCaM2off), Ca2CaM + Ca4CaN <--> Ca2CaMCa4CaN
        (kcanCaM4on, kcanCaM4off), Ca4CaM + Ca4CaN <--> Ca4CaMCa4CaN
        # CaMKII
        (kib2, kb2i), Ca2CaM + Pi <--> Pb2
        (kb24 * Ca^2, kb42), Pb2 <--> Pb
        (kib, kbi), Ca4CaM + Pi <--> Pb
        kbt(Pb, Pt, Pt2, Pa), Pb --> Pt
        kpp1*PP1tot*Pt/(Kmpp1+CaMKIItot*Pt), Pt => Pb
        (kt42 * Ca^2, kt24), Pt <--> Pt2
        (kta, kat), Pt <--> Ca4CaM + Pa
        (kt2a, kat2), Pt2 <--> Ca2CaM + Pa
        kpp1*PP1tot*Pt2/(Kmpp1+CaMKIItot*Pt2), Pt2 => Pb2
        kpp1*PP1tot*Pa/(Kmpp1+CaMKIItot*Pa), Pa => Pi
    end

    # Calculate Ca/CaM parameters
    if Mg <= 1
        Kd02 = 0.0025*(1+K/0.94-Mg/0.012)*(1+K/8.1+Mg/0.022)  # [uM^2]
        Kd24 = 0.128*(1+K/0.64+Mg/0.0014)*(1+K/13.0-Mg/0.153) # [uM^2]
    else
        Kd02 = 0.0025*(1+K/0.94-1/0.012+(Mg-1)/0.060)*(1+K/8.1+1/0.022+(Mg-1)/0.068)   # [uM^2]
        Kd24 = 0.128*(1+K/0.64+1/0.0014+(Mg-1)/0.005)*(1+K/13.0-1/0.153+(Mg-1)/0.150)  # [uM^2]
    end

    k20 = 10*1e-3          # [ms^-1]
    k02 = k20/Kd02         # [uM^-2 ms^-1]
    k42 = 500*1e-3         # [ms^-1]
    k24 = k42/Kd24         # [uM^-2 ms^-1]

    # CaM buffering (B) parameters
    k0Boff = 0.0014*1e-3   # [ms^-1]
    k0Bon = k0Boff/0.2     # [uM^-1 ms^-1] kon = koff/Kd
    k2Boff = k0Boff/100    # [ms^-1]
    k2Bon = k0Bon          # [uM^-1 ms^-1]
    k4Boff = k2Boff        # [ms^-1]
    k4Bon = k0Bon          # [uM^-1 ms^-1]

    # using thermodynamic constraints
    k20B = k20/100 # [ms^-1] thermo constraint on loop 1
    k02B = k02     # [uM^-2 ms^-1]
    k42B = k42     # [ms^-1] thermo constraint on loop 2
    k24B = k24     # [uM^-2 ms^-1]

    # CaMKII parameters
    # Wi Wa Wt Wp
    kbi = 2.2*1e-3      # [ms^-1] (Ca4CaM dissocation from Wb)
    kib = kbi/33.5e-3   # [uM^-1 ms^-1]
    kib2 = kib
    kb2i = kib2*5
    kb24 = k24
    kb42 = k42*33.5e-3/5
    kpp1 = 1.72*1e-3    # [ms^-1] (PP1-dep dephosphorylation rates)
    Kmpp1 = 11.5        # [uM]
    kta = kbi/1000      # [ms^-1] (Ca4CaM dissociation from Wt)
    kat = kib           # [uM^-1 ms^-1] (Ca4CaM reassociation with Wa)
    kt42 = k42*33.5e-6/5
    kt24 = k24
    kat2 = kib
    kt2a = kib*5

    # CaN parameters
    kcanCaoff = 1e-3           # [ms^-1]
    kcanCaon = kcanCaoff/0.5   # [uM^-1 ms^-1]
    kcanCaM4on = 46*1e-3       # [uM^-1 ms^-1]
    kcanCaM4off = 1.3e-6       # [ms^-1]
    kcanCaM2on = kcanCaM4on
    kcanCaM2off = 2508*kcanCaM4off
    kcanCaM0on = kcanCaM4on
    kcanCaM0off = 165*kcanCaM2off
    k02can = k02
    k20can = k20/165
    k24can = k24
    k42can = k20/2508

    setdefaults!(rn, [
        :k20 => k20,
        :k02 => k02,
        :k42 => k42,
        :k24 => k24,
        :k0Boff => k0Boff,
        :k0Bon => k0Bon,
        :k2Boff => k2Boff,
        :k2Bon => k2Bon,
        :k4Boff => k4Boff,
        :k4Bon => k4Bon,
        :k20B => k20B,
        :k02B => k02B,
        :k42B => k42B,
        :k24B => k24B,
        :kbi => kbi,
        :kib => kib,
        :kib2 => kib2,
        :kb2i => kb2i,
        :kb24 => kb24,
        :kb42 => kb42,
        :kpp1 => kpp1,
        :Kmpp1 => Kmpp1,
        :kta => kta,
        :kat => kat,
        :kt42 => kt42,
        :kt24 => kt24,
        :kat2 => kat2,
        :kt2a => kt2a,
        :kcanCaoff => kcanCaoff,
        :kcanCaon => kcanCaon,
        :kcanCaM4on => kcanCaM4on,
        :kcanCaM4off => kcanCaM4off,
        :kcanCaM2on => kcanCaM2on,
        :kcanCaM2off => kcanCaM2off,
        :kcanCaM0on => kcanCaM0on,
        :kcanCaM0off => kcanCaM0off,
        :k02can => k02can,
        :k20can => k20can,
        :k24can => k24can,
        :k42can => k42can,
        :CaMKIItot => CaMKII
    ])

    return rn
end

# adjusting Ca buffering in EC coupling model
# I don'treally want to repeat this but I'm still waiting for options
function ca_buf_cam_eq(sys)
    @unpack Ca, k02, k20, k24, k42, CaM, Ca2CaM, Ca4CaM = sys
    rcn02 = k02*Ca^2*CaM - k20*Ca2CaM
    rcn24 = k24*Ca^2*Ca2CaM - k42*Ca4CaM

    @unpack k02B, k20B, k24B, k42B, CaMB, Ca2CaMB, Ca4CaMB = sys
    rcn02B = k02B*Ca^2*CaMB - k20B*Ca2CaMB
    rcn24B = k24B*Ca^2*Ca2CaMB - k42B*Ca4CaMB

    @unpack kcanCaon, kcanCaoff, Ca2CaN, Ca4CaN = sys
    rcnCa4CaN = kcanCaon*Ca^2*Ca2CaN - kcanCaoff*Ca4CaN

    @unpack k02can, k20can, k24can, k42can, CaMCa4CaN, Ca2CaMCa4CaN, Ca4CaMCa4CaN = sys
    rcn02CaN = k02can*Ca^2*CaMCa4CaN - k20can*Ca2CaMCa4CaN
    rcn24CaN = k24can*Ca^2*Ca2CaMCa4CaN - k42can*Ca4CaMCa4CaN

    @unpack kb24, kb42, Pb2, Pb, kt42, kt24, Pt, Pt2 = sys
    rcnCKb2b = kb24*Ca^2*Pb2 - kb42*Pb
    rcnCKtt2 = kt42*Pt - kt24*Ca^2*Pt2

    @unpack CaMKIItot = sys

    return -2*(CaMKIItot*(-rcnCKtt2+rcnCKb2b) - (rcn02+rcn24+rcn02B+rcn24B+rcnCa4CaN+rcn02CaN+rcn24CaN))
end


cam_rn = make_cam_rn()

ca_buf_cam_eq(cam_rn)
