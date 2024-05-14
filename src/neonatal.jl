# Model of ECC of rat neonatal ventricular myocyte 2009
# Korhonen et al. "Model of excitation-contraction coupling of rat neonatal
# ventricular myocytes" Biophys J. 2009, Feb; 96(3):1189-1209
using ModelingToolkit

function build_neonatal_model(;
    dx=0.1,  # [um]
    rSR_true=6, # [um]
    rSL_true=10.5, # [um]
)

    @constants begin
        Acap = 4 * pi * rSL_true^2 * 1e-8 # cm^2
        VSR = 0.043 * 1.5 * 1.4 # pL
        VNSR = 0.9 * VSR # pL
        VJSR = 0.1 * VSR # pL
        V_sub_SR = 4 / 3 * pi * (rSR_true + dx)^3 / 1000 - 4 / 3 * pi * (rSR_true)^3 / 1000 # pL
        V_sub_SL = 4 / 3 * pi * rSL_true^3 / 1000 - 4 / 3 * pi * (rSL_true - dx)^3 / 1000 # pL
        Vmyo = 4 / 3 * pi * rSL_true^3 / 1000 - 4 / 3 * pi * rSR_true^3 / 1000 # pL
        ACAP_VMYO_F = Acap / F * 1e6 / Vmyo
    end

    @parameters begin
        Cao = 1796   # [uM]
        Nao = 154578 # [uM]
        Ko = 5366    # [uM]
        Mg = 1000    # [uM]
        Ki = 135000  # [uM]
    end

end
