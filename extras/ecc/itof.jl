# I_to: Transient Outward K Current (slow and fast components)
# Only the fast one (I_tof) is present on mouse

using ModelingToolkit

function get_itof_eqs(vm, ek, Kcoeff=1, CKIIflag=0)
    xtofs = expit((vm+3)/13)
    ytofs = expit(-(vm+48)/5)
    tauxtof = 0.08+0.7*exp(-((vm+25)/30)^2)
    tauytof = 10 - 5 * (CKIIflag) +32*exp(-((vm+55)/16)^2)+(8 + 4 * CKIIflag)/(1+exp(-(vm+60)/8))

    @parameters GtoFast = 0.44 * (1 - 1/3*CKIIflag) # [mS/uF] changed from rabbit (0.02)
    @variables t iTofx(t) iTofy(t) I_tof(t) I_to(t)
    D = Differential(t)

    return [
        D(iTofx) ~ (xtofs - iTofx) / tauxtof
        D(iTofy) ~ (ytofs - iTofy) / tauytof
        I_tof ~ Kcoeff * GtoFast * iTofx * iTofy * (vm - ek)
    ]
end
