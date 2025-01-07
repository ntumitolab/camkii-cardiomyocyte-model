# # Pacing response
using ModelingToolkit
using Latexify
using CaMKIIModel
using CaMKIIModel: get_ik_sys, get_ica_sys

@independent_variables t
@variables Vm(t)
@parameters IKUR_PKAp LCCb_PKAp k_i k_o na_i na_o ca_i ca_o

sys = get_ik_sys(k_i, k_o, na_i, na_o, V_m; IKUR_PKAp)

Latexify.copy_to_clipboard(true)

latexify( get_ica_sys(na_i, ca_i, na_o, ca_o, Vm; LCCb_PKAp=0), env=:inline, fmt=FancyNumberFormatter())
