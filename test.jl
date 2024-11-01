using ProgressLogging
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM, get_camkii_fast_ca_binding_sys
Plots.default(lw=2)

sys = get_camkii_fast_ca_binding_sys(simplify=true)

using Latexify

write("eqs.md", latexify(equations(sys)))
