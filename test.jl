using OrdinaryDiffEq
using DiffEqCallbacks
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots
using CaMKIIModel: build_neonatal_ecc_sys

sys = build_neonatal_ecc_sys()

simp = structural_simplify(sys)
