# # Effects of pacing durations
# Experiments vs simulations.
using ModelingToolkit
using OrdinaryDiffEq, SteadyStateDiffEq, DiffEqCallbacks
using Plots
using CSV
using DataFrames
using CurveFit
using CaMKIIModel
using CaMKIIModel: second
Plots.default(lw=1.5)
