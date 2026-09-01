#===
# Figure 12: Temporal pattern discrimination by CaMKII
#
# Part of the CaMKII temporal-encoding analysis. Run order:
#   1. temporal_encoding_utils.jl   (helpers; included automatically)
#   2. figure11_temporal_map.jl
#   3. figure12_temporal_pattern.jl (this file)
#
# Panels produced:
#   A  WT active CaMKII traces          4 protocols, each delivering 60 stimuli
#   B  R275H-like traces                same protocols, y-axis shared with A
#   C  Peak active CaMKII fraction      bars, WT vs mutant x 4 protocols
#   D  Active CaMKII 30 s after last pulse ("memory")
#
# Supplementary figure: AUC bars, tau bars, and the frequency-vs-ROS comparison.
#
# ## The comparison that matters
# Protocols 1 and 4 are the tight test: BOTH are 1 Hz and BOTH deliver 60 beats,
# differing only in whether the beats are continuous or split by a 30 s pause.
# Protocols 2 and 3 vary frequency as well, so they mainly re-test the frequency
# dependence already established in Figure 11. The script prints the continuous-
# vs-split comparison explicitly for both genotypes.
#
# ## AUC window
# The protocols end their final pulse at different wall-clock times (60 s vs 90 s
# from the first pulse), so peak and AUC are computed over ONE shared absolute
# window (same start, same length, set by the longest protocol). Integrating each
# protocol from its own last pulse would give the split protocol a 30 s longer
# window and inflate its AUC independent of any real effect. tau and the 30 s
# memory readout stay referenced to each protocol's own last pulse, which is the
# correct comparison for those.
#
# Only 10 ODE solves; runs in under a minute.
===#
using Model
using Model: second, Hz, μM
using CSV
using DataFrames
using CurveFit
using DiffEqCallbacks
using ModelingToolkit
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using Plots
using StatsPlots
Plots.default(lw=2)

include(joinpath(@__DIR__, "temporal_encoding_utils.jl"))

# ## Configuration
const OUT = joinpath(@__DIR__, "figure12")
const PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
const PALETTE_SUP = ["#4a3aa7", "#e34948"]
const WT_COLOR, MUT_COLOR = "#2a78d6", "#eb6834"
const GENOTYPE_LABELS = Dict(:WT => "WT", :R275H => "R275H-like")

sys = Model.DEFAULT_SYS
alg = KenCarp47()
@unpack Istim = sys
const STIMSTART = 30.0second
const POSTSTIM_WINDOW = 90.0
const MEMORY_DELAY_S = 30.0

# IMPORTANT: `build_stim_callbacks` treats both starttime and endtime as pulse times
# (`starttime:period:endtime` is inclusive). Therefore the endpoint for N stimuli is
# start + (N-1)*period, not start + N*period. These definitions deliver EXACTLY
# 60 stimuli in every protocol.
const PROTOCOLS = [
    "60 stimuli @ 1 Hz" => [(STIMSTART, STIMSTART + 59.0second, 1.0second)],
    "60 stimuli @ 2 Hz" => [(STIMSTART, STIMSTART + 29.5second, 0.5second)],
    "60 stimuli @ 0.5 Hz" => [(STIMSTART, STIMSTART + 118.0second, 2.0second)],
    # Pulses 0...29 s, 30 s without a stimulus (30...59 s), then pulses 60...89 s.
    "30+30 stimuli @ 1 Hz, 30 s gap" => [(STIMSTART, STIMSTART + 29.0second, 1.0second),
        (STIMSTART + 60.0second, STIMSTART + 89.0second, 1.0second)],
]

const PROTOCOL_SHORT = Dict(
    "60 stimuli @ 1 Hz" => "1 Hz\ncont.",
    "60 stimuli @ 2 Hz" => "2 Hz\ncont.",
    "60 stimuli @ 0.5 Hz" => "0.5 Hz\ncont.",
    "30+30 stimuli @ 1 Hz, 30 s gap" => "1 Hz\nsplit",
)

const CONTINUOUS = "60 stimuli @ 1 Hz"
const SPLIT = "30+30 stimuli @ 1 Hz, 30 s gap"

const MAX_STIMEND = maximum(maximum(b[2] for b in blocks) for (_, blocks) in PROTOCOLS)
const COMMON_TEND_S = (MAX_STIMEND - STIMSTART) / second + POSTSTIM_WINDOW
println("Common recording window for peak/AUC: $(COMMON_TEND_S) s from first pulse.")

prob0 = ODEProblem(sys, [], (0.0, STIMSTART + COMMON_TEND_S * second + 10second))

function run_protocol(blocks; mutant=false, iso=0.0μM, ros=0.0μM, record_tend_s=nothing)
    stimend = maximum(b[2] for b in blocks)
    default_tend_s = (stimend - STIMSTART) / second + POSTSTIM_WINDOW
    tend_s = record_tend_s === nothing ? default_tend_s : max(record_tend_s, default_tend_s)
    tend = STIMSTART + tend_s * second + 5second
    cb = build_schedule_callback(Istim, blocks)
    prob = build_condition_prob(prob0, sys; mutant, iso, ros, tend)
    sol = solve(prob, alg; callback=cb, reltol=1e-8, abstol=1e-8)
    summ = summarize_camkii(sol, sys, STIMSTART, STIMSTART + tend_s * second, 0.0)
    tau = get_decay_tau(sol, sys, stimend; window=min(50.0, POSTSTIM_WINDOW))
    memory = sol(stimend + MEMORY_DELAY_S * second, idxs=sys.CaMKAct)
    return (; sol, stimend, summ, tau, memory)
end

# ## Run all protocols x both genotypes
rows = NamedTuple[]
sols = Dict{Tuple{String,Symbol},Any}()
for (name, blocks) in PROTOCOLS, (gname, is_mut) in pairs((WT=false, R275H=true))
    out = run_protocol(blocks; mutant=is_mut, record_tend_s=COMMON_TEND_S)
    push!(rows, (; protocol=name, genotype=gname, peak=out.summ.peak, auc=out.summ.auc,
        tau=out.tau, memory_30s=out.memory))
    sols[(name, gname)] = (out.sol, out.stimend)
    println("$name / $gname: peak=$(round(out.summ.peak, digits=3)), auc=$(round(out.summ.auc, digits=2)), tau=$(round(out.tau, digits=2)) s, memory@30s=$(round(out.memory, digits=4))")
end
results = DataFrame(rows)
CSV.write(OUT * "-results.csv", results)

const NAMES_ORDER = first.(PROTOCOLS)
const SHORT_LABELS = [PROTOCOL_SHORT[n] for n in NAMES_ORDER]
const YMAX = maximum(results.peak) * 1.08

# ## The matched-frequency, matched-beat-count comparison
println("\n=== Continuous vs split (both 1 Hz, both exactly 60 stimuli) ===")
for g in (:WT, :R275H)
    c = only(results[(results.protocol.==CONTINUOUS).&(results.genotype.==g), :])
    s = only(results[(results.protocol.==SPLIT).&(results.genotype.==g), :])
    pct(a, b) = round(100 * (a - b) / b, digits=1)
    println("$(GENOTYPE_LABELS[g]):")
    println("  peak      $(round(c.peak, digits=3)) vs $(round(s.peak, digits=3))   ($(pct(c.peak, s.peak))%)")
    println("  memory30  $(round(c.memory_30s, digits=4)) vs $(round(s.memory_30s, digits=4))   ($(pct(c.memory_30s, s.memory_30s))%)")
    println("  tau       $(round(c.tau, digits=2)) vs $(round(s.tau, digits=2)) s   ($(pct(c.tau, s.tau))%)")
    println("  AUC       $(round(c.auc, digits=2)) vs $(round(s.auc, digits=2))   ($(pct(c.auc, s.auc))%)")
end

# ## Panels A and B: traces, one genotype each, shared y-limits
function trace_panel(gname, letter; showlegend=false)
    p = plot(xlabel="Time relative to first pulse (s)", ylabel="Active CaMKII fraction",
        title=letter * "   " * GENOTYPE_LABELS[gname], titlelocation=:left, titlefontsize=12,
        guidefontsize=10, tickfontsize=9, ylims=(0, YMAX),
        legend=showlegend ? :outertopright : false, legendfontsize=8,
        left_margin=10Plots.mm, bottom_margin=9Plots.mm, top_margin=5Plots.mm)
    for (i, (name, _)) in enumerate(PROTOCOLS)
        sol, stimend = sols[(name, gname)]
        plot!(p, sol, idxs=((sys.t - STIMSTART) / 1000, sys.CaMKAct),
            tspan=(STIMSTART, stimend + POSTSTIM_WINDOW * second), lab=name, color=PALETTE[i])
    end
    return p
end

# ## Panels C and D: bars
function bar_panel(col, ylab, letter; showlegend=false)
    wt = [only(results[(results.protocol.==n).&(results.genotype.==:WT), col]) for n in NAMES_ORDER]
    mut = [only(results[(results.protocol.==n).&(results.genotype.==:R275H), col]) for n in NAMES_ORDER]
    return groupedbar(SHORT_LABELS, [wt mut], label=["WT" "R275H-like"], color=[WT_COLOR MUT_COLOR],
        ylabel=ylab, title=letter * "   " * ylab, titlelocation=:left, titlefontsize=12,
        legend=showlegend ? :topleft : false,
        guidefontsize=10, tickfontsize=9, legendfontsize=9,
        left_margin=10Plots.mm, bottom_margin=9Plots.mm, top_margin=5Plots.mm)
end

panelA = trace_panel(:WT, "A"; showlegend=true)
panelB = trace_panel(:R275H, "B")
panelC = bar_panel(:peak, "Peak active CaMKII fraction", "C"; showlegend=true)
panelD = bar_panel(:memory_30s, "Active CaMKII, 30 s after last pulse", "D")

for (p, tag, sz) in ((panelA, "A-wt-traces", (1050, 520)), (panelB, "B-mutant-traces", (1050, 520)),
    (panelC, "C-peak", (900, 520)), (panelD, "D-memory", (900, 520)))
    savefig(plot(p, size=sz), OUT * "-" * tag * ".png")
    savefig(plot(p, size=sz), OUT * "-" * tag * ".pdf")
end

composite = plot(panelA, panelB, panelC, panelD, layout=(2, 2), size=(1900, 1100),
    left_margin=8Plots.mm, right_margin=8Plots.mm, bottom_margin=8Plots.mm, top_margin=5Plots.mm)
savefig(composite, OUT * "-composite.png"); savefig(composite, OUT * "-composite.pdf")

# ## Supplementary: AUC, tau, and the frequency-vs-ROS dissociation
# Note on the frequency/ROS pair: these were chosen to attempt matched peaks, but
# peaks are NOT matched in practice (0.64 vs 0.26) because frequency dominates peak
# far more strongly than ROS does. Read it as a dissociation - peak tracks
# frequency, tau tracks ROS - rather than a matched-peak test.
supA = bar_panel(:auc, "AUC, fixed $(round(Int, COMMON_TEND_S))-s window", "S1"; showlegend=true)
supB = bar_panel(:tau, "Post-pacing decay \u03c4 (s)", "S2")

combo = [
    "2 Hz, 60 s, low ROS (25 μM)" => (blocks=[(STIMSTART, STIMSTART + 60.0second, 0.5second)], ros=25.0μM),
    "0.5 Hz, 60 s, high ROS (100 μM)" => (blocks=[(STIMSTART, STIMSTART + 60.0second, 2.0second)], ros=100.0μM),
]
rows_sup = NamedTuple[]
supC = plot(xlabel="Time relative to first pulse (s)", ylabel="Active CaMKII fraction",
    title="S3   Frequency vs. ROS", titlelocation=:left, titlefontsize=12,
    guidefontsize=10, tickfontsize=9, legend=:outertopright, legendfontsize=8,
    left_margin=10Plots.mm, bottom_margin=9Plots.mm, top_margin=5Plots.mm)
for (i, (name, cfg)) in enumerate(combo)
    out = run_protocol(cfg.blocks; mutant=false, ros=cfg.ros)
    push!(rows_sup, (; protocol=name, peak=out.summ.peak, auc=out.summ.auc, tau=out.tau, memory_30s=out.memory))
    plot!(supC, out.sol, idxs=((sys.t - STIMSTART) / 1000, sys.CaMKAct),
        tspan=(STIMSTART, out.stimend + POSTSTIM_WINDOW * second), lab=name, color=PALETTE_SUP[i])
end
CSV.write(OUT * "-supplementary-freq-vs-ros.csv", DataFrame(rows_sup))
println("\nFrequency vs ROS:\n", DataFrame(rows_sup))

savefig(plot(supA, supB, supC, layout=(1, 3), size=(1900, 520),
        left_margin=10Plots.mm, bottom_margin=10Plots.mm, top_margin=5Plots.mm),
    OUT * "-S-supplementary.png")
savefig(plot(supA, supB, supC, layout=(1, 3), size=(1900, 520),
        left_margin=10Plots.mm, bottom_margin=10Plots.mm, top_margin=5Plots.mm),
    OUT * "-S-supplementary.pdf")

println("\nWrote figure12-{A-wt-traces, B-mutant-traces, C-peak, D-memory, composite, S-supplementary}.png/.pdf")
println("Wrote figure12-results.csv, figure12-supplementary-freq-vs-ros.csv")
