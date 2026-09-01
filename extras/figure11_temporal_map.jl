#===
# Figure 11: CaMKII temporal response map
#
# Part of the CaMKII temporal-encoding analysis. Run order:
#   1. temporal_encoding_utils.jl   (helpers; included automatically)
#   2. figure11_temporal_map.jl     (this file)
#   3. figure12_temporal_pattern.jl
#
# Panels produced:
#   A  Peak active CaMKII fraction        heatmaps, 2 genotypes x 4 conditions
#   B  AUC of active CaMKII               heatmaps, 2 genotypes x 4 conditions
#   C  Decay tau vs. pacing duration      WT solid/filled, mutant dashed/open
#   D  Mutant/WT dynamic-range ratio (%)  peak, AUC, tau x 4 conditions
#
# Also written: a supplementary tau heatmap (not part of Figure 11), the
# dynamic-range summary CSV, and each panel as a standalone file.
#
# ## Simulation cache
# The 288-cell sweep is checkpointed to `temporal_map_results.csv`: each cell is
# appended as soon as it is solved, and any cell already present is skipped on a
# re-run. Deleting that file forces a full re-solve (tens of minutes); keeping it
# means this script reaches the plotting stage in seconds.
#
# ## Tau exclusion (important)
# All tau summaries exclude the 15 s pacing duration. Single-exponential fits to
# the short, low-amplitude decay following 15 s of pacing are poorly constrained
# and produced non-monotonic estimates (WT/Control at 0.5 Hz: 11.2 s after 15 s of
# pacing vs 20.7 s after 30 s). Every tau_min in the unfiltered summary came from
# that row, which inflated the apparent WT-vs-mutant difference. Peak and AUC are
# direct measurements rather than fits, so they use the full grid.
# `TAU_MIN_DURATION_S` controls the cutoff; the 15 s points are still PLOTTED in
# panel C, with a dotted rule marking where the statistics window begins.
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
Plots.default(lw=1.5)

include(joinpath(@__DIR__, "temporal_encoding_utils.jl"))

# ## Configuration
const SIM_CACHE = joinpath(@__DIR__, "temporal_map_results.csv")
const OUT = joinpath(@__DIR__, "figure11")
const TAU_MIN_DURATION_S = 30.0
const QUICK_TEST = false   # true = 3x3 grid smoke test; false = full 6x6 grid

const FREQS_HZ = QUICK_TEST ? [0.5, 1.0, 2.0] : [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
const DURATIONS_S = QUICK_TEST ? [15.0, 60.0, 120.0] : [15.0, 30.0, 60.0, 90.0, 120.0, 180.0]

const CONDITION_ORDER = (:Control, :ISO, :ROS, :ISO_ROS)
const CONDITION_LABELS = Dict(:Control => "Control", :ISO => "ISO (0.1 μM)",
    :ROS => "ROS (50 μM H₂O₂)", :ISO_ROS => "ISO + ROS")
const CONDITION_SHORT = Dict(:Control => "Control", :ISO => "ISO",
    :ROS => "ROS", :ISO_ROS => "ISO+ROS")
const GENOTYPE_ORDER = (:WT, :R275H)
const GENOTYPE_LABELS = Dict(:WT => "WT", :R275H => "R275H-like")

const WT_COLOR, MUT_COLOR, TAU_COLOR = "#2a78d6", "#eb6834", "#1baf7a"
const FREQ_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
const FREQ_MARKERS = [:circle, :square, :diamond, :utriangle]

const GENOTYPES = (WT=false, R275H=true)
const CONDITIONS = (
    Control=(iso=0.0μM, ros=0.0μM),
    ISO=(iso=0.1μM, ros=0.0μM),
    ROS=(iso=0.0μM, ros=50.0μM),
    ISO_ROS=(iso=0.1μM, ros=50.0μM),
)

# ## Sweep (skipped entirely when the cache is complete)
function run_sweep()
    sys = Model.DEFAULT_SYS
    alg = KenCarp47()
    @unpack Istim = sys
    stimstart = 30.0second
    poststim_window = 60.0

    tend_max = stimstart + maximum(DURATIONS_S) * second + poststim_window * second + 5second
    prob0 = ODEProblem(sys, [], (0.0, tend_max))

    done = Set{Tuple{Symbol,Symbol,Float64,Float64}}()
    if isfile(SIM_CACHE)
        cached = CSV.read(SIM_CACHE, DataFrame)
        for r in eachrow(cached)
            push!(done, (Symbol(r.genotype), Symbol(r.condition), Float64(r.freq_Hz), Float64(r.dur_s)))
        end
        println("Cache: $(nrow(cached)) rows already computed.")
    end

    n_new = n_fail = 0
    for (gname, is_mut) in pairs(GENOTYPES), (cname, cond) in pairs(CONDITIONS)
        for dur_s in DURATIONS_S, freq_Hz in FREQS_HZ
            (gname, cname, freq_Hz, dur_s) in done && continue
            try
                stimend = stimstart + dur_s * second
                tend = stimend + poststim_window * second + 5second
                cb = build_stim_callbacks(Istim, stimend; period=inv(freq_Hz * Hz), starttime=stimstart)
                prob = build_condition_prob(prob0, sys; mutant=is_mut, iso=cond.iso, ros=cond.ros, tend=tend)
                sol = solve(prob, alg; callback=cb, reltol=1e-8, abstol=1e-8)
                summ = summarize_camkii(sol, sys, stimstart, stimend, poststim_window)
                tau = get_decay_tau(sol, sys, stimend)
                row = (; genotype=gname, condition=cname, freq_Hz, dur_s,
                    peak=summ.peak, auc=summ.auc, tau=tau)
                CSV.write(SIM_CACHE, DataFrame([row]); append=isfile(SIM_CACHE))
                push!(done, (gname, cname, freq_Hz, dur_s))
                n_new += 1
            catch e
                ## An interrupt means the user (or the runtime) stopped the run.
                ## Never swallow it: rethrow so the script halts rather than
                ## silently producing an incomplete grid.
                e isa InterruptException && rethrow()
                n_fail += 1
                println("  FAILED $gname/$cname f=$freq_Hz d=$dur_s: ", sprint(showerror, e))
            end
        end
    end
    println("Sweep: $n_new new, $n_fail failed.")
    return CSV.read(SIM_CACHE, DataFrame)
end

# ## Panels A and B: heatmap grids
function heatmap_panels(results, metric; letter="")
    vmin, vmax = extrema(results[!, metric])
    panels = Any[]
    for (row, gname) in enumerate(GENOTYPE_ORDER), (col, cname) in enumerate(CONDITION_ORDER)
        sub = results[(results.genotype.==string(gname)).&(results.condition.==string(cname)), :]
        ## Use the full configured axes rather than whatever happens to be present,
        ## so an incomplete sweep still produces a correctly-shaped grid with the
        ## missing cells drawn as gaps (NaN) instead of throwing.
        fs = sort(FREQS_HZ)
        ds = sort(DURATIONS_S)
        Z = [begin
            v = sub[(sub.freq_Hz.==f).&(sub.dur_s.==d), metric]
            isempty(v) ? NaN : first(v)
        end for d in ds, f in fs]
        is_bottom = row == length(GENOTYPE_ORDER)
        is_left = col == 1
        is_right = col == length(CONDITION_ORDER)
        ttl = GENOTYPE_LABELS[gname] * " — " * CONDITION_LABELS[cname]
        (row == 1 && col == 1 && !isempty(letter)) && (ttl = letter * "   " * ttl)
        push!(panels, heatmap(fs, ds, Z;
            clims=(vmin, vmax), color=:viridis, colorbar=is_right,
            title=ttl, titlefontsize=11,
            xlabel=is_bottom ? "Pacing frequency (Hz)" : "",
            ylabel=is_left ? "Pacing duration (s)" : "",
            guidefontsize=10, tickfontsize=8,
            left_margin=is_left ? 35Plots.mm : 3Plots.mm,
            bottom_margin=is_bottom ? 10Plots.mm : 3Plots.mm,
            top_margin=5Plots.mm,
            right_margin=is_right ? 8Plots.mm : 3Plots.mm))
    end
    return panels
end

# ## Panel C: tau vs pacing duration
# Frequency is encoded by BOTH color and marker shape (so the panel survives
# greyscale printing); genotype by line style and marker fill. Two legend-only
# entries carry the genotype key, keeping the legend to six entries.
function tau_duration_panel(results; condition=:Control, show_freqs=[0.5, 1.0, 2.0, 3.0], letter="C")
    p = plot(xlabel="Pacing duration (s)", ylabel="Post-pacing decay \u03c4 (s)",
        title=letter * "   \u03c4 vs. pacing duration (" * CONDITION_LABELS[condition] * ")",
        titlelocation=:left, titlefontsize=12, guidefontsize=10, tickfontsize=9,
        legend=:outertopright, legendfontsize=9, ylims=(0, 25),
        left_margin=12Plots.mm, bottom_margin=9Plots.mm, top_margin=5Plots.mm)

    vline!(p, [TAU_MIN_DURATION_S], color=:gray40, linestyle=:dot, lw=2,
        lab="\u03c4 stats use \u2265 $(Int(TAU_MIN_DURATION_S)) s")

    for (i, f) in enumerate(show_freqs)
        col = FREQ_COLORS[mod1(i, length(FREQ_COLORS))]
        mk = FREQ_MARKERS[mod1(i, length(FREQ_MARKERS))]
        for (gname, ls, mfill, lab) in ((:WT, :solid, col, "$(f) Hz"),
            (:R275H, :dash, :white, false))
            sub = results[(results.genotype.==string(gname)).&(results.condition.==string(condition)).&(results.freq_Hz.==f), :]
            isempty(sub) && continue
            sort!(sub, :dur_s)
            plot!(p, sub.dur_s, sub.tau; linestyle=ls, color=col, marker=mk, markersize=6,
                markercolor=mfill, markerstrokecolor=col, lab=lab)
        end
    end

    plot!(p, Float64[], Float64[]; color=:black, linestyle=:solid, marker=:circle,
        markersize=6, markercolor=:black, markerstrokecolor=:black, lab="WT (solid, filled)")
    plot!(p, Float64[], Float64[]; color=:black, linestyle=:dash, marker=:circle,
        markersize=6, markercolor=:white, markerstrokecolor=:black, lab="R275H-like (dashed, open)")
    return p
end

# ## Dynamic-range summary
function range_summary_table(results)
    peak_auc = combine(groupby(results, [:genotype, :condition]),
        :peak => (x -> maximum(x) - minimum(x)) => :peak_range,
        :peak => minimum => :peak_min, :peak => maximum => :peak_max,
        :auc => (x -> maximum(x) - minimum(x)) => :auc_range,
        :auc => minimum => :auc_min, :auc => maximum => :auc_max, nrow => :n_peak_auc)
    tau_only = combine(groupby(results[results.dur_s.>=TAU_MIN_DURATION_S, :], [:genotype, :condition]),
        :tau => (x -> maximum(x) - minimum(x)) => :tau_range,
        :tau => minimum => :tau_min, :tau => maximum => :tau_max, nrow => :n_tau)
    return innerjoin(peak_auc, tau_only, on=[:genotype, :condition])
end

# ## Panel D: mutant/WT range ratio, expressed in one common unit (%)
function ratio_panel(range_summary; letter="D")
    getv(col, g, c) = only(range_summary[(range_summary.genotype.==string(g)).&(range_summary.condition.==string(c)), col])
    ratios = Dict(m => [100 * getv(m, :R275H, c) / getv(m, :WT, c) for c in CONDITION_ORDER]
                  for m in (:peak_range, :auc_range, :tau_range))
    p = groupedbar([CONDITION_SHORT[c] for c in CONDITION_ORDER],
        [ratios[:peak_range] ratios[:auc_range] ratios[:tau_range]],
        label=["Peak" "AUC" "\u03c4"], color=[WT_COLOR MUT_COLOR TAU_COLOR],
        ylabel="Mutant range as % of WT", ylims=(0, 140),
        title=letter * "   Dynamic-range ratio (R275H-like / WT)",
        titlelocation=:left, titlefontsize=12,
        guidefontsize=10, tickfontsize=9, legendfontsize=9, legend=:topleft,
        left_margin=12Plots.mm, bottom_margin=9Plots.mm, top_margin=5Plots.mm)
    hline!(p, [100], color=:black, linestyle=:dash, lw=1.5, lab="No change (100%)")
    for m in (:peak_range, :auc_range, :tau_range)
        println("  $m (mutant/WT %): ", round.(ratios[m], digits=1))
    end
    return p
end

# ## Main
results = run_sweep()
println("Rows for plotting: ", nrow(results))

# ## Completeness check
# An interrupted or partially failed sweep yields a grid with holes. Those cells
# are drawn as gaps rather than crashing the plot, but any summary statistic
# computed from an incomplete grid (the dynamic ranges below) would be wrong, so
# the gaps are reported loudly here. Simply re-running this script fills them:
# completed cells are cached and skipped.
let expected = length(GENOTYPES) * length(CONDITIONS) * length(FREQS_HZ) * length(DURATIONS_S)
    missing_cells = NamedTuple[]
    for (gname, _) in pairs(GENOTYPES), (cname, _) in pairs(CONDITIONS)
        sub = results[(results.genotype.==string(gname)).&(results.condition.==string(cname)), :]
        for d in DURATIONS_S, f in FREQS_HZ
            isempty(sub[(sub.freq_Hz.==f).&(sub.dur_s.==d), :]) &&
                push!(missing_cells, (; genotype=gname, condition=cname, freq_Hz=f, dur_s=d))
        end
    end
    if isempty(missing_cells)
        println("Grid complete: $(nrow(results))/$expected cells.")
    else
        println("\n" * "="^70)
        println("WARNING: grid INCOMPLETE — $(length(missing_cells)) of $expected cells missing.")
        println("Range statistics below are computed from an incomplete grid and")
        println("should NOT be quoted. Re-run this script to fill the gaps;")
        println("cached cells are skipped, so only the missing ones are solved.")
        println("Missing:")
        for m in missing_cells
            println("  $(m.genotype)/$(m.condition)  f=$(m.freq_Hz) Hz  d=$(m.dur_s) s")
        end
        println("="^70 * "\n")
    end
end

range_summary = range_summary_table(results)
println("\nDynamic-range summary (tau restricted to durations >= $(TAU_MIN_DURATION_S) s):")
println(range_summary)
CSV.write(OUT * "-dynamic-range.csv", range_summary)
println("\nRatios:")

panelC = tau_duration_panel(results)
panelD = ratio_panel(range_summary)

# Individual panel files (recommended for assembling the submitted figure; the
# composite squashes the heatmap rows vertically). No in-panel letter here, since
# each file carries its own title.
figA = plot(heatmap_panels(results, :peak)..., layout=(2, 4), size=(2000, 900),
    plot_title="A   Peak active CaMKII fraction", plot_titlefontsize=15)
savefig(figA, OUT * "-A-peak.png"); savefig(figA, OUT * "-A-peak.pdf")

figB = plot(heatmap_panels(results, :auc)..., layout=(2, 4), size=(2000, 900),
    plot_title="B   Area under the active-CaMKII curve (AU\u00b7s)", plot_titlefontsize=15)
savefig(figB, OUT * "-B-auc.png"); savefig(figB, OUT * "-B-auc.pdf")

savefig(plot(panelC, size=(1100, 560)), OUT * "-C-tau-vs-duration.png")
savefig(plot(panelC, size=(1100, 560)), OUT * "-C-tau-vs-duration.pdf")
savefig(plot(panelD, size=(1000, 560)), OUT * "-D-range-ratio.png")
savefig(plot(panelD, size=(1000, 560)), OUT * "-D-range-ratio.pdf")

# Composite (letters are in-panel here, since there are no per-section titles)
composite = plot(vcat(heatmap_panels(results, :peak; letter="A"),
        heatmap_panels(results, :auc; letter="B"), [panelC, panelD])...,
    layout=@layout([grid(4, 4); [c{0.5w} d{0.5w}]]),
    size=(2100, 2000), left_margin=5Plots.mm, right_margin=5Plots.mm, bottom_margin=5Plots.mm)
savefig(composite, OUT * "-composite.png"); savefig(composite, OUT * "-composite.pdf")

# Supplementary tau heatmap (not part of Figure 11)
figS = plot(heatmap_panels(results, :tau)..., layout=(2, 4), size=(2000, 900),
    plot_title="Supplementary: post-pacing decay \u03c4 (s)", plot_titlefontsize=15)
savefig(figS, OUT * "-S-tau-heatmap.png"); savefig(figS, OUT * "-S-tau-heatmap.pdf")

println("\nWrote figure11-{A-peak, B-auc, C-tau-vs-duration, D-range-ratio, composite, S-tau-heatmap}.png/.pdf")
println("Wrote figure11-dynamic-range.csv")
