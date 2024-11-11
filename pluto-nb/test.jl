### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ f4d83570-9d1c-11ef-0654-53e2d6f12d6e
begin
    import Pkg
    # careful: this is _not_ a reproducible environment
    # activate the global environment
    Pkg.activate(Base.current_project())
    using ProgressLogging
	using ModelingToolkit
	using OrdinaryDiffEq
	using DiffEqCallbacks
	using Plots
	using LsqFit
	using CaMKIIModel
	using CaMKIIModel: μM, Hz
	Plots.default(lw=1.5, fmt=:png)
end

# ╔═╡ fe7e7c0b-1ae0-4a04-9549-23aba60efcf5
Base.current_project()

# ╔═╡ 366ecfd5-358c-4bf6-8e8d-c494c0a3f557
sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true);

# ╔═╡ fb697567-fa32-4963-a490-0668657b22f5
tend = 100.0

# ╔═╡ 638be248-21ad-4833-bbc8-f2f9a805a900
ps = [
	sys.ktrCaSR => 50Hz,
	sys.kRyR => 20Hz,
	sys.GCaL => 1.3e-4, # 6.3e-5
]

# ╔═╡ 78fb22d4-2f06-49b0-a818-377eedf7ca8d
prob = ODEProblem(sys, [], tend, ps)

# ╔═╡ 0d97c77a-fdd9-463c-9c53-d40be74b266b
@unpack Istim = sys;

# ╔═╡ 12ded313-721a-4d7b-a2e3-25b313c1a582
callback = build_stim_callbacks(Istim, tend; period=1)

# ╔═╡ 761a0baf-b5ec-4caa-ae56-5808a2840987
sol = solve(prob, FBDF(); callback, abstol=1e-6, reltol=1e-6, progress=true)

# ╔═╡ 8a2eb462-92ee-4ed3-a315-3c58a9a194f4
plot(sol, idxs=sys.vm*1000, lab="Membrane potential", tspan=(90, 92))

# ╔═╡ 152536f6-3cb9-48dd-9467-f1a2b64b3dba
# Calcium transient is smaller than paper (550 nM vs 800 nM)
plot(sol, idxs=[sys.Cai_sub_SL*1E6, sys.Cai_sub_SR*1E6], tspan=(90, 92), ylabel="nM", lab=["Ca (sub SL)" "Ca (sub SR)"], ylims=(100, 800))

# ╔═╡ cd11a3d5-bdd8-415c-8b8d-58c4228f8877
# SR Ca release is reletively small (see JSR)
plot(sol, idxs=[sys.CaJSR, sys.CaNSR], tspan=(90, 92), ylabel="mM")

# ╔═╡ 58073dd2-9670-48d1-876e-82ba2a1f5611
begin
	# ICaL and INaCa are OK
	# RyR flux is tiny in NRVM
	@unpack INaCa, ICaL = sys
	plot(sol, idxs=[INaCa, ICaL], tspan=(90, 92), ylabel="uA/uF or uM/ms")
end

# ╔═╡ b563ef4c-1dc3-4b3e-9f58-bace24d03d9e
plot(sol, idxs=sys.CaMKAct, ylim=(0.18, 0.28))

# ╔═╡ Cell order:
# ╠═fe7e7c0b-1ae0-4a04-9549-23aba60efcf5
# ╠═f4d83570-9d1c-11ef-0654-53e2d6f12d6e
# ╠═366ecfd5-358c-4bf6-8e8d-c494c0a3f557
# ╠═fb697567-fa32-4963-a490-0668657b22f5
# ╠═638be248-21ad-4833-bbc8-f2f9a805a900
# ╠═78fb22d4-2f06-49b0-a818-377eedf7ca8d
# ╠═0d97c77a-fdd9-463c-9c53-d40be74b266b
# ╠═12ded313-721a-4d7b-a2e3-25b313c1a582
# ╠═761a0baf-b5ec-4caa-ae56-5808a2840987
# ╠═8a2eb462-92ee-4ed3-a315-3c58a9a194f4
# ╠═152536f6-3cb9-48dd-9467-f1a2b64b3dba
# ╠═cd11a3d5-bdd8-415c-8b8d-58c4228f8877
# ╠═58073dd2-9670-48d1-876e-82ba2a1f5611
# ╠═b563ef4c-1dc3-4b3e-9f58-bace24d03d9e