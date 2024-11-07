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
	using CaMKIIModel: μM, Hz, μm
	Plots.default(lw=1.5)
end

# ╔═╡ fe7e7c0b-1ae0-4a04-9549-23aba60efcf5
Base.current_project()

# ╔═╡ 366ecfd5-358c-4bf6-8e8d-c494c0a3f557
sys = build_neonatal_ecc_sys(simplify=true, reduce_iso=true);

# ╔═╡ fb697567-fa32-4963-a490-0668657b22f5
tend = 300.0

# ╔═╡ 78fb22d4-2f06-49b0-a818-377eedf7ca8d
prob = ODEProblem(sys, [], tend, [sys.ktrCaSR => 5])

# ╔═╡ 0d97c77a-fdd9-463c-9c53-d40be74b266b
@unpack Istim = sys;

# ╔═╡ 12ded313-721a-4d7b-a2e3-25b313c1a582
callback = build_stim_callbacks(Istim, tend; period=1)

# ╔═╡ 761a0baf-b5ec-4caa-ae56-5808a2840987
sol = solve(prob, FBDF(); callback, abstol=1e-6, reltol=1e-6, progress=true)

# ╔═╡ 152536f6-3cb9-48dd-9467-f1a2b64b3dba
# Calcium transient is smaller than paper (550 nM vs 800 nM)
plot(sol, idxs=[sys.Cai_sub_SL*1E6, sys.Cai_sub_SR*1E6], tspan=(100, 101), ylabel="nM", lab=["Ca (sub SL)" "Ca (sub SR)"])

# ╔═╡ 58073dd2-9670-48d1-876e-82ba2a1f5611
begin
	# ICaL and INaCa are OK
	# RyR flux is tiny in NRVM
	@unpack INaCa, ICaL, ICaT, ICab, JCa_SL, JCa_SR = sys
	plot(sol, idxs=[INaCa, ICaL, JCa_SL, JCa_SR], tspan=(100, 101), ylabel="uA/uF or uM/ms")
end

# ╔═╡ 0920584e-d6c5-494e-9c92-178fdc3815f5
begin
	@unpack Jup, Jrel, Jtr = sys
	plot(sol, idxs=[Jup, Jtr], tspan=(100, 102), ylabel="uM/ms")
end

# ╔═╡ cd11a3d5-bdd8-415c-8b8d-58c4228f8877
# SR Ca release is reletively small (see JSR)
plot(sol, idxs=[sys.CaJSR, sys.CaNSR], tspan=(100, 102), ylabel="mM")

# ╔═╡ 89bb2f80-8f4c-44fe-90dc-d22e9314e5aa
plot(sol, idxs=[sys.betaSR], tspan=(100, 102))

# ╔═╡ Cell order:
# ╠═fe7e7c0b-1ae0-4a04-9549-23aba60efcf5
# ╠═f4d83570-9d1c-11ef-0654-53e2d6f12d6e
# ╠═366ecfd5-358c-4bf6-8e8d-c494c0a3f557
# ╠═fb697567-fa32-4963-a490-0668657b22f5
# ╠═78fb22d4-2f06-49b0-a818-377eedf7ca8d
# ╠═0d97c77a-fdd9-463c-9c53-d40be74b266b
# ╠═12ded313-721a-4d7b-a2e3-25b313c1a582
# ╠═761a0baf-b5ec-4caa-ae56-5808a2840987
# ╠═152536f6-3cb9-48dd-9467-f1a2b64b3dba
# ╠═58073dd2-9670-48d1-876e-82ba2a1f5611
# ╠═0920584e-d6c5-494e-9c92-178fdc3815f5
# ╠═cd11a3d5-bdd8-415c-8b8d-58c4228f8877
# ╠═89bb2f80-8f4c-44fe-90dc-d22e9314e5aa
