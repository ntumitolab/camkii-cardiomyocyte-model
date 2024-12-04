### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 34de6de0-b14c-11ef-3c48-a3452d7d7f85
md"""
# Common units

For electrophysiology
"""

# ╔═╡ 600921d1-9338-48a5-8c11-486cf14836f4
# Microseconds
ms = 1

# ╔═╡ 69259624-2f29-4828-9fb1-f8bf5e9997af
second = 1000ms

# ╔═╡ 82854a56-9131-4525-9986-cd0cac3869b5
Hz = 1/second

# ╔═╡ a77a1d2f-968f-4557-a132-a2489d78feb5
mmol = 1

# ╔═╡ 4cdf9934-bc30-4456-a06b-76ed19a2ba0f
mol = 1000mmol

# ╔═╡ 2f51d1ba-b241-4eb4-9ca2-f67ea1f26bc8
metre = 1

# ╔═╡ 08033864-3fc2-4f70-b3ff-2dc110b39885
μM = mmol/metre^3

# ╔═╡ 8a6e71dc-ad8f-4246-947c-26a98ce79a8b
μM/ms

# ╔═╡ 25c59861-7469-4f5f-838d-c2f137841d34
mV = 1

# ╔═╡ 7e46345b-8af8-4d9e-9951-e0132f0995e9
# J / mol·K
R = 8.314

# ╔═╡ 197e0bf3-3856-422a-87da-c4a0b1d36dfa
Columb = 1

# ╔═╡ 627b72f4-4fbb-41bc-a32f-aafed6d84563
μA = Columb / second / 10^6

# ╔═╡ 9b01cb43-6da7-42ac-afe5-ecce0c986003
Faraday = 96485 * Columb / mol

# ╔═╡ 2fda1066-de1e-435c-bd2d-a4ff66818bcf
# Kelvin
Temp = 310

# ╔═╡ 5f2f7bcc-c2d4-4b5a-b021-dc9504ace92a
# Thermal voltage (mV)
vt = R * Temp / Faraday

# ╔═╡ 741f88c4-b93b-42eb-bbd4-e23ca9edb5c4
# Current density
μA_μF = mV / ms

# ╔═╡ 90544cd6-d305-4d80-b8c8-4efe29b35896
# Conductance
mS_μF = μA_μF / mV

# ╔═╡ 27faf21d-fab3-4259-9677-7da769d3e2c7
mm = 1e-3metre

# ╔═╡ 5b45cfdd-2a49-45ad-b1a5-2698c573f4bb
μm = 1e-6metre

# ╔═╡ 40dc4501-5771-4977-9922-eadcffd6e7e8
litre = 1e-3metre^3

# ╔═╡ c758f8ef-c68a-49c6-8d21-993b90b079d6
# picoliter
pL = 1e-12litre

# ╔═╡ 3295ec28-e23c-4f7c-aa9f-e5868824f121


# ╔═╡ Cell order:
# ╠═34de6de0-b14c-11ef-3c48-a3452d7d7f85
# ╠═600921d1-9338-48a5-8c11-486cf14836f4
# ╠═69259624-2f29-4828-9fb1-f8bf5e9997af
# ╠═82854a56-9131-4525-9986-cd0cac3869b5
# ╠═a77a1d2f-968f-4557-a132-a2489d78feb5
# ╠═4cdf9934-bc30-4456-a06b-76ed19a2ba0f
# ╠═2f51d1ba-b241-4eb4-9ca2-f67ea1f26bc8
# ╠═08033864-3fc2-4f70-b3ff-2dc110b39885
# ╠═8a6e71dc-ad8f-4246-947c-26a98ce79a8b
# ╠═25c59861-7469-4f5f-838d-c2f137841d34
# ╠═7e46345b-8af8-4d9e-9951-e0132f0995e9
# ╠═197e0bf3-3856-422a-87da-c4a0b1d36dfa
# ╠═627b72f4-4fbb-41bc-a32f-aafed6d84563
# ╠═9b01cb43-6da7-42ac-afe5-ecce0c986003
# ╠═2fda1066-de1e-435c-bd2d-a4ff66818bcf
# ╠═5f2f7bcc-c2d4-4b5a-b021-dc9504ace92a
# ╠═741f88c4-b93b-42eb-bbd4-e23ca9edb5c4
# ╠═90544cd6-d305-4d80-b8c8-4efe29b35896
# ╠═27faf21d-fab3-4259-9677-7da769d3e2c7
# ╠═5b45cfdd-2a49-45ad-b1a5-2698c573f4bb
# ╠═40dc4501-5771-4977-9922-eadcffd6e7e8
# ╠═c758f8ef-c68a-49c6-8d21-993b90b079d6
# ╠═3295ec28-e23c-4f7c-aa9f-e5868824f121
