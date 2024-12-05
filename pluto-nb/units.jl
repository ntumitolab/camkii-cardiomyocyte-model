### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 34de6de0-b14c-11ef-3c48-a3452d7d7f85
md"""
# Common units

For electrophysiology models

- time as millisecond (ms)
- concentration as micromolar (μM)
- voltage as millivolt (mV)
- amount of substance as millimolar (mmol)
- current density as μA/μF

Thus, 

- Faraday constant = 96.485 columb / mmol
- Length unit = **meter (m)** since μM = mmol * $m^{-3}$
- Current density = μA/μF = mV/ms
"""

# ╔═╡ 69259624-2f29-4828-9fb1-f8bf5e9997af
second = 10^3

# ╔═╡ 4cdf9934-bc30-4456-a06b-76ed19a2ba0f
mol = 1000

# ╔═╡ c74d07d7-ae92-4970-8cf5-6f49155e8e07
metre = 1

# ╔═╡ e70930d7-b131-40d3-9d3e-c5448ef96cd3
Ampere = 1

# ╔═╡ 93b7ce14-dd13-4a54-8cac-df854406013b
Joule = 10^6

# ╔═╡ e267f81e-6878-4dd4-8d43-d096a7b7293e
Kelvin = 1

# ╔═╡ 600921d1-9338-48a5-8c11-486cf14836f4
ms = second/10^3

# ╔═╡ 82854a56-9131-4525-9986-cd0cac3869b5
Hz = inv(second)

# ╔═╡ a77a1d2f-968f-4557-a132-a2489d78feb5
mmol = mol/10^3

# ╔═╡ f29a4041-f2f3-4bbb-bb90-eb3ff5ab6124
μM = mmol/metre^3

# ╔═╡ c64a8b61-3f23-45c9-8e4a-77f603c3d0e6
@show μM/ms

# ╔═╡ 197e0bf3-3856-422a-87da-c4a0b1d36dfa
Columb = Ampere * second

# ╔═╡ e9b68a96-cf94-4dee-b1d6-0ceae87209ab
Volt = Joule / Columb

# ╔═╡ 25c59861-7469-4f5f-838d-c2f137841d34
mV = Volt/10^3

# ╔═╡ 6ed5c1ec-f7e5-4111-8a80-2ce567d9a7a5
@show mV/ms

# ╔═╡ 7e46345b-8af8-4d9e-9951-e0132f0995e9
# J / mol·K
R = 8.314Joule/Kelvin/mol

# ╔═╡ 627b72f4-4fbb-41bc-a32f-aafed6d84563
μA = Ampere / 10^6

# ╔═╡ 3c566405-7cd7-4e44-8b27-3b833352d495
Farad = Columb / Volt

# ╔═╡ 1d5c936e-b424-4377-b41d-f4059acadcb9
μF = Farad/10^6

# ╔═╡ a6b5e066-e737-45bb-af18-e40b76e7b7fc
# Current density
@show μA/μF

# ╔═╡ 9b01cb43-6da7-42ac-afe5-ecce0c986003
Faraday = 96485 * Columb / mol

# ╔═╡ 2fda1066-de1e-435c-bd2d-a4ff66818bcf
Temp = 310Kelvin

# ╔═╡ 5f2f7bcc-c2d4-4b5a-b021-dc9504ace92a
# Thermal voltage (mV)
vt = R * Temp / Faraday

# ╔═╡ 4e12a340-4700-428f-b23e-5d954b4a6994
Siemens = Ampere/Volt

# ╔═╡ 90544cd6-d305-4d80-b8c8-4efe29b35896
# Conductance
mS_μF = Siemens / 10^3 / (Farad / 10^6)

# ╔═╡ 27faf21d-fab3-4259-9677-7da769d3e2c7
mm = 1e-3metre

# ╔═╡ fa990df3-d2c1-4797-9416-b755de09b51b
cm = 0.01metre

# ╔═╡ 5b45cfdd-2a49-45ad-b1a5-2698c573f4bb
μm = 1e-6metre

# ╔═╡ 40dc4501-5771-4977-9922-eadcffd6e7e8
litre = 1000cm^3

# ╔═╡ c758f8ef-c68a-49c6-8d21-993b90b079d6
# picoliter
pL = 1e-12litre

# ╔═╡ 3295ec28-e23c-4f7c-aa9f-e5868824f121
cm² = cm^2

# ╔═╡ 5a156633-4635-4b4f-8456-e0764db626dc
@show μA/cm² μF/cm² 

# ╔═╡ bdd32a87-a075-4b81-803a-6787e8834678
md"""
Conversion factor between current and flux

$\frac{A_{cap} * C_m}{z * Vol * F}$

Volume should be mm^3 (or μL) to match SI unit?
"""

# ╔═╡ ebee4c1f-4b0a-4b7c-a30b-e3284b95b965
(cm² * μF/cm²) / (10^(-9) * metre^3 * Columb / mmol)

# ╔═╡ c5972a01-aa29-41de-9309-8568f0eef314
rSL_true = 10.5μm

# ╔═╡ c240c987-50e2-41a7-8a62-7e6b135d41cc
rSR_true = 6μm

# ╔═╡ dccd732a-8fcf-4711-ab88-e69c873c198d
Acap = 4π * rSL_true^2

# ╔═╡ 6c054eb3-91a7-4553-96ac-3893420014d8
Cm = 1μF / cm²

# ╔═╡ 79f70c73-81f6-42e6-ae05-1d99dda07962
Vmyo = 4 / 3 * π * (rSL_true^3 - rSR_true^3)

# ╔═╡ 0a58a619-cfae-41e6-81fa-7a4f627dd26a
Acap * Cm / Vmyo / Faraday

# ╔═╡ Cell order:
# ╠═34de6de0-b14c-11ef-3c48-a3452d7d7f85
# ╠═69259624-2f29-4828-9fb1-f8bf5e9997af
# ╠═4cdf9934-bc30-4456-a06b-76ed19a2ba0f
# ╠═c74d07d7-ae92-4970-8cf5-6f49155e8e07
# ╠═e70930d7-b131-40d3-9d3e-c5448ef96cd3
# ╠═93b7ce14-dd13-4a54-8cac-df854406013b
# ╠═e267f81e-6878-4dd4-8d43-d096a7b7293e
# ╠═600921d1-9338-48a5-8c11-486cf14836f4
# ╠═82854a56-9131-4525-9986-cd0cac3869b5
# ╠═a77a1d2f-968f-4557-a132-a2489d78feb5
# ╠═f29a4041-f2f3-4bbb-bb90-eb3ff5ab6124
# ╠═c64a8b61-3f23-45c9-8e4a-77f603c3d0e6
# ╠═197e0bf3-3856-422a-87da-c4a0b1d36dfa
# ╠═e9b68a96-cf94-4dee-b1d6-0ceae87209ab
# ╠═25c59861-7469-4f5f-838d-c2f137841d34
# ╠═6ed5c1ec-f7e5-4111-8a80-2ce567d9a7a5
# ╠═7e46345b-8af8-4d9e-9951-e0132f0995e9
# ╠═627b72f4-4fbb-41bc-a32f-aafed6d84563
# ╠═3c566405-7cd7-4e44-8b27-3b833352d495
# ╠═1d5c936e-b424-4377-b41d-f4059acadcb9
# ╠═a6b5e066-e737-45bb-af18-e40b76e7b7fc
# ╠═9b01cb43-6da7-42ac-afe5-ecce0c986003
# ╠═2fda1066-de1e-435c-bd2d-a4ff66818bcf
# ╠═5f2f7bcc-c2d4-4b5a-b021-dc9504ace92a
# ╠═4e12a340-4700-428f-b23e-5d954b4a6994
# ╠═90544cd6-d305-4d80-b8c8-4efe29b35896
# ╠═27faf21d-fab3-4259-9677-7da769d3e2c7
# ╠═fa990df3-d2c1-4797-9416-b755de09b51b
# ╠═5b45cfdd-2a49-45ad-b1a5-2698c573f4bb
# ╠═40dc4501-5771-4977-9922-eadcffd6e7e8
# ╠═c758f8ef-c68a-49c6-8d21-993b90b079d6
# ╠═3295ec28-e23c-4f7c-aa9f-e5868824f121
# ╠═5a156633-4635-4b4f-8456-e0764db626dc
# ╠═bdd32a87-a075-4b81-803a-6787e8834678
# ╠═ebee4c1f-4b0a-4b7c-a30b-e3284b95b965
# ╠═c5972a01-aa29-41de-9309-8568f0eef314
# ╠═c240c987-50e2-41a7-8a62-7e6b135d41cc
# ╠═dccd732a-8fcf-4711-ab88-e69c873c198d
# ╠═6c054eb3-91a7-4553-96ac-3893420014d8
# ╠═79f70c73-81f6-42e6-ae05-1d99dda07962
# ╠═0a58a619-cfae-41e6-81fa-7a4f627dd26a
