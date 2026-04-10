# CaMKII effects

## LCC

```julia
k_ckLCC = 0.4                   # [s^-1]
k_pp1LCC = 0.1103               # [s^-1]
k_pp2aLCC = 10.1                # [s^-1]
KmCK_LCC = 12                   # [uM]
KmPKA_LCC = 21                  # [uM]
KmPP2A_LCC = 47                 # [uM]
KmPP1_LCC = 9                   # [uM]
LCCtotSL = 0.0846          # [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
PP1_SL = 0.57              # [uM] - Total Subsarcolemmal [PP1]
CaMKIItotSL = 120 * 8.293e-4      # [uM]
CaMKIIact_SL = CaMKIItotSL .* (Pb_sl + Pt_sl + Pt2_sl + Pa_sl)
D(LCC_CKslp) ~ (LCCSL_PHOS - LCCSL_DEPHOS)
LCC_CKsln = LCCtotSL - LCC_CKslp
LCCSL_PHOS = (k_ckLCC * CaMKIIact_SL * LCC_CKsln) / (KmCK_LCC + LCC_CKsln)
LCCSL_DEPHOS = (k_pp1LCC * PP1_SL * LCC_CKslp) / (KmPP1_LCC + LCC_CKslp)
```

## RyR

Does NVRM have dyadic structures?

Increasing SR leak.

```julia
k_ckRyR = 0.4                   # [s^-1]
k_pp1RyR = 1.07                 # [s^-1]
k_pp2aRyR = 0.481               # [s^-1]
kb_2815 = 0.35                  # [uM/s] - CaMKII site
KmCK_RyR = 12                   # [uM]
KmPP1_RyR = 9                   # [uM]
KmPP2A_RyR = 47                 # [uM]
RyRtot = 382.6             # [uM] - Total RyR (in Dyad)
PP1_dyad = 95.7            # [uM] - Total dyadic [PP1]
PP2A_dyad = 95.76          # [uM] - Total dyadic PP2A
CaMKIItotDyad = 120             # [uM]
PP1_PLB_avail = 1 - I1p_PP1 / PP1_PLBtot + 0.081698
PP1_PLB = PP1_dyad * PP1_PLB_avail
D(RyR2815p) ~ (RyR_BASAL - RyR_PP1_DEPHOS - RyR_PP2A_DEPHOS)
RyR_BASAL = kb_2815 * RyR2815n
RyR_PHOS = (k_ckRyR*CaMKIIact_Dyad*RyR2815n)/(KmCK_RyR+RyR2815n)
RyR_PP1_DEPHOS = (k_pp1RyR * PP1_dyad * RyR2815p) / (KmPP1_RyR + RyR2815p)
RyR_PP2A_DEPHOS = (k_pp2aRyR * PP2A_dyad * RyR2815p) / (KmPP2A_RyR + RyR2815p)
RyR_CKp = RyR2815p / RyRtot
kleak = (1 / 2 + 5 * RyR_CKp / 2) * 5e-6Hz
```

## PLB

Decreasing Km for SERCA

```julia
k_ckPLB = 8e-3                  # [s^-1]
k_pp1PLB = 0.0428               # [s^-1]
PP1_PLBtot = 0.89               # [uM] - [umol/L cytosol]
KmCK_PLB = 12 # [uM]
KmPP1_PLB = 9 # [uM]
PLB_PHOS = (k_ckPLB*PLBT17n*CaMKIIact_Dyad)/(KmCK_PLB+PLBT17n)
PLB_DEPHOS = (k_pp1PLB * PP1_PLB * PLBT17p) / (KmPP1_PLB + PLBT17p)
PLB_CKp = PLBT17p / PLBtot
fCKII_PLB = (1 - 0.5 * PLB_CKp)  # Max effect: fCKII_PLB=0.5
```
