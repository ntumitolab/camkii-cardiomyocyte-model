# Parameters

## General parameters

| Parameter | Value    | Units  | Description                      |
| --------- | -------- | ------ | -------------------------------- |
| r_SR      | 6        | μm     | Radius of SR                     |
| r_SL      | 10.5     | μm     | Radius of sarcolemma             |
| V_SR      | 0.0903   | pL     | SR volume                        |
| V_NSR     | 0.9V_SR  | pL     | Network SR volume                |
| V_JSR     | 0.1V_JSR | pL     | Junctional SR volume             |
| V_subSR   | 0.046    | pL     | Sub-SR volume                    |
| V_subSL   | 0.137    | pL     | Sub-sarcolemma volume            |
| V_myo     | 3.944    | pL     | Cytosolic volume                 |
| A_cap     | 1385.44  | μm²    | Cell membrane area               |
| C_m       | 1        | μFcm⁻² | Cell membrane capacitance        |
| ca_o      | 1796     | μM     | External calcium concentration   |
| na_o      | 154578   | μM     | External sodium concentration    |
| k_o       | 5366     | μM     | External potassium concentration |
| `[ATP]`   | 5000     | μM     | ATP concentration                |
| V_T       | 26.7     | mV     | Thermal voltage                  |

## Calcium diffusion

| Parameter  | Value    | Units   | Description                        |
| ---------- | -------- | ------- | ---------------------------------- |
| ΣTrpn      | 35       | μM      | Total troponin content             |
| Km_Trpn    | 0.5      | μM      | Half-saturation Ca concentration   |
| ΣCmdn      | 30       | μM      | Total calmodulin content           |
| Km_Cmdn    | 2.38     | μM      | Half-saturation Ca concentration   |
| D_ca       | 7        | μm²ms⁻¹ | Calcium diffusion rate             |
| dx         | 0.1      | μm      | Discretization distance            |
| fracTnIp_0 | 0.062698 | -       | Baseline effect of PKA on Troponin |

## Endoplasmic reticulum

| Parameter    | Value          | Units | Description                            |
| ------------ | -------------- | ----- | -------------------------------------- |
| k_RyR        | 20             | 1/s   | RyR permeability                       |
| kapos_RyR    | 1000           | 1/s   | RyR state transition rate              |
| kaneg_RyR    | 160            | 1/s   | RyR state transition rate              |
| Vmax_SR      | 999.6          | μM/s  | SERCA reaction rate                    |
| Kmf_SR       | 0.5            | μM    | Calcium affinity for SERCA             |
| Kmr_SR       | 7000$Kmf_{SR}$ | μM    | Calcium affinity for SERCA             |
| kSR_leak     | 0.005          | 1/s   | SR leak rate                           |
| ktrCa_SR     | 50             | 1/s   | Calcium diffusion rate from NSR to JSR |
| ΣCSQN$       | 24750          | μM    | Calsequestrin concentration            |
| Km_csqn      | 800            | μM    | Calcium affinity for calsequestrin     |
| fracPKA_PLBo | 00.920245      | -     |                                        |

## Sarcolemmal sodium channels

| Parameter | Value  | Units | Description                            |
| --------- | ------ | ----- | -------------------------------------- |
| G_Na      | 12.8   | mS/μF | Fast sodium channels conductance       |
| G_Nab     | 0.0026 | mS/μF | Background sodium channels conductance |

## Sarcolemmal potassium channels

| Parameter | Value  | Units | Description                                      |
| --------- | ------ | ----- | ------------------------------------------------ |
| G_K1      | 0.0515 | mS/μF | Potassium channels conductance                   |
| G_t       | 0.1    | mS/μF | Transient outward potassium channels conductance |
| G_Ks      | 0.05   | mS/μF | Potassium channels conductance                   |
| τ_nKs     | 750    | ms    | Potassium channels time scale                    |
| G_Kr      | 0.06   | mS/μF | Potassium channels conductance                   |
| k_fIKr    | 23.761 | 1/s   | Potassium channels transition rate               |
| k_bIKr    | 36.778 | 1/s   | Potassium channels transition rate               |
| G_f       | 0.021  | mS/μF | Funny current conductance                        |
| f_Na      | 0.021  | -     | Funny current sodium fraction                    |
| f_is      | 0.706  | -     | Transient outward gating variable                |

## Sarcolemmal calcium channels

| Parameter    | Value         | Units      | Description |
| ------------ | ------------- | ---------- | ----------- |
| fNaCa        | 1             | -          |             |
| kNaCa        | 2.268 * 10⁻¹⁶ | μAμF⁻¹μM⁻⁴ |             |
| dNaCa        | 10⁻¹⁶         | μM⁻⁴       |             |
| gNaCa        | 0.5           | -          |             |
| G_CaL        | 6.3 * 10⁻⁵    | m³s⁻¹F⁻¹   |             |
| τ_fca        | 10            | ms         |             |
| g_CaT        | 0.2           | mSμF⁻¹     |             |
| g_Cab        | 0.0008        | mSμF⁻¹     |             |
| ICascale_0   | 0.95          | -          |             |
| fracLCCbp_0  | 0.250657      | -          |             |
| fracLCCbpISO | 0.525870      | -          |             |

## Sarcolemmal Na-K pump

| Parameter | Value | Units | Description                              |
| --------- | ----- | ----- | ---------------------------------------- |
| Imax_NaK  | 2.7   | μA/μF | Maximal rate of Na-K pump                |
| KmNai_NaK | 18600 | μM    |                                          |
| KmKo_NaK  | 1500  | μM    |                                          |
| nNaK      | 3.2   | -     | Hill coefficient for sodium of Na-K pump |

## B1AR

| Parameter  | Value     | Units | Description                            |
| ---------- | --------- | ----- | -------------------------------------- |
| PKACI_0    | 0.0734    | -     | Basal PKACI activity                   |
| PKACI_act  | 0.1995    | -     | Activated PKACI activity               |
| PKACI_KM   | 0.0139    | μM    | PKACI sensitivity to ISO               |
| PKACII_0   | 0.1840    | -     | Basal PKACII activity                  |
| PKACII_act | 0.3444    | -     | Activated PKACII activity              |
| PKACII_KM  | 0.0103    | μM    | PKACII sensitivity to ISO              |
| PP1_0      | 0.8927    | -     | Basal PP1 activity                     |
| PP1_act    | 0.0492    | -     | Activated PP1 activity                 |
| PP1_KI     | 0.00637   | μM    | PP1 sensitivity to ISO                 |
| PLBp_0     | 0.0824    | -     | Basal PLB phosphorylation              |
| PLBp_act   | 0.7961    | -     | Activated PLB phosphorylation          |
| PLBp_KM    | 0.00597   | μM    | PLB phosphorylation sensitivity to ISO |
| PLBp_n     | 1.8167    | -     | Hill coefficient for ISO               |
| PLMp_0     | 0.1172    | -     | Basal PLMp phosphorylation             |
| PLMp_act   | 0.6645    | -     | Activated PLMp phosphorylation         |
| PLMp_KM    | 0.00823   | μM    | PLM phosphorylation sensitivity to ISO |
| PLMp_n     | 1.35784   | -     | Hill coefficient for ISO               |
| TnIp_0     | 0.0669    | -     |                                        |
| TnIp_act   | 0.7524    | -     |                                        |
| TnIp_KM    | 0.007913  | μM    |                                        |
| TnIp_n     | 1.6736    | -     |                                        |
| LCCap_0    | 0.2205    | -     |                                        |
| LCCap_act  | 0.2339    | -     |                                        |
| LCCap_KM   | 0.00726   | μM    |                                        |
| LCCbp_0$   | 0.2517    | -     |                                        |
| LCCbp_act  | 0.2461    | -     |                                        |
| LCCbp_KM   | 0.00695   | μM    |                                        |
| KURp_0     | 0.4390    | -     |                                        |
| KURp_act   | 0.2563    | -     |                                        |
| KURp_KM    | 0.00557   | μM    |                                        |
| RyRp_0     | 0.2054    | -     |                                        |
| RyRp_act   | 0.2399    | -     |                                        |
| RyRp_KM    | 0.0075135 | μM    |                                        |

## CaMKII

| Parameter  | Value  | Units | Description                                               |
| ---------- | ------ | ----- | --------------------------------------------------------- |
| v_CaMK     | 3      | Hz    | CaMK-CaM binding rate                                     |
| r_CaMKO    | 0      | -     | Oxidized CaMK-CaM binding ratio                           |
| r_CaMKP    | 0      | -     | Phosphorylated CaMK-CaM binding ratio                     |
| kb_CaMKP   | 1/3    | Hz    | Dissociation rate of CaMKP                                |
| kfa2_CaMK  | 0.2650 | -     | Maximal CaM-Ca2 binding ratio                             |
| kfa4_CaMK  | 0.1636 | -     | Maximal CaM-Ca4 binding ratio                             |
| kfb_CaMK   | 0.001  | -     | Basal CaMK-CaM binding ratio                              |
| kmCa2_CaMK | 0.7384 | μM    | Half-saturation calcium concentration for CaM-Ca2 binding |
| kmCa4_CaMK | 1.2513 | μM    | Half-saturation calcium concentration for CaM-Ca4 binding |
| kphos_CaMK | 5      | Hz    | Autophosphorylation rate                                  |
| kdeph_CaMK | 1/6    | Hz    | Dephosphorylation rate                                    |
| k_P1P2     | 1/60   | Hz    | Second autophosphorylation rate                           |
| k_P2P1     | 1/15   | Hz    | Second dephosphorylation rate                             |
| kox_CaMK   | 291    | Hz/mM | Oxidation rate                                            |
| krd_CaMK   | 1/45   | Hz    | Reduction rate                                            |
