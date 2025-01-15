# Model descriptions

Hill equation.

$$
H(x, k) := \frac{x}{x + k}
$$

$$
H(x, k, n) := \frac{x^n}{x^n + k^n}
$$

Relative exponential function. The definition here is reciprocal to the Python one.

$$
exprel(x) := \frac{x}{ \exp(x) - 1}
$$

Logistic function

$$
expit(x) := \frac{1}{1 + \exp(-x)}
$$

Thermal voltage.

$$
V_T = \frac{RT}{F}
$$

GHK flux equation.

$$
GHK(p, z, V_m, S_i, S_o) := pzF \cdot exprel(-zV_m/V_T) \cdot (S_i - S_o \exp(-zV_m/V_T))
$$

Ion concentrations:

$$
\begin{align}
\mathrm{na_x} &= \mathrm{[Na^+]_x} \\
\mathrm{k_x} &= \mathrm{[K^+]_x} \\
\mathrm{ca_x} &= \mathrm{[Ca^{2+}]_x} \\
\end{align}
$$

Reversal potentials:

$$
\begin{align}
E_{Na} &= V_T \ln \frac{\mathrm{na_o}}{\mathrm{na_i}} \\
E_{K } &= V_T \ln\frac{\mathrm{k_o}}{\mathrm{k_i}}  \\
E_{Ca} &= 0.5V_T \ln \frac{\mathrm{ca_o}}{\mathrm{ca_{sl}}} \\
E_{Kr}  &= V_T \ln\left( \frac{0.98 \mathrm{k_o} + \mathrm{na_o}}{0.98 \mathrm{k_i} + 0.02\mathrm{na_i}} \right) \\
\end{align}
$$

## Calcium diffusion

Cytosolic calcium is diffused between sub-sarcolemma (SL) and sub-sarcoplasmic (SR) spaces.

Calcium buffering in each compartment:

$$
\beta_{Ca,x} = (1 + \Sigma Trpn \cdot Km_{Trpn, 2} / (ca_x + Km_{Trpn, 2})^2 + \Sigma Cmdn \cdot Km_{Cmdn} / (ca_x + Km_{Cmdn, 2})^2 )^{-1}
$$




## Endoplasmic reticulum

### Ryanodine receptor (Jrel)

### SERCA (Jup)

## Sarcolemmal ion channels

### Sodium channels

Including fast sodium ($\mathrm{I_{Na}}$) and background sodium ($\mathrm{I_{Na,b}}$) currents.

$$
\begin{align}
\mathrm{I_{Na}} &= \bar G_{Na} \cdot m_{Na}^{3} \cdot h_{Na} \cdot j_{Na} (V_m-E_{Na}) \\
\mathrm{I_{Na,b}} &= \bar G_{Na,b} (V_m - E_{Na}) \\
\frac{dm_{Na}}{dt} &= \alpha_{m} - m_{Na}(\alpha_{m} + \beta_{m}) \\
\frac{dh_{Na}}{dt} &= \alpha_{h} - h_{Na}(\alpha_{h} + \beta_{h}) \\
\frac{dj_{Na}}{dt} &= \alpha_{j} - m_{Na}(\alpha_{j} + \beta_{j}) \\
\alpha_{m} &= 3.2 \mathrm{ms}^{-1} exprel(-(V_m + 47.13) / 10)  \\
\beta_{m} &= 0.08 \mathrm{ms}^{-1} \exp(-V_m / 11) \\
\\
\text{For} \ V_m & \ge -40\text{mV} \\
\alpha_{h} &= \alpha_{j} = 0 \\
\beta_{h} &= 7.6923 \mathrm{ms}^{-1} expit((V_m+10.66)/11.1) \\
\beta_{j} &= 0.3 \mathrm{ms}^{-1} \exp(-2.535 \cdot 10^{-7}V_m) expit(0.1(V_m + 32)) \\
\\
\text{For} \ V_m & < -40\text{mV} \\
\alpha_{h} &= 0.135 \text{ms}^{-1} \exp(-(V_m+80)/6.8) \\
\alpha_{j} &= (-127140 \exp(0.2444 V_m)-3.474 \cdot 10^{-5}\exp(-0.04391 V_m))\frac{V_m + 37.78}{1 + \exp(0.311( V_m + 79.23))} / \text{ms} \\
\beta_{h} &= (3.56\exp(0.079 V_m) + 3.1 \cdot 10^{5} \exp(0.35 V_m)) / \text{ms} \\
\beta_{j} &= \frac{0.1212 \mathrm{ms}^{-1} \exp(-0.01052 V_m)}{1 + \exp(-0.1378(V_m + 40.14))} \\
\end{align}
$$

| Parameter  | Value  | Units | Description                            |
| ---------- | ------ | ----- | -------------------------------------- |
| $G_{Na}$   | 12.8   | mS/μF | Fast sodium channels conductance       |
| $G_{Na,b}$ | 0.0026 | mS/μF | Background sodium channels conductance |

### Potassium currents

$$
\begin{align}
\mathrm{I_{K1}}  &= \mathrm{G_{K1}} \cdot H(\mathrm{k_o}, 210)  \frac{ ( V_m - 6.1373 - E_K )}{0.1653 + \exp(0.0319 (V_m - 6.1373 - E_K))} \\
\mathrm{I_{to}}  &= G_t \cdot i_r (( 1 - f_{is} ) i_{s,slow} + f_{is} i_s )(
 V_m - E_K) \\
\mathrm{I_{Ks}}  &= 2 G_{Ks} \cdot i_{nKs}^{2}( 0.68804 + 0.71283 \mathrm{IKUR_{PKAp}}) (V_m - E_K) \\
\mathrm{I_{Kr}}  &= G_{Kr} \cdot i_{OK} (V_m - E_{Kr})   \\
\mathrm{I_{fNa}} &= f_{Na} \cdot G_f \cdot i_y \cdot ( V_m - E_{Na}) \\
\mathrm{I_{fK}}  &= ( 1 - f_{Na}) \cdot G_f \cdot i_y \cdot ( V_m - E_K ) \\
\mathrm{I_{f}}  &= \mathrm{I_{fK}} + \mathrm{I_{fNa}} \\
\frac{d i_r }{dt} &= \frac{ r_∞  - i_r }{τ_r} \\
\frac{d i_s }{dt} &= \frac{ s_∞  - i_s }{τ_s} \\
\frac{d i_{s,slow} }{dt} &= \frac{ slow_∞ - i_{s,slow} }{τ_{s,slow} }
\\
\frac{d i_{nKs} }{dt} &= \frac{ nks_∞ - i_{nKs} }{τ_{nKs}} \\
\frac{d i_{CK1} }{dt} &= k_{b, IKr} i_{CK2}  - k_{f, IKr} i_{CK1}  + 0.022348 e^{0.01176 V_m } i_{CK0}  - 0.047002 e^{ - 0.0631 V_m } i_{CK1}  \\
\frac{d i_{CK2} }{dt} &= - k_{b, IKr} i_{CK2}  + k_{f, IKr} i_{CK1}  - 0.013733 i_{CK2}  e^{0.038198 V_m} + 6.89 \cdot 10^{-5} e^{ - 0.04178 V_m} i_{OK}  \\
\frac{d i_{OK} }{dt} &= 0.006497 i_{IK} e^{-0.03268 V_m} + 0.013733 i_{CK2}  e^{0.038198 V_m} - 6.89 \cdot 10^{-5} e^{-0.04178 V_m} i_{OK} - 0.090821 e^{0.023391 V_m} i_{OK}  \\
\frac{d i_{IK} }{dt} &=  - 0.006497 i_{IK}  e^{ - 0.03268 V_m } + 0.090821 e^{0.023391 V_m } i_{OK}  \\
\frac{d i_y}{dt} &= \frac{y_∞ - i_y}{τ_y} \\
s_∞ &= expit((V_m + 31.97156) / -4.64291) \\
r_∞ &= expit((V_m - 3.55716) / 14.61299) \\
slow_∞ &= s_∞  \\
nks_∞ &= \frac{\alpha_{nks}}{\alpha_{nks} + \beta_{nks}} \\
\alpha_{nks} &= 0.00000481333 / 0.128 * exprel(-0.128 * (V + 26.5)) \\
\beta_{nks} &= 0.0000953333 * \exp(-0.038 * (V + 26.5)) \\
τ_r &= \frac{1000 \text{ms}}{45.16 \exp(0.03577( V_m + 50)) + 98.9 \exp( - 0.1 ( V_m + 38))}
\\
τ_s &= 8.1\text{ms} + 350\text{ms} \cdot \exp( - \frac{1}{225}(V_m + 70)^{2}) \\
τ_{s,slow} &= 72.4\text{ms} + 3700\text{ms} \cdot \exp( - \frac{1}{900}( V_m + 70 )^{2}) \\
1 &= i_{IK}  + i_{CK2}  + i_{OK}  + i_{CK0}  + i_{CK1}  \\
y_∞ &= expit(-0.15798 (V_m + 78.65)) \\
τ_y  &= \frac{1000 \text{ms}}{0.56236 \exp(-0.070472 (V_m + 75)) + 0.11885 \exp(0.035249 ( V_m + 75))} \\
\end{align}
$$

| Parameter    | Value    | Units | Description                                      |
| ------------ | -------- | ----- | ------------------------------------------------ |
| $G_{K1}$     | 0.0515   | mS/μF | Potassium channels conductance                   |
| $G_t$        | 0.1      | mS/μF | Transient outward potassium channels conductance |
| $G_{Ks}$     | 0.05     | mS/μF | Potassium channels conductance                   |
| $τ_{nKs}$    | 750      | ms    | Potassium channels time scale                    |
| $G_{Kr}$     | 0.06     | mS/μF | Potassium channels conductance                   |
| $k_{f, IKr}$ | 0.023761 | 1/ms  | Potassium channels transition rate               |
| $k_{b, IKr}$ | 0.036778 | 1/ms  | Potassium channels transition rate               |
| $G_f$        | 0.021    | mS/μF | Funny current conductance                        |
| $f_{Na}$     | 0.021    | -     | Funny current sodium fraction                    |
| $f_{is}$      | 0.706    | -     |                                                  |

### Calcium currents

L-type calcium channels, T-type calcium channels, and background calcium currents.

$$
\begin{align}
\mathrm{I_{CaL}}  &= \mathrm{ICa_{scale}} \cdot G_{CaL} \cdot i_d \cdot i_f \cdot i_fca \cdot GHK(G_{CaL}, 2, V_m, \mathrm{ca_{sl}}, 0.341 \mathrm{ca_o}) \\
\mathrm{I_{CaT}} &= \mathrm{gCaT} \cdot i_b \cdot i_g ( V_m + 106.5 - \mathrm{E_Ca}) \\
\mathrm{I_{Cab}} &= \mathrm{gCab} (V_m - \mathrm{E_Ca}) \\
\frac{d i_d }{dt} &= \frac{d_∞  - i_d}{τ_d} \\
\frac{d i_f }{dt} &= \frac{f_∞  - i_f}{τ_f} \\
\frac{d i_{fca} }{dt} &= \frac{(fca_∞ - i_{fca}) \left( 1 - \left( fca_∞  > i_{fca}  \right) \left( V_m  > -60 \text{mV} \right) \right)}{τ_{fca}} \\
\frac{d i_b }{dt} &= \frac{b_∞ - i_b}{τ_b} \\
\frac{d i_g }{dt} &= \frac{g_∞ - i_g}{τ_g} \\
\mathrm{ICa_{scale}}  &= \mathrm{ICa_{scale, 0}} \left( 1 + \frac{0.56}{1 - \frac{\mathrm{fracLCCbpISO}}{\mathrm{fracLCCbp0}}} \right) \\
\mathrm{I_{NaCa}}  &= \mathrm{kNaCa} \cdot \mathrm{ICa_{scale}} \frac{\mathrm{na_i}^{3} \mathrm{ca_o} \exp( \mathrm{gNaCa} V_m /V_T ) - \mathrm{na_o}^{3} \mathrm{ca_{sl}} \mathrm{fNaCa} \exp (( \mathrm{gNaCa} - 1 ) V_m F /RT )} {1 + \left(\mathrm{na_i}^3 \mathrm{ca_o} + \mathrm{na_o}^{3} \mathrm{ca_{sl}} \mathrm{fNaCa} \right) \mathrm{dNaCa}} \\
d_∞ &= expit((V + 11.1) / 7.2) \\
τ_d &= (\alpha_d \beta_d + \gamma_d) \\
\alpha_d &= 1.4 expit((V_m + 35) / 13) + 0.25 \\
\beta_d &= 1.4 expit(-(V_m + 5) / 5) \\
\gamma_d &= expit((V_m - 50) / 20) \\
f_∞ &= expit(-(V_m + 23.3) / 5.4) \\
τ_f &= 120 + 165 * expit((V_m - 25) / 10) + 1125 \exp( -(V_m + 27)^{2} / 240) \\
fca_∞  &= (\alpha_{fca} + \beta_{fca} + \gamma_{fca} + 0.23) / 1.46 \\
\alpha_{fca} &= H(0.4875, \mathrm{ca_{sl}}, 8) \\
\beta_{fca} &= 0.1 expit(-(\mathrm{ca_{sl}} - 0.5) / 0.1) \\
\gamma_{fca} &= 0.2 expit(-(\mathrm{ca_{sl}} - 0.75) / 0.8) \\
b_∞ &= expit((V_m + 37.49098) / 5.40634) \\
τ_b &= 0.6 + 5.4 expit(-0.03 (V_m + 100)) \\
g_∞ &= expit(-(V_m + 66) / 6) \\
τ_g &= 1 + 40 expit(-0.08 (V_m + 65)) \\
\end{align}
$$

| Parameter | Value         | Units      | Description |
| --------- | ------------- | ---------- | ----------- |
| fNaCa     | 1             | -          |             |
| kNaCa     | 2.268 * 10⁻¹⁶ | μAμF⁻¹μM⁻⁴ |             |
| dNaCa     | 10⁻¹⁶         | μM⁻⁴       |             |
| gNaCa     | 0.5           | -          |             |
| $G_{CaL}$ | 6.3 * 10⁻⁵    | m³s⁻¹F⁻¹   |             |
| $τ_{fca}$ | 10            | ms         |             |
| $g_{CaT}$ | 0.2           | mSμF⁻¹     |             |
| $g_{Cab}$ | 0.0008        | mSμF⁻¹     |             |

### Na-K pump

$$
\begin{align}
\mathrm{I_{NaK}} &= I_{NaK}^{max} fNaK \frac{\mathrm{k_o}}{ \mathrm{k_o} + KmKo_{NaK} } \frac{\mathrm{na_i}^{nNaK}}{ \mathrm{na_i}^{nNaK} + KmNai_{NaK}^{nNaK} } \\
fNaK &= (1 + 0.1245 \exp(-0.1 V_m/ V_T) + 0.0365 \sigma_{NaK} \exp(V_m/V_T))^{-1} \\
\sigma_{NaK} &= (\exp(\mathrm{na_i} / 67.3 \text{mM}) - 1) / 7
\end{align}
$$

| Parameter       | Value | Units | Description               |
| --------------- | ----- | ----- | ------------------------- |
| $I_{NaK}^{max}$ | 2.7   | μA/μF | Maximal rate of Na-K pump |
| $KmNai_{NaK}$   | 18.6  | mM    |                           |
| $KmKo_{NaK}$    | 1.5   | mM    |                           |
| $nNaK$          | 3.2   | -     |                           |

## Beta-adrenergic system

Activities are fitted to the steady-state activities in the Morroti model.

$$
\begin{align}
f_{PKACI} &= PKACI_0 + PKACI_{act} H(ISO, PKACI_{KM}) \\
f_{PKACII} &= PKACII_0 + PKACII_{act} H(ISO, PKACII_{KM}) \\
f_{PP1} &= PP1_0 + PP1_{act} H(PP1_{KI}, ISO) \\
f_{PLBp} &= PLBp_0 + PLBp_{act} H(ISO, PLBp_{KM}, PLBp_{nHill}) \\
f_{PLMp} &= PLMp_0 + PLMp_{act} H(ISO, PLMp_{KM}, PLMp_{nHill}) \\
TnI_{PKAp} &= TnIp_0 + TnIp_{act} H(ISO, TnIp_{KM}, TnIp_{nHill}) \\
LCCa_{PKAp} &= LCCap_0 + LCCap_{act} H(ISO, LCCap_{KM}) \\
LCCb_{PKAp} &= LCCbp_0 + LCCbp_{act} H(ISO, LCCbp_{KM}) \\
KUR_{PKAp} &= KURp_0 + KURp_{act} H(ISO, KURp_{KM}) \\
RyR_{PKAp} &= RyRp_0 + RyRp_{act} H(ISO, RyRp_{KM}) \\
\end{align}
$$

| Parameter      | Value     | Units | Description               |
| -------------- | --------- | ----- | ------------------------- |
| $PKACI_0$      | 0.0734    | -     | Basal PKACI activity      |
| $PKACI_{act}$  | 0.1995    | -     | Activated PKACI activity  |
| $PKACI_{KM}$   | 0.0139    | μM    | PKACI sesitivity to iso   |
| $PKACII_0$     | 0.1840    | -     | Basal PKACII activity     |
| $PKACII_{act}$ | 0.3444    | -     | Activated PKACII activity |
| $PKACII_{KM}$  | 0.0103    | μM    | PKACII sesitivity to iso  |
| $PP1_0$        | 0.8927    | -     | Basal PP1 activity        |
| $PP1_{act}$    | 0.0492    | -     | Activated PP1 activity    |
| $PP1_{KI}$     | 0.00637   | μM    | PP1 sesitivity to iso     |
| $PLBp_0$       | 0.0824    | -     |                           |
| $PLBp_{act}$   | 0.7961    | -     |                           |
| $PLBp_{KM}$    | 0.00597   | μM    |                           |
| $PLBp_{nHill}$ | 1.8167    | -     |                           |
| $PLMp_0$       | 0.1172    | -     |                           |
| $PLMp_{act}$   | 0.6645    | -     |                           |
| $PLMp_{KM}$    | 0.00823   | μM    |                           |
| $PLMp_{nHill}$ | 1.35784   | -     |                           |
| $TnIp_0$       | 0.0669    | -     |                           |
| $TnIp_{act}$   | 0.7524    | -     |                           |
| $TnIp_{KM}$    | 0.007913  | μM    |                           |
| $TnIp_{nHill}$ | 1.6736    | -     |                           |
| $LCCap_0$      | 0.2205    | -     |                           |
| $LCCap_{act}$  | 0.2339    | -     |                           |
| $LCCap_{KM}$   | 0.00726   | μM    |                           |
| $LCCbp_0$      | 0.2517    | -     |                           |
| $LCCbp_{act}$  | 0.2461    | -     |                           |
| $LCCbp_{KM}$   | 0.00695   | μM    |                           |
| $KURp_0$       | 0.4390    | -     |                           |
| $KURp_{act}$   | 0.2563    | -     |                           |
| $KURp_{KM}$    | 0.00557   | μM    |                           |
| $RyRp_0$       | 0.2054    | -     |                           |
| $RyRp_{act}$   | 0.2399    | -     |                           |
| $RyRp_{KM}$    | 0.0075135 | μM    |                           |

## CaMKII system



## ODE system

$$
\begin{align}
\frac{d\mathrm{na_i}}{dt} &= -(I_{Na} + 3I_{NaCa} + 3I_{NaK})\frac{A_{cap}}{V_{myo}F} + (V_{NHE} - 3V_{NaCa}) \frac{V_{mito}}{V_{myo}} \\
\frac{d[K^+]_i}{dt} &= -(I_{Ks} + I_{Kr} + I_{K1} + I_{Kp} + I_{Ca,K}-2I_{NaK})\frac{A_{cap}}{V_{myo}F} \\
C_m\frac{dV_m}{dt} &= -(I_{Na} + I_{CaL} + I_{Kr} + I_{Ks} + I_{K1} + I_{Kp} + I_{NaCa} + I_{NaK} + I_{pCa} + I_{Ca, b} + I_{K_{ATP}} + I_{stim}) \\
β_i &= \frac{(K_m^{CMDN} + [Ca^{2+}]_i)^2}{ (K_m^{CMDN} + [Ca^{2+}]_i)^2 + K_m^{CMDN}  \cdot  [CMDN]_{tot}} \\
β_{SR} &= \frac{(K_m^{CSQN} + [Ca^{2+}]_{SR})^2}{(K_m^{CSQN} + [Ca^{2+}]_{SR})^2 + K_m^{CSQN}  \cdot  [CSQN]_{tot}} \\
\frac{d[Ca^{2+}]_i}{dt} &= \beta_i(J_{xfer}\frac{V_{ss}}{V_{myo}} - J_{up} - J_{trpn} - (I_{Ca,b} -2I_{NaCa} + I_{pCa})\frac{A_{cap}}{2V_{myo}F} + (V_{NaCa} - V_{uni})\frac{V_{mito}}{V_{myo}}) \\
\frac{d[Ca^{2+}]_{SR}}{dt} &= \beta_{SR}(J_{up}\frac{V_{myo}}{V_{SR}} - J_{rel}\frac{V_{ss}}{V_{SR}}) \\
\end{align}
$$
