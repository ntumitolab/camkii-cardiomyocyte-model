# Model descriptions

Hill equation.

$$
H(x, k, n) := \frac{x^n}{x^n + k^n}
$$

Relative exponential function. The definition here is reciprocal to the Python one.

$$
exprel(x) := \frac{x}{e^x - 1} \approx \frac{ln(1 + 2^{-x})}{ln2}
$$

GHK flux equation

$$
GHK(p, z, V_m, S_i, S_o) := pzF \cdot exprel(-zV_mF/RT) \cdot (S_i - S_o \exp(-zV_mF/RT))
$$

## Calcium diffusion

## Endoplasmic reticulum

### Ryanodine receptor (Jrel)

### SERCA (Jup)

## Sarcolemmal ion channels

Reversal potentials

$$
\begin{align}
E_{Na} &= \frac{RT}{F} \ln \frac{[Na^{+}]_o}{[Na^{+}]_i} \\
E_{K } &= \frac{RT}{F} \ln\frac{[K^{+}]_o}{[K^{+}]_i}  \\
E_{Ca} &= \frac{RT}{2F} \ln \frac{[Ca^{2+}]_o}{[Ca^{2+}]_i} \\
E_{Kr}  &= \frac{RT}{F} \ln\left( \frac{0.98 [K^{+}]_o + 0.02 [Na^{+}]_o}{0.98 [K^{+}]_i + 0.02[Na^{+}]_i} \right) \\
\end{align}
$$

### Sodium channels

Including fast sodium ($I_{Na}$) and background sodium ($I_{Na,b}$) currents.

$$
\begin{align}
I_{Na} &= \bar G_{Na} m_{Na}^{3} h_{Na} j_{Na} (V_m-E_{Na}) \\
I_{Na,b} &= \bar G_{Na,b} (V_m - E_{Na}) \\
\frac{dm_{Na}}{dt} &= \alpha_{m} - m_{Na}(\alpha_{m} + \beta_{m}) \\
\frac{dh_{Na}}{dt} &= \alpha_{h} - h_{Na}(\alpha_{h} + \beta_{h}) \\
\frac{dj_{Na}}{dt} &= \alpha_{j} - m_{Na}(\alpha_{j} + \beta_{j}) \\
\alpha_{m} &= 3.2 \text{ms}^{-1} exprel(-(V_m + 47.13) / 10)  \\
\beta_{m} &= 0.08 \text{ms}^{-1} \exp(-V_m / 11) \\
\\
\text{For} \ V_m & \ge -40\text{mV} \\
\alpha_{h} &= \alpha_{j} = 0 \text{ms}^{-1} \\
\beta_{h} &= (0.13 \text{ms} (1+\exp(-(V_m+10.66)/11.1)))^{-1} \\
\beta_{j} &= 0.3 \text{ms}^{-1} \frac{\exp(-2.535 \cdot 10^{-7}V_m)}{1 + \exp(-0.1(V_m + 32))} \\
\\
\text{For} \ V_m & < -40\text{mV} \\
\alpha_{h} &= 0.135 \text{ms}^{-1} \exp(-(V_m+80)/6.8) \\
\alpha_{j} &= (-127140 \exp(0.2444 V_m)-3.474 \cdot 10^{-5}\exp(-0.04391 V_m))\frac{V_m + 37.78}{1 + \exp(0.311( V_m + 79.23))} / \text{ms} \\
\beta_{h} &= (3.56\exp(0.079 V_m) + 3.1 \cdot 10^{5} \exp(0.35 V_m)) / \text{ms} \\
\beta_{j} &= \frac{0.1212 \exp(-0.01052 V_m)}{1 + \exp(-0.1378(V_m + 40.14))} / \text{ms} \\
\end{align}
$$

| Parameter  | Value  | Units | Description                            |
| ---------- | ------ | ----- | -------------------------------------- |
| $G_{Na}$   | 12.8   | mS/μF | Fast sodium channels conductance       |
| $G_{Na,b}$ | 0.0026 | mS/μF | Background sodium channels conductance |

### Potassium currents

$$
\begin{align}
\frac{d i_r }{dt} &= \frac{ r_∞  - i_r }{τ_r } \\
\frac{d i_s }{dt} &= \frac{ s_∞  - i_s }{τ_s } \\
\frac{d i_{s,slow} }{dt} &= \frac{ slow_∞ - i_{s,slow} }{τ_{s,slow} }
\\
\frac{d i_{nKs} }{dt} &= \frac{ nks_∞ - i_{nKs} }{τ_{nKs}} \\
\frac{d i_{CK1} }{dt} &= k_{b, IKr} i_{CK2}  - k_{f, IKr} i_{CK1}  + 0.022348 e^{0.01176 V_m } i_{CK0}  - 0.047002 e^{ - 0.0631 V_m } i_{CK1}  \\
\frac{d i_{CK2} }{dt} &= - k_{b, IKr} i_{CK2}  + k_{f, IKr} i_{CK1}  - 0.013733 i_{CK2}  e^{0.038198 V_m} + 6.89 \cdot 10^{-5} e^{ - 0.04178 V_m} i_{OK}  \\
\frac{d i_{OK} }{dt} &= 0.006497 i_{IK} e^{-0.03268 V_m} + 0.013733 i_{CK2}  e^{0.038198 V_m} - 6.89 \cdot 10^{-5} e^{-0.04178 V_m} i_{OK} - 0.090821 e^{0.023391 V_m} i_{OK}  \\
\frac{d i_{IK} }{dt} &=  - 0.006497 i_{IK}  e^{ - 0.03268 V_m } + 0.090821 e^{0.023391 V_m } i_{OK}  \\
\frac{d i_y}{dt} &= \frac{y_∞ - i_y}{τ_y} \\
I_{K1}  &= \frac{\mathtt{gK1} ( V_m - 6.1373 - E_K )}{0.1653 + \exp(0.0319 (V_m - 6.1373 - E_K))} \\
s_∞ &= \frac{1}{1 + \exp(0.21538 (31.972 + V_m))} \\
r_∞ &= \frac{1}{1 + e^{ - 0.068432 \left( -3.5572 + V_m  \right)}} \\
slow_∞ &= s_∞  \\
nks_∞ &= \frac{4.8133 \cdot 10^{-6} \left( 26.5 + V_m  \right)}{ - expm1\left(  - 0.128 \left( 26.5 + V_m  \right) \right) \left( \frac{4.8133 \cdot 10^{-6} \left( 26.5 + V_m  \right)}{ - expm1\left(  - 0.128 \left( 26.5 + V_m  \right) \right)} + 9.5333 \cdot 10^{-5} e^{ - 0.038 \left( 26.5 + V_m  \right)} \right)} \\
τ_r &= \frac{1000}{45.16 e^{0.03577 \left( 50 + V_m  \right)} + 98.9 e^{ - 0.1 \left( 38 + V_m  \right)}}
\\
τ_s &= 8.1 + 350 e^{ - \frac{1}{225} \left( 70 + V_m  \right)^{2}} \\
τ_{s,slow} &= 72.4 + 3700 e^{ - \frac{1}{900} \left( 70 + V_m  \right)^{2}} \\
I_{to}  &= G_t i_r (( 1 - f_{is} ) i_{s,slow} + f_{is} i_s )(
 V_m - E_K) \\
I_{Ks}  &= 2 G_{Ks} i_{nKs}^{2}( 0.68804 + 0.71283 \mathtt{IKUR\_PKAp}) (V_m - E_K) \\
I_{Kr}  &= G_{Kr} i_{OK} ( V_m - E_{Kr})   \\
1 &= i_{IK}  + i_{CK2}  + i_{OK}  + i_{CK0}  + i_{CK1}  \\
y_∞ &= \frac{1}{1 + e^{ - 0.15798 \left( -78.65 - V_m  \right)}} \\
τ_y  &= \frac{1000}{0.56236 e^{0.070472 \left( -75 - V_m  \right)} + 0.11885 e^{0.035249 ( 75 + V_m )}} \\
I_{fNa} &= f_{Na} G_f i_y ( V_m - E_{Na}) \\
I_{fK}  &= ( 1 - f_{Na}) G_f i_y ( V_m - E_K ) \\
I_{f}  &= I_{fK} + I_{fNa}
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

L-type calcium channels, T-type calcium channels, and background calcium currents

$$
\begin{align}
\frac{d \mathtt{i\_d} }{dt} &= \frac{\mathtt{dinf}  - \mathtt{i\_d} }{τ_d } \\
\frac{d \mathtt{i\_f} }{dt} &= \frac{\mathtt{finf}  - \mathtt{i\_f} }{τ_f } \\
\frac{d \mathtt{i\_fca} }{dt} &= \frac{\left( \mathtt{fcainf}  - \mathtt{i\_fca}  \right) \left( 1 - \left( \mathtt{fcainf}  > \mathtt{i\_fca}  \right) \left( \mathtt{Vm}  > -60 \right) \right)}{τ_{fca}} \\
\frac{d \mathtt{i\_b} }{dt} &= \frac{ - \mathtt{i\_b}  + \mathtt{binf} }{\mathtt{taub} } \\
\frac{d \mathtt{i\_g} }{dt} &= \frac{\mathtt{ginf}  - \mathtt{i\_g} }{\mathtt{taug} } \\
\mathtt{ICa\_scale}  &= \mathtt{ICa\_scale0} \left( 1 + \frac{-0.56}{-1 + \frac{\mathtt{fracLCCbpISO}}{\mathtt{fracLCCbp0}}} \right) \\
\mathtt{E\_Ca}  &= 13.356 \ln\left( \frac{\mathtt{ca\_o}}{\mathtt{ca\_i}} \right) \\
\mathtt{INaCa}  &= \frac{\left([Na^{+}]_i^{3} \mathtt{ca\_o} e^{0.037436 \mathtt{gNaCa} \mathtt{Vm} } - [Na^{+}]_o^{3} \mathtt{ca\_i} \mathtt{fNaCa} e^{0.037436 \left( -1 + \mathtt{gNaCa} \right) \mathtt{Vm} } \right) \mathtt{kNaCa} \mathtt{ICa\_scale} }{1 + \left([Na^{+}]_i^{3} \mathtt{ca\_o} + [Na^{+}]_o^{3} \mathtt{ca\_i} \mathtt{fNaCa} \right) \mathtt{dNaCa}} \\
\mathtt{ICaL}  &= \frac{14448 \mathtt{GCaL} \left(  - 0.341 \mathtt{ca\_o} + \mathtt{ca\_i} \left( 1 + expm1\left( 0.074872 \mathtt{Vm}  \right) \right) \right) \mathtt{ICa\_scale}  \mathtt{Vm}  \mathtt{i\_f}  \mathtt{i\_fca}  \mathtt{i\_d} }{expm1\left( 0.074872 \mathtt{Vm}  \right)} \\
\mathtt{dinf}  &= \frac{1}{1 + e^{ - 0.13889 \left( 11.1 + \mathtt{Vm}  \right)}} \\
τ_d  &= \frac{1.4 \left( 0.25 + \frac{1.4}{1 + e^{ - \frac{1}{13} \left( 35 + \mathtt{Vm}  \right)}} \right)}{1 + e^{ - \frac{1}{5} \left( -5
- \mathtt{Vm}  \right)}} + \frac{1}{1 + e^{ - \frac{1}{20} \left( -50 + \mathtt{Vm}  \right)}} \\
\mathtt{finf}  &= \frac{1}{1 + e^{ - 0.18519 \left( -23.3 - \mathtt{Vm}  \right)}} \\
τ_f  &= 120 + \frac{165}{1 + e^{ - \frac{1}{10} \left( -25 + \mathtt{Vm}  \right)}} + 1125 e^{ - \frac{1}{240} \left( 27 + \mathtt{Vm}  \right)^{2}} \\
\mathtt{fcainf}  &= 0.68493 \left( 0.23 + \frac{0.00319}{0.00319 + \mathtt{ca\_i}^{8}} + \frac{0.1}{1 + e^{ - 10 \left( 0.5 - \mathtt{ca\_i} \right)}} + \frac{0.2}{1 + e^{ - 1.25 \left( 0.75 - \mathtt{ca\_i} \right)}} \right) \\
\mathtt{binf}  &= \frac{1}{1 + e^{ - 0.18497 \left( 37.491 + \mathtt{Vm}  \right)}} \\
\mathtt{taub}  &= 0.6 + \frac{5.4}{1 + e^{ - 0.03 \left( -100 - \mathtt{Vm}  \right)}} \\
\mathtt{ginf}  &= \frac{1}{1 + e^{ - \frac{1}{6} \left( -66 - \mathtt{Vm}  \right)}} \\
\mathtt{taug}  &= 1 + \frac{40}{1 + e^{ - 0.08 \left( -65 - \mathtt{Vm}  \right)}} \\
\mathtt{ICaT}  &= \mathtt{gCaT} \left( 106.5 - \mathtt{E\_Ca}  + \mathtt{Vm}  \right) \mathtt{i\_b}  \mathtt{i\_g}  \\
\mathtt{ICab}  &= \mathtt{gCab} \left(  - \mathtt{E\_Ca}  + \mathtt{Vm}  \right)
\end{align}
$$

| Parameter | Value         | Units      | Description |
| --------- | ------------- | ---------- | ----------- |
| fNaCa     | 1             | -          |             |
| kNaCa     | 2.268 * 10⁻¹⁶ | μAμF⁻¹μM⁻⁴ |             |
| dNaCa     | 10⁻¹⁶         | μM⁻⁴       |             |
| gNaCa     | 0.5           |            |             |
| $G_{CaL}$ | 6.3 * 10⁻⁵    | m³s⁻¹F⁻¹   |             |
| $τ_{fca}$ | 10            | ms         |             |
| $g_{CaT}$ | 0.2           | mS/μF      |             |
| $g_{Cab}$ | 0.0008        | mS/μF      |             |

### Na-K pump

$$
\begin{align}
I_{NaK} &= I_{NaK}^{max} fNaK \frac{[K^+]_o}{ [K^+]_o + KmKo_{NaK} } \frac{[Na^+]_i^{nNaK}}{ [Na^+]_i^{nNaK} + KmNai_{NaK}^{nNaK} } \\
fNaK &= (1 + 0.1245 \exp(-0.1 V_m F / RT) + 0.0365 sigma \exp(V_m F / RT))^{-1} \\
sigma &= 1 / 7 \cdot (\exp([Na^+]_o / 67.3 \text{mM}) - 1)
\end{align}
$$

| Parameter       | Value | Units | Description               |
| --------------- | ----- | ----- | ------------------------- |
| $I_{NaK}^{max}$ | 2.7   | μA/μF | Maximal rate of Na-K pump |
| $KmNai_{NaK}$   | 18.6  | mM    |                           |
| $KmKo_{NaK}$    | 1.5   | mM    |                           |
| $nNaK$          | 3.2   | -     |                           |

## Beta-adrenergic system

$$
\begin{align}
f_{PKACI} &= PKACI_0 + \frac{ PKACI_{act} ISO }{ISO + PKACI_{KM}} \\
f_{PKACII} &= PKACII_0 + \frac{ PKACII_{act} ISO }{ISO + PKACII_{KM}} \\
\end{align}
$$

## CaMKII system

## ODE system

$$
\begin{align}
\frac{d[Na^+]_i}{dt} &= -(I_{Na} + 3I_{NaCa} + 3I_{NaK})\frac{A_{cap}}{V_{myo}F} + (V_{NHE} - 3V_{NaCa}) \frac{V_{mito}}{V_{myo}} \\
\frac{d[K^+]_i}{dt} &= -(I_{Ks} + I_{Kr} + I_{K1} + I_{Kp} + I_{Ca,K}-2I_{NaK})\frac{A_{cap}}{V_{myo}F} \\
C_m\frac{dV_m}{dt} &= -(I_{Na} + I_{CaL} + I_{Kr} + I_{Ks} + I_{K1} + I_{Kp} + I_{NaCa} + I_{NaK} + I_{pCa} + I_{Ca, b} + I_{K_{ATP}} + I_{stim}) \\
β_i &= \frac{(K_m^{CMDN} + [Ca^{2+}]_i)^2}{ (K_m^{CMDN} + [Ca^{2+}]_i)^2 + K_m^{CMDN}  \cdot  [CMDN]_{tot}} \\
β_{SR} &= \frac{(K_m^{CSQN} + [Ca^{2+}]_{SR})^2}{(K_m^{CSQN} + [Ca^{2+}]_{SR})^2 + K_m^{CSQN}  \cdot  [CSQN]_{tot}} \\
\frac{d[Ca^{2+}]_i}{dt} &= \beta_i(J_{xfer}\frac{V_{ss}}{V_{myo}} - J_{up} - J_{trpn} - (I_{Ca,b} -2I_{NaCa} + I_{pCa})\frac{A_{cap}}{2V_{myo}F} + (V_{NaCa} - V_{uni})\frac{V_{mito}}{V_{myo}}) \\
\frac{d[Ca^{2+}]_{SR}}{dt} &= \beta_{SR}(J_{up}\frac{V_{myo}}{V_{SR}} - J_{rel}\frac{V_{ss}}{V_{SR}}) \\
\end{align}
$$
