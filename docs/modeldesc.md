# Model descriptions

## Calcium diffusion

## Endoplasmic reticulum

### Ryanodine receptor (Jrel)

### SERCA (Jup)

## Sarcolemmal ion channels

### Sodium channels

$$
\begin{align}
I_{Na} &= \bar G_{Na} m_{Na}^{3} h_{Na} j_{Na} (V_m-E_{Na}) \\
I_{Na,b} &= \bar G_{Na,b} (V_m - E_{Na}) \\
E_{Na} &= \frac{RT}{F} \ln \frac{[Na^{+}]_o}{[Na^{+}]_i} \\
\frac{dm_{Na}}{dt} &= \alpha_{m} - m_{Na}(\alpha_{m} + \beta_{m}) \\
\frac{dh_{Na}}{dt} &= \alpha_{h} - h_{Na}(\alpha_{h} + \beta_{h}) \\
\frac{dj_{Na}}{dt} &= \alpha_{j} - m_{Na}(\alpha_{j} + \beta_{j}) \\
\alpha_{m} &= 0.32 \text{ms}^{-1} \frac{V_m + 47.13}{1 - \exp(-0.1(V_m+47.13))} \\
\beta_{m} &= 0.08 \text{ms}^{-1} \exp(-V_m / 11) \\
\\
\text{For} \ V_m & \ge -40\text{mV} \\
\alpha_{h} &= \alpha_{j} = 0 \\
\beta_{h} &= (0.13 \text{ms} (1+\exp(-(V_m+10.66)/11.1)))^{-1} \\
\beta_{j} &= 0.3 \text{ms}^{-1} \frac{\exp(-2.535 \times 10^{-7}V_m)}{1 + \exp(-0.1(V_m + 32))} \\
\\
\text{For} \ V_m & < -40\text{mV} \\
\alpha_{h} &= 0.135 \text{ms}^{-1} \exp(-(V_m+80)/6.8) \\
\alpha_{j} &= (-127140 \exp(0.2444 V_m)-3.474 \times 10^{-5}\exp(-0.04391 V_m))\frac{V_m + 37.78}{1 + \exp(0.311( V_m + 79.23))} / \text{ms} \\
\beta_{h} &= (3.56\exp(0.079 V_m) + 3.1 \times 10^{5} \exp(0.35 V_m)) / \text{ms} \\
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
\frac{\mathrm{d} \mathtt{i\_r}\left( t \right)}{\mathrm{d}t} &= \frac{\mathtt{rinf}\left( t \right) - \mathtt{i\_r}\left( t \right)}{\mathtt{taur}\left( t \right)} \\
\frac{\mathrm{d} \mathtt{i\_s}\left( t \right)}{\mathrm{d}t} &= \frac{\mathtt{sinf}\left( t \right) - \mathtt{i\_s}\left( t \right)}{\mathtt{taus}\left( t \right)} \\
\frac{\mathrm{d} \mathtt{i\_sslow}\left( t \right)}{\mathrm{d}t} &= \frac{ - \mathtt{i\_sslow}\left( t \right) + \mathtt{slowinf}\left( t \right)}{\mathtt{tausslow}\left( t \right)}
\\
\frac{\mathrm{d} \mathtt{i\_nKs}\left( t \right)}{\mathrm{d}t} &= \frac{\frac{4.8133 \cdot 10^{-6} \left( 26.5 + \mathtt{V_m}\left( t \right) \right)}{ - expm1\left(  - 0.128 \left( 26.5 + \mathtt{V_m}\left( t \right) \right) \right) \left( \frac{4.8133 \cdot 10^{-6} \left( 26.5 + \mathtt{V_m}\left( t \right) \right)}{ - expm1\left(  - 0.128 \left( 26.5 + \mathtt{V_m}\left( t \right) \right) \right)} + 9.5333 \cdot 10^{-5} e^{ - 0.038 \left( 26.5 + \mathtt{V_m}\left( t \right) \right)} \right)} - \mathtt{i\_nKs}\left( t \right)}{\mathtt{nKstau}} \\
\frac{\mathrm{d} \mathtt{i\_CK1}\left( t \right)}{\mathrm{d}t} &= \mathtt{kb} \mathtt{i\_CK2}\left( t \right) - \mathtt{kf} \mathtt{i\_CK1}\left( t \right) + 0.022348 e^{0.01176 \mathtt{V_m}\left( t \right)} \mathtt{CK0}\left( t \right) - 0.047002 e^{ - 0.0631 \mathtt{V_m}\left( t \right)} \mathtt{i\_CK1}\left( t \right) \\
\frac{\mathrm{d} \mathtt{i\_CK2}\left( t \right)}{\mathrm{d}t} &=  - \mathtt{kb} \mathtt{i\_CK2}\left( t \right) + \mathtt{kf} \mathtt{i\_CK1}\left( t \right) - 0.013733 \mathtt{i\_CK2}\left( t \right) e^{0.038198 \mathtt{V_m}\left( t \right)} + 6.89 \cdot 10^{-5} e^{ - 0.04178 \mathtt{V_m}\left( t \right)} \mathtt{i\_OK}\left( t \right) \\
\frac{\mathrm{d} \mathtt{i\_OK}\left( t \right)}{\mathrm{d}t} &= 0.006497 \mathtt{i\_IK}\left( t \right) e^{ - 0.03268 \mathtt{V_m}\left( t \right)} + 0.013733 \mathtt{i\_CK2}\left( t \right) e^{0.038198 \mathtt{V_m}\left( t \right)} - 6.89 \cdot 10^{-5} e^{ - 0.04178 \mathtt{V_m}\left( t \right)} \mathtt{i\_OK}\left( t \right) - 0.090821 e^{0.023391 \mathtt{V_m}\left( t \right)} \mathtt{i\_OK}\left( t \right) \\
\frac{\mathrm{d} \mathtt{i\_IK}\left( t \right)}{\mathrm{d}t} &=  - 0.006497 \mathtt{i\_IK}\left( t \right) e^{ - 0.03268 \mathtt{V_m}\left( t \right)} + 0.090821 e^{0.023391 \mathtt{V_m}\left( t \right)} \mathtt{i\_OK}\left( t \right) \\
\frac{\mathrm{d} \mathtt{i\_y}\left( t \right)}{\mathrm{d}t} &= \frac{\mathtt{yinf}\left( t \right) - \mathtt{i\_y}\left( t \right)}{\mathtt{tauy}\left( t \right)} \\
\mathtt{E\_Na}\left( t \right) &= \frac{RT}{F} \ln\left( \frac{\mathtt{na\_o}}{\mathtt{na\_i}} \right) \\
\mathtt{E\_K}\left( t \right) &= \frac{RT}{F} \ln\left( \frac{\mathtt{k\_o}}{\mathtt{k\_i}} \right) \\
\mathtt{IK1}\left( t \right) &= \frac{\mathtt{gK1} \left( -6.1373 - \mathtt{E\_K}\left( t \right) + \mathtt{V_m}\left( t \right) \right)}{0.1653 + e^{0.0319 \left( -6.1373 - \mathtt{E\_K}\left( t \right) + \mathtt{V_m}\left( t \right) \right)}} \\
\mathtt{sinf}\left( t \right) &= \frac{1}{1 + e^{0.21538 \left( 31.972 + \mathtt{V_m}\left( t \right) \right)}} \\
\mathtt{rinf}\left( t \right) &= \frac{1}{1 + e^{ - 0.068432 \left( -3.5572 + \mathtt{V_m}\left( t \right) \right)}} \\
\mathtt{slowinf}\left( t \right) &= \mathtt{sinf}\left( t \right) \\
\mathtt{taur}\left( t \right) &= \frac{1000}{45.16 e^{0.03577 \left( 50 + \mathtt{V_m}\left( t \right) \right)} + 98.9 e^{ - 0.1 \left( 38 + \mathtt{V_m}\left( t \right) \right)}}
\\
\mathtt{taus}\left( t \right) &= 8.1 + 350 e^{ - \frac{1}{225} \left( 70 + \mathtt{V_m}\left( t \right) \right)^{2}} \\
\mathtt{tausslow}\left( t \right) &= 72.4 + 3700 e^{ - \frac{1}{900} \left( 70 + \mathtt{V_m}\left( t \right) \right)^{2}} \\
\mathtt{Ito}\left( t \right) &= \left( \left( 1 - \mathtt{f\_is} \right) \mathtt{i\_sslow}\left( t \right) + \mathtt{f\_is} \mathtt{i\_s}\left( t \right) \right) \mathtt{gt} \left(
- \mathtt{E\_K}\left( t \right) + \mathtt{V_m}\left( t \right) \right) \mathtt{i\_r}\left( t \right) \\
\mathtt{IKs}\left( t \right) &= 2 \left( \mathtt{i\_nKs}\left( t \right) \right)^{2} \mathtt{GKs} \left( 0.68804 + 0.71283 \mathtt{IKUR\_PKAp} \right) \left(  - \mathtt{E\_K}\left( t \right) + \mathtt{V_m}\left( t \right) \right) \\
\mathtt{E\_Kr}\left( t \right) &= \frac{RT}{F} \ln\left( \frac{0.98 \mathtt{k\_o} + 0.02 \mathtt{na\_o}}{0.98 \mathtt{k\_i} + 0.02 \mathtt{na\_i}} \right) \\
\mathtt{IKr}\left( t \right) &= \mathtt{GKr} \left(  - \mathtt{E\_Kr}\left( t \right) + \mathtt{V_m}\left( t \right) \right) \mathtt{i\_OK}\left( t \right) \\
1 &= \mathtt{i\_IK}\left( t \right) + \mathtt{i\_CK2}\left( t \right) + \mathtt{i\_OK}\left( t \right) + \mathtt{CK0}\left( t \right) + \mathtt{i\_CK1}\left( t \right) \\
\mathtt{yinf}\left( t \right) &= \frac{1}{1 + e^{ - 0.15798 \left( -78.65 - \mathtt{V_m}\left( t \right) \right)}} \\
\mathtt{tauy}\left( t \right) &= \frac{1000}{0.56236 e^{0.070472 \left( -75 - \mathtt{V_m}\left( t \right) \right)} + 0.11885 e^{0.035249 \left( 75 + \mathtt{V_m}\left( t \right) \right)}} \\
\mathtt{IfNa}\left( t \right) &= \mathtt{fNa} \mathtt{gf} \left( \mathtt{V_m}\left( t \right) - \mathtt{E\_Na}\left( t \right) \right) \mathtt{i\_y}\left( t \right) \\
\mathtt{IfK}\left( t \right) &= \left( 1 - \mathtt{fNa} \right) \mathtt{gf} \left(  - \mathtt{E\_K}\left( t \right) + \mathtt{V_m}\left( t \right) \right) \mathtt{i\_y}\left( t \right) \\
\mathtt{If}\left( t \right) &= \mathtt{IfK}\left( t \right) + \mathtt{IfNa}\left( t \right)
\end{align}
$$

| Parameter    | Value    | Units | Description                                         |
| ------------ | -------- | ----- | --------------------------------------------------- |
| $G_{K1}$     | 0.0515   | mS/μF | Potassium channels conductance                      |
| $G_t$        | 0.1      | mS/μF | Transient outward potassium channels conductance    |
| $G_{Ks}$     | 0.05     | mS/μF | Potassium channels conductance                      |
| $\tau_{nKs}$ | 750      | ms    | Potassium channels time scale                       |
| $G_{Kr}$     | 0.06     | mS/μF | Potassium channels conductance                      |
| $k_{f, IKr}$ | 0.023761 | 1/ms  | Potassium channels transition rate                  |
| $k_{b, IKr}$ | 0.036778 | 1/ms  | Potassium channels transition rate                  |
| $G_f$        | 0.021    | mS/μF | Hyperpolarization activated current conductance     |
| $f_{Na}$     | 0.021    | -     | Hyperpolarization activated current sodium fraction |

### Calcium currents

L-type calcium channels, T-type calcium channels, and background calcium currents

$$
\begin{align}
\frac{\mathrm{d} \mathtt{i\_d}\left( t \right)}{\mathrm{d}t} &= \frac{\mathtt{dinf}\left( t \right) - \mathtt{i\_d}\left( t \right)}{\tau_d\left( t \right)} \\
\frac{\mathrm{d} \mathtt{i\_f}\left( t \right)}{\mathrm{d}t} &= \frac{\mathtt{finf}\left( t \right) - \mathtt{i\_f}\left( t \right)}{\tau_f\left( t \right)} \\
\frac{\mathrm{d} \mathtt{i\_fca}\left( t \right)}{\mathrm{d}t} &= \frac{\left( \mathtt{fcainf}\left( t \right) - \mathtt{i\_fca}\left( t \right) \right) \left( 1 - \left( \mathtt{fcainf}\left( t \right) > \mathtt{i\_fca}\left( t \right) \right) \left( \mathtt{Vm}\left( t \right) > -60 \right) \right)}{\tau_{fca}} \\
\frac{\mathrm{d} \mathtt{i\_b}\left( t \right)}{\mathrm{d}t} &= \frac{ - \mathtt{i\_b}\left( t \right) + \mathtt{binf}\left( t \right)}{\mathtt{taub}\left( t \right)} \\
\frac{\mathrm{d} \mathtt{i\_g}\left( t \right)}{\mathrm{d}t} &= \frac{\mathtt{ginf}\left( t \right) - \mathtt{i\_g}\left( t \right)}{\mathtt{taug}\left( t \right)} \\
\mathtt{ICa\_scale}\left( t \right) &= \mathtt{ICa\_scale0} \left( 1 + \frac{-0.56}{-1 + \frac{\mathtt{fracLCCbpISO}}{\mathtt{fracLCCbp0}}} \right) \\
\mathtt{E\_Ca}\left( t \right) &= 13.356 \ln\left( \frac{\mathtt{ca\_o}}{\mathtt{ca\_i}} \right) \\
\mathtt{INaCa}\left( t \right) &= \frac{\left( \mathtt{na\_i}^{3} \mathtt{ca\_o} e^{0.037436 \mathtt{gNaCa} \mathtt{Vm}\left( t \right)} - \mathtt{na\_o}^{3} \mathtt{ca\_i} \mathtt{fNaCa} e^{0.037436 \left( -1 + \mathtt{gNaCa} \right) \mathtt{Vm}\left( t \right)} \right) \mathtt{kNaCa} \mathtt{ICa\_scale}\left( t \right)}{1 + \left( \mathtt{na\_i}^{3} \mathtt{ca\_o} + \mathtt{na\_o}^{3} \mathtt{ca\_i} \mathtt{fNaCa} \right) \mathtt{dNaCa}} \\
\mathtt{ICaL}\left( t \right) &= \frac{14448 \mathtt{GCaL} \left(  - 0.341 \mathtt{ca\_o} + \mathtt{ca\_i} \left( 1 + expm1\left( 0.074872 \mathtt{Vm}\left( t \right) \right) \right) \right) \mathtt{ICa\_scale}\left( t \right) \mathtt{Vm}\left( t \right) \mathtt{i\_f}\left( t \right) \mathtt{i\_fca}\left( t \right) \mathtt{i\_d}\left( t \right)}{expm1\left( 0.074872 \mathtt{Vm}\left( t \right) \right)} \\
\mathtt{dinf}\left( t \right) &= \frac{1}{1 + e^{ - 0.13889 \left( 11.1 + \mathtt{Vm}\left( t \right) \right)}} \\
\tau_d\left( t \right) &= \frac{1.4 \left( 0.25 + \frac{1.4}{1 + e^{ - \frac{1}{13} \left( 35 + \mathtt{Vm}\left( t \right) \right)}} \right)}{1 + e^{ - \frac{1}{5} \left( -5
- \mathtt{Vm}\left( t \right) \right)}} + \frac{1}{1 + e^{ - \frac{1}{20} \left( -50 + \mathtt{Vm}\left( t \right) \right)}} \\
\mathtt{finf}\left( t \right) &= \frac{1}{1 + e^{ - 0.18519 \left( -23.3 - \mathtt{Vm}\left( t \right) \right)}} \\
\tau_f\left( t \right) &= 120 + \frac{165}{1 + e^{ - \frac{1}{10} \left( -25 + \mathtt{Vm}\left( t \right) \right)}} + 1125 e^{ - \frac{1}{240} \left( 27 + \mathtt{Vm}\left( t \right) \right)^{2}} \\
\mathtt{fcainf}\left( t \right) &= 0.68493 \left( 0.23 + \frac{0.00319}{0.00319 + \mathtt{ca\_i}^{8}} + \frac{0.1}{1 + e^{ - 10 \left( 0.5 - \mathtt{ca\_i} \right)}} + \frac{0.2}{1 + e^{ - 1.25 \left( 0.75 - \mathtt{ca\_i} \right)}} \right) \\
\mathtt{binf}\left( t \right) &= \frac{1}{1 + e^{ - 0.18497 \left( 37.491 + \mathtt{Vm}\left( t \right) \right)}} \\
\mathtt{taub}\left( t \right) &= 0.6 + \frac{5.4}{1 + e^{ - 0.03 \left( -100 - \mathtt{Vm}\left( t \right) \right)}} \\
\mathtt{ginf}\left( t \right) &= \frac{1}{1 + e^{ - \frac{1}{6} \left( -66 - \mathtt{Vm}\left( t \right) \right)}} \\
\mathtt{taug}\left( t \right) &= 1 + \frac{40}{1 + e^{ - 0.08 \left( -65 - \mathtt{Vm}\left( t \right) \right)}} \\
\mathtt{ICaT}\left( t \right) &= \mathtt{gCaT} \left( 106.5 - \mathtt{E\_Ca}\left( t \right) + \mathtt{Vm}\left( t \right) \right) \mathtt{i\_b}\left( t \right) \mathtt{i\_g}\left( t \right) \\
\mathtt{ICab}\left( t \right) &= \mathtt{gCab} \left(  - \mathtt{E\_Ca}\left( t \right) + \mathtt{Vm}\left( t \right) \right)
\end{align}
$$

| Parameter    | Value         | Units      | Description |
| ------------ | ------------- | ---------- | ----------- |
| fNaCa        | 1             | -          |             |
| kNaCa        | 2.268 * 10⁻¹⁶ | μAμF⁻¹μM⁻⁴ |             |
| dNaCa        | 10⁻¹⁶         | μM⁻⁴       |             |
| gNaCa        | 0.5           |            |             |
| $G_{CaL}$    | 6.3 * 10⁻⁵    | m³s⁻¹F⁻¹   |             |
| $\tau_{fca}$ | 10            | ms         |             |
| $g_{CaT}$    | 0.2           | mS/μF      |             |
| $g_{Cab}$    | 0.0008        | mS/μF      |             |

### Na-K pump

$$
\begin{align}
I_{NaK} &= I_{NaK}^{max} fNaK \frac{[K^+]_o}{ [K^+]_o + KmKo_{NaK} } \frac{[Na^+]_i^{nNaK}}{ [Na^+]_i^{nNaK} + KmNai_{NaK}^{nNaK} } \\
fNaK &= (1 + 0.1245 \exp(-0.1 V_m F / RT) + 0.0365 sigma \exp(V_m F / RT))^{-1} \\
sigma &= 1 / 7 \times (\exp([Na^+]_o / 67.3 \text{mM}) - 1)
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
