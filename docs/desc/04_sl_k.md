# Sarcolemmal Potassium channels

IK1:

$$
\begin{align}
\mathrm{I_{K1}}  &= \mathrm{G_{K1}} \cdot \frac{k_o}{k_o + 210 \mathrm{\mu M}}  \frac{ ( V_m - 6.1373 - E_K )}{0.1653 + \exp(0.0319 (V_m - 6.1373 - E_K))} \\
E_K &= V_T \ln\frac{k_o}{k_i}  \\
\end{align}
$$

Ito:

$$
\begin{align}
\mathrm{I_{to}} &= G_t \cdot i_r ((1 - f_{is} ) i_{s,slow} + f_{is} i_s )(V_m - E_K) \\
\frac{d}{dt} i_r &= \frac{r_∞  - i_r }{τ_r} \\
\frac{d}{dt} i_s &= \frac{s_∞  - i_s }{τ_s} \\
\frac{d}{dt} i_{s,slow} &= \frac{slow_∞ - i_{s,slow} }{τ_{s,slow}} \\
s_∞ &= (1 + \exp((V_m + 31.97156) / 4.64291))^{-1} \\
r_∞ &= (1 + \exp(-(V_m - 3.55716) / 14.61299))^{-1 } \\
slow_∞ &= s_∞ \\
τ_r &= \frac{1000 \text{ms}}{45.16 \exp(0.03577( V_m + 50)) + 98.9 \exp( -0.1 ( V_m + 38))}
\\
τ_s &= 8.1\text{ms} + 350\text{ms} \cdot \exp( - \frac{1}{225}(V_m + 70)^{2}) \\
τ_{s,slow} &= 72.4\text{ms} + 3700\text{ms} \cdot \exp( - \frac{1}{900}( V_m + 70 )^{2}) \\
\end{align}
$$

IKs:

$$
\begin{align}
\mathrm{I_{Ks}}  &= 2 G_{Ks} \cdot i_{nKs}^{2}(0.68804 + 0.71283 \mathrm{IKUR_{PKAp}}) (V_m - E_K) \\
\frac{d}{dt} i_{nKs} &= \frac{nks_∞ - i_{nKs}}{τ_{nKs}} \\
nks_∞ &= \frac{\alpha_{nks}}{\alpha_{nks} + \beta_{nks}} \\
\alpha_{nks} &= 0.00000481333 / 0.128 \cdot exprel(-0.128 (V_m + 26.5)) \\
\beta_{nks} &= 0.0000953333 \cdot \exp(-0.038 (V_m + 26.5)) \\
\end{align}
$$

IKr:

$$
\begin{align}
\mathrm{I_{Kr}} &= G_{Kr} \cdot i_{OK} (V_m - E_{Kr}) \\
E_{Kr} &= V_T \ln\left( \frac{0.98 k_o + 0.02na_o}{0.98 k_i + 0.02a_i} \right) \\
\frac{d i_{CK1} }{dt} &= k_{b, IKr} i_{CK2} - k_{f, IKr} i_{CK1} + k_{01} i_{CK0} - k_{10} i_{CK1}  \\
\frac{d i_{CK2} }{dt} &= - k_{b, IKr} i_{CK2} + k_{f, IKr} i_{CK1} - k_{2o} i_{CK2} + k_{o2} i_{OK}  \\
\frac{d i_{OK} }{dt} &= k_{io} i_{IK} + k_{2o} i_{CK2} - k_{o2} i_{OK} - k_{oi} i_{OK}  \\
\frac{d i_{IK} }{dt} &= -k_{io} i_{IK} + k_{oi} i_{OK}  \\
k_{01} &= 0.022348\mathrm{ms}^{-1} \exp(0.01176 V_m ) \\
k_{10} &= 0.047002\mathrm{ms}^{-1} \exp(-0.0631 V_m ) \\
k_{2o} &= 0.013733\mathrm{ms}^{-1} \exp(0.038198 V_m) \\
k_{o2} &= 6.89 \cdot 10^{-5}\mathrm{ms}^{-1} \exp(-0.04178 V_m) \\
k_{oi} &= 0.090821\mathrm{ms}^{-1} \exp(0.023391 V_m) \\
k_{io} &= 0.006497\mathrm{ms}^{-1} \exp(-0.03268 V_m) \\
1 &= i_{CK0} + i_{CK1} + i_{CK2} + i_{OK} + i_{IK} \\
\end{align}
$$


If:

$$
\begin{align}
\mathrm{I_{fNa}} &= f_{Na} \cdot G_f \cdot i_y \cdot ( V_m - E_{Na}) \\
\mathrm{I_{fK}}  &= (1 - f_{Na}) \cdot G_f \cdot i_y \cdot ( V_m - E_K ) \\
\mathrm{I_{f}}  &= \mathrm{I_{fK}} + \mathrm{I_{fNa}} \\
\frac{d}{dt} i_y &= \frac{y_∞ - i_y}{τ_y} \\
y_∞ &= expit(-0.15798 (V_m + 78.65)) \\
τ_y  &= \frac{1000 \text{ms}}{0.56236 \exp(-0.070472 (V_m + 75)) + 0.11885 \exp(0.035249 ( V_m + 75))} \\
\end{align}
$$
