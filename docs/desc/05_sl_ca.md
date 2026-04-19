# Sarcolemmal Calcium channels

L-type calcium channels, T-type calcium channels, and background calcium currents.

LCC

$$
\begin{align}
J_{CaSL} &= (2 I_{NaCa} - I_{CaL} - I_{CaT} - I_{Cab}) \frac{A_{CAP} C_m}{2 F V_{subSL}} \\
\mathrm{I_{CaL}}  &= \mathrm{ICa_{scale}} \cdot G_{CaL,k} \cdot i_d \cdot i_f \cdot i_fca \cdot GHK(G_{CaL}, 2, V_m, \mathrm{ca_{sl}}, 0.341 \mathrm{ca_o}) \\
\cdot G_{CaL,k} &= G_{CaL} (1 + 0.1CaMK_{Act}) \\
\frac{d i_d }{dt} &= \frac{d_∞  - i_d}{τ_d} \\
\frac{d i_f }{dt} &= \frac{f_∞  - i_f}{τ_f} \\
\frac{d i_{fca} }{dt} &= \frac{(fca_∞ - i_{fca}) \left( 1 - \left( fca_∞  > i_{fca}  \right) \left( V_m  > -60 \text{mV} \right) \right)}{τ_{fca}} \\
\mathrm{ICa_{scale}}  &= \mathrm{ICa_{scale, 0}} \left( 1 + \frac{0.56}{1 - \frac{\mathrm{fracLCCbpISO}}{\mathrm{fracLCCbp0}}} \right) \\
d_∞ &= (1 + \exp(-(V + 11.1) / 7.2))^{-1} \\
τ_d &= (\alpha_d \beta_d + \gamma_d) \\
\alpha_d &= 1.4 (1 + \exp(-(V_m + 35) / 13))^{-1} + 0.25 \\
\beta_d &= 1.4 (1 + \exp((V_m + 5) / 5))^{-1} \\
\gamma_d &= (1 + \exp(-(V_m - 50) / 20))^{-1} \\
f_∞ &= (1 + \exp((V_m + 23.3) / 5.4))^{-1} \\
τ_f &= 120 + \frac{165}{1 + \exp(-(V_m - 25) / 10)} + 1125 \exp( -(V_m + 27)^{2} / 240) \\
fca_∞  &= (\alpha_{fca} + \beta_{fca} + \gamma_{fca} + 0.23) / 1.46 \\
\alpha_{fca} &= \frac{0.4875^8}{0.4875 ^8 + \mathrm{ca_{sl}}^8} \\
\beta_{fca} &= 0.1 {1 + \exp((\mathrm{ca_{sl}} - 0.5) / 0.1)}^{-1} \\
\gamma_{fca} &= 0.2 {1 + \exp((\mathrm{ca_{sl}} - 0.75) / 0.8)}^{-1} \\
\end{align}
$$

TCC

$$
\begin{align}
\mathrm{I_{CaT}} &= \mathrm{gCaT} \cdot i_b \cdot i_g ( V_m + 106.5 - \mathrm{E_Ca}) \\
E_{Ca} &= 0.5V_T \ln \frac{ca_o}{ca_{sl}} \\
\frac{d}{dt} i_b &= \frac{b_∞ - i_b}{τ_b} \\
\frac{d}{dt} i_g &= \frac{g_∞ - i_g}{τ_g} \\
b_∞ &= (1 + \exp(-(V_m + 37.49098) / 5.40634))^{-1} \\
τ_b &= 0.6 + 5.4 (1 + \exp(0.03 (V_m + 100)))^{-1} \\
g_∞ &= (1 + \exp((V_m + 66) / 6))^{-1} \\
τ_g &= 1 + 40 (1 + \exp(0.08 (V_m + 65)))^{-1} \\
\end{align}
$$

Background

$$
\begin{align}
\mathrm{I_{Cab}} &= \mathrm{gCab} (V_m - \mathrm{E_Ca}) \\
\end{align}
$$

NCX
$$
\begin{align}
\mathrm{I_{NaCa}}  &= \mathrm{kNaCa} \cdot \mathrm{ICa_{scale}} \frac{\mathrm{na_i}^{3} \mathrm{ca_o} \exp( \mathrm{gNaCa} V_m /V_T ) - \mathrm{na_o}^{3} \mathrm{ca_{sl}} \mathrm{fNaCa} \exp (( \mathrm{gNaCa} - 1 ) V_m F /RT )} {1 + \left(\mathrm{na_i}^3 \mathrm{ca_o} + \mathrm{na_o}^{3} \mathrm{ca_{sl}} \mathrm{fNaCa} \right) \mathrm{dNaCa}} \\
\end{align}
$$
