# Sarcolemmal sodium channels

$$
\begin{align}
\mathrm{I_{Na}} &= \bar G_{Na} \cdot m_{Na}^{3} \cdot h_{Na} \cdot j_{Na} (V_m-E_{Na}) \\
\mathrm{I_{Na,b}} &= \bar G_{Na,b} (V_m - E_{Na}) \\
E_{Na} &= V_T \ln \frac{na_o}{na_i} \\
\frac{dm_{Na}}{dt} &= \alpha_{m} - m_{Na}(\alpha_{m} + \beta_{m}) \\
\frac{dh_{Na}}{dt} &= \alpha_{h} - h_{Na}(\alpha_{h} + \beta_{h}) \\
\frac{dj_{Na}}{dt} &= \alpha_{j} - m_{Na}(\alpha_{j} + \beta_{j}) \\
\alpha_{m} &= 3.2 \mathrm{ms}^{-1} exprel(-(V_m + 47.13) / 10)  \\
\beta_{m} &= 0.08 \mathrm{ms}^{-1} \exp(-V_m / 11) \\
\alpha_{h} &= 0.135 \text{ms}^{-1} \exp(-(V_m+80)/6.8) \\
\beta_{h} &= \frac{7.6923}{ 1 + \exp(-(V_m+10.66)/11.1) } \mathrm{ms}^{-1} \\
\alpha_{j} &= (-127140 \exp(0.2444 V_m)-3.474 \cdot 10^{-5}\exp(-0.04391 V_m))\frac{V_m + 37.78}{1 + \exp(0.311( V_m + 79.23))} \mathrm{ms}^{-1} \\
\beta_{j} &= 0.3 \frac{\exp(-2.535 \cdot 10^{-7}V_m)}{1 + \exp(-0.1(V_m + 32))} \mathrm{ms}^{-1} \\
\end{align}
$$
