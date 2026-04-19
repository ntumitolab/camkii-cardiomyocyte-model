# Sarcolemmal Na-K pump

$$
\begin{align}
\mathrm{I_{NaK}} &= I_{NaK}^{max} fNaK \frac{\mathrm{k_o}}{ \mathrm{k_o} + KmKo_{NaK} } \frac{\mathrm{na_i}^{nNaK}}{ \mathrm{na_i}^{nNaK} + KmNai_{NaK}^{nNaK} } \\
fNaK &= (1 + 0.1245 \exp(-0.1 V_m/ V_T) + 0.0365 \sigma_{NaK} \exp(V_m/V_T))^{-1} \\
\sigma_{NaK} &= (\exp(\mathrm{na_i} / 67.3 \text{mM}) - 1) / 7
\end{align}
$$

ODEs for sarcolemmal ions

$$
\begin{align}
C_m \frac{d}{dt} \mathrm{V_m} &= -(I_{Nab} + I_{NaCa} + I_{CaL} + I_{CaT} + I_{f} + I_{to} + I_{K1} + I_{Ks} + I_{Kr} + I_{Na} + I_{NaK} + I_{Cab} + I_{stim}) \\
\frac{d}{dt} \mathrm{na_i} &= -(I_{fNa} + I_{Nab} + I_{Na} + 3 I_{NaCa} + 3 I_{NaK}) \frac{A_{cap} C_m}{ F V_{myo}} \\
\frac{d}{dt} \mathrm{k_i} &= -(I_{fK} + I_{to} + I_{K1} + I_{Ks} + I_{Kr} + I_{stim} - 2 I_{NaK}) \frac{A_{cap} C_m}{ F V_{myo}}  \\
\end{align}
$$
