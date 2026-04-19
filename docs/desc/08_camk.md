# CaMKII system
$$
\begin{align}
CaMK_{act} &= 1 - CaMK \\
1 &=  CaMK + CaMKB + CaMKBOX + CaMKP + CaMKPOX + CaMKA + CaMKA2 + CaMKAOX + CaMKOX \\
\frac{d}{dt} CaMKB &= -v_{IB} - v_{BP} - v_{BBo} \\
\frac{d}{dt} CaMKP &= v_{AP} + v_{BP} - v_{PPo} \\
\frac{d}{dt} CaMKA &= -v_{AP} - v_{AI} - v_{A1A2} + v_{AoA} \\
\frac{d}{dt} CaMKA2 &= v_{A1A2} \\
\frac{d}{dt} CaMKBOX &= v_{IoBo} - v_{BoPo} + v_{BBo} \\
\frac{d}{dt} CaMKPOX &= v_{AoPo} + v_{BoPo} + v_{PPo} \\
\frac{d}{dt} CaMKAOX &= -v_{AoPo} - v_{AoA} - v_{AoIo} \\
\frac{d}{dt} CaMKOX &= -v_{IoBo} + v_{AoIo} - v_{IoI} \\
v_{IB} &= k_f \cdot CaMK - k_b \cdot CaMKB \\
v_{IoBo} &= k_f \cdot r_{CaMKO} \cdot CaMKOX - k_b \cdot CaMKBOX \\
v_{AP} &= k_f \cdot r_{CaMKP} \cdot CaMKA - k_b \cdot CaMKP \\
v_{AoPo} &= k_f \cdot r_{CaMKP} \cdot CaMKAOX - k_b \cdot CaMKPOX \\
camkb_\infty &=  kfa2_{CaMK} \frac{ca_{avg}^2}{ca_{avg}^2 + kmCa2_{CaMK}^2} + kfa4_{CaMK} \frac{ca_{avg}^4}{ca_{avg}^4 + kmCa4_{CaMK}^4} + kfb_{CaMK} \\
k_f &= v_{CaMK} \cdot camkb_\infty \\
k_b &= v_{CaMK} \cdot (1 - camkb_\infty) \\
kph &= kphos_{CaMK} \cdot aMK_{act}  \\
v_{BP} &= kph \cdot CaMKB - kdeph_{CaMK} \cdot CaMKP \\
v_{BoPo} &= kph \cdot CaMKBOX - kdeph_{CaMK} \cdot CaMKPOX \\
v_{A1A2} &= k_{P1P2} \cdot CaMKA - k_{P2P1} \cdot CaMKA2 \\
v_{AI} &= kdeph_{CaMK} \cdot CaMKA \\
v_{AoIo} &= kdeph_{CaMK} \cdot CaMKAOX \\
v_{BBo} &= kox_{CaMK} \cdot \mathtt{H_2O_2} \cdot CaMKB - krd_{CaMK} \cdot CaMKBOX \\
v_{PPo} &= kox_{CaMK} \cdot \mathtt{H_2O_2} \cdot CaMKP - krd_{CaMK} \cdot CaMKPOX \\
v_{IoI} &= krd_{CaMK} \cdot CaMKOX \\
v_{AoA} &= krd_{CaMK} \cdot CaMKAOX \\
\end{align}
$$
