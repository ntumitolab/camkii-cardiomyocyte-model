# Original CaMKII model

Reference model by Chang and coworkers for fitting CaMKII-CaMCa binding.

$$
\begin{align}
ca &= \mathrm{[Ca]^{2+}_{avg}} \\
v_1 &= k_1 \cdot ca^2 \cdot CaM_{Apo} \\
k_1 &= \frac{CaM1C_{on} CaM2C_{on}}{CaM1C_{off} + CaM2C_{on} \cdot ca} \\
v_2 &= k_2 \cdot Ca2CaM_C \\
k_2 &= \frac{CaM1C_{off} CaM2C_{off}}{CaM1C_{off} + CaM2C_{on} \cdot ca} \\
v_3 &= k_3 \cdot ca^2 \cdot CaM_{Apo} \\
k_3 &= \frac{CaM1N_{on} CaM2N_{on}}{CaM1N_{off} + CaM2N_{on} \cdot ca} \\
v_4 &= k_4 \cdot Ca2CaM_N \\
k_4 &= \frac{CaM1N_{off} CaM2N_{off}}{CaM1N_{off} + CaM2N_{on} \cdot ca} \\
v_5 &= k_5 \cdot ca^2 \cdot Ca2CaM_C \\
k_5 &= k_3 \\
v_6 &= k_6 \cdot Ca4CaM \\
k_6 &= k_4 \\
v_7 &= k_7 \cdot ca^2 \cdot Ca2CaM_N \\
k_7 &= k_1 \\
v_8 &= k_8 \cdot Ca4CaM \\
k_8 &= k_2 \\
v_9 &= k_9 \cdot ca^2 \cdot KCaM_{Apo} \\
k_9 &= \frac{KCaM1C_{on} KCaM2C_{on}}{KCaM1C_{off} + KCaM2C_{on} \cdot ca} \\
v_{10} &= k_{10} \cdot KCa2CaM_C \\
k_{10} &= \frac{KCaM1C_{off} KCaM2C_{off}}{KCaM1C_{off} + KCaM2C_{on} \cdot ca} \\
v_{11} &= k_{11} \cdot ca^2 \cdot KCaM_{Apo} \\
k_{11} &= \frac{KCaM1N_{on} KCaM2N_{on}}{KCaM1N_{off} + KCaM2N_{on} \cdot ca} \\
v_{12} &= k_{12} \cdot KCa2CaM_N \\
k_{12} &= \frac{KCaM1N_{off} KCaM2N_{off}}{KCaM1N_{off} + KCaM2N_{on} \cdot ca} \\
v_{13} &= k_{13} \cdot ca^2 \cdot KCa2CaM_C \\
k_{13} &= k_{11} \\
v_{14} &= k_{14} \cdot KCa4CaM \\
k_{14} &= k_{12} \\
v_{15} &= k_{15} \cdot ca^2 \cdot KCa2CaM_N \\
k_{15} &= k_{9} \\
v_{16} &= k_{16} \cdot KCa4CaM \\
k_{16} &= k_{10} \\
v_{17} &= k_{17} \cdot CaMK \cdot CaM_{Apo} \\
v_{18} &= k_{18} KCaM_{Apo} \\
v_{19} &= k_{19} \cdot CaMK \cdot Ca2CaM_C \\
v_{20} &= k_{20} KCa2CaM_C \\
v_{21} &= k_{21} \cdot CaMK \cdot Ca2CaM_N \\
v_{22} &= k_{22} KCa2CaM_N \\
v_{23} &= k_{23} \cdot CaMK \cdot Ca4CaM \\
v_{24} &= k_{24} KCa4CaM \\
\end{align}
$$

$$
\begin{align}
\frac{d}{dt} Ca2CaM_C &= v_1 - v_2 - v_5 + v_6 - v_{19} + v_{20} \\
\frac{d}{dt} Ca2CaM_N &= v_3 - v_4 - v_7 + v_8 - v_{21} + v_{22} \\
\frac{d}{dt} Ca4CaM &= v_5 - v_6 + v_7 - v_8 - v_{23} + v_{24} \\
\frac{d}{dt} KCaM_{Apo} &= -v_9 + v_{10} - v_{11} + v_{12} + v_{17} - v_{18} \\
\frac{d}{dt} KCa2CaM_C &= v_9 - v_{10} - v_{13} + v_{14} + v_{19} - v_{20} \\
\frac{d}{dt} KCa2CaM_N &= v_{11} - v_{12} - v_{15} + v_{16} + v_{21} - v_{22} \\
\frac{d}{dt} KCa4CaM &= v_{13} - v_{14} + v_{15} - v_{16} + v_{23} - v_{24}\\
\end{align}
$$

$$
\begin{align}
CaM_T &= CaM_{Apo} + Ca2CaM_C + Ca2CaM_N + Ca4CaM + KCaM_{Apo} + KCa2CaM_C + KCa2CaM_N + KCa4CaM \\
CaMK_T &= CaMK + KCaM_{Apo} + KCa2CaM_C + KCa2CaM_N + KCa4CaM
\end{align}
$$

## Parameters

k17: Chang's parameter (3.8) is different from Pepke's (0.0038), probably a unit error.

| Parameter  | Value  | Units  | Description |
| ---------- | ------ | ------ | ----------- |
| CaM_T      | 30     | μM     |             |
| CaMK_T     | 70     | μM     |             |
| CaM1C_on   | 5      | 1/μM/s |             |
| CaM1C_off  | 50     | 1/s    |             |
| CaM2C_on   | 10     | 1/μM/s |             |
| CaM2C_off  | 10     | 1/s    |             |
| CaM1N_on   | 100    | 1/μM/s |             |
| CaM1N_off  | 2000   | 1/s    |             |
| CaM2N_on   | 200    | 1/μM/s |             |
| CaM2N_off  | 500    | 1/s    |             |
| KCaM1C_on  | 44     | 1/μM/s |             |
| KCaM1C_off | 33     | 1/s    |             |
| KCaM2C_on  | 44     | 1/μM/s |             |
| KCaM2C_off | 0.8    | 1/s    |             |
| KCaM1N_on  | 76     | 1/μM/s |             |
| KCaM1N_off | 300    | 1/s    |             |
| KCaM2N_on  | 76     | 1/μM/s |             |
| KCaM2N_off | 20     | 1/s    |             |
| k_17       | 0.0038 | 1/μM/s |             |
| k_18       | 5.5    | 1/s    |             |
| k_19       | 0.92   | 1/μM/s |             |
| k_20       | 6.8    | 1/s    |             |
| k_21       | 0.12   | 1/μM/s |             |
| k_22       | 1.7    | 1/s    |             |
| k_23       | 30     | 1/μM/s |             |
| k_24       | 1.5    | 1/s    |             |
