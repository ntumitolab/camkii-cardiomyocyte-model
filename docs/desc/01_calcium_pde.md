## Cytosolic calcium diffusion

Cytosolic calcium ions freely diffuse between sub-sarcolemma (SL) and sub-sarcoplasmic (SR) spaces.

Calcium buffering in each compartment:

$$
\begin{align}
\beta_{Ca,x} &= \left(1 + \frac{\Sigma Trpn \cdot Km_{Trpn, 2}}{(ca_x + Km_{Trpn, 2})^2}  + \frac{\Sigma Cmdn \cdot Km_{Cmdn}}{(ca_x + Km_{Cmdn, 2})^2} \right)^{-1} \\
Km_{Trpn, 2} &= \frac{Km_{Trpn}}{fPKA_{TnI}}  \\
fPKA_{TnI} &= 1.61 - 0.61 \frac{1 - TnI_{PKAp}}{1 - fracTnIp_0}
\end{align}
$$

Calcium diffusion space is divided into $(r_{SL} - r_{SR}) / dx$ concentric compartments.

For i = 2 to $(r_{SL} - r_{SR}) / dx - 1$

$$
\begin{align}
\frac{d}{dt}ca_i &= \frac{D_{ca} \cdot \beta_{Ca,i}}{dx^2 \cdot j_i} ((j_i + 1)ca_{i+1} - 2j_i \cdot ca_{i} + (j_i - 1)ca_{i-1}) \\
j_i &= r_{SR} / dx + i - 1
\end{align}
$$

For the boundaries:

$$
\begin{align}
\frac{d}{dt}ca_1 &= \frac{D_{ca} \cdot \beta_{Ca,1}}{dx^2 \cdot j_1} ((j_1 + 1)ca_{2} - 2j_1 \cdot ca_{1} + (j_1 - 1)ca_{1} + J_{CaSR})  \\
\frac{d}{dt}ca_n &= \frac{D_{ca} \cdot \beta_{Ca,n}}{dx^2 \cdot j_n} ((j_n + 1)ca_{n} - 2j_n \cdot ca_{n} + (j_n - 1)ca_{n-1} + J_{CaSL})  \\
j_1 &= r_{SR} / dx \\
j_n &= r_{SL} / dx \\
ca_{sr} &= ca_1 \\
ca_{sl} &= ca_n \\
\end{align}
$$

Average calcium concentration:

$$
ca_{avg} = \frac{\sum^N_{i=1} ca_i}{(r_{SL} - r_{SR}) / dx}
$$
