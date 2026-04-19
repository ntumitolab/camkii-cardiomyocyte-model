# General functions

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

GHK flux equation.

$$
GHK(p, z, V, S_i, S_o) := pzF\delta\frac{\exp(\delta)S_i - S_o}{\exp(\delta) - 1} ; \quad  \delta = \frac{zVF}{RT}
$$

Ion concentrations:

$$
\begin{align}
na_i &= \mathrm{[Na^+]_i} \\
k_i &= \mathrm{[K^+]_i} \\
ca_i &= \mathrm{[Ca^{2+}]_i} \\
na_o &= \mathrm{[Na^+]_o} \\
k_o &= \mathrm{[K^+]_o} \\
ca_o &= \mathrm{[Ca^{2+}]_o} \\
\end{align}
$$

Reversal potentials:

$$
\begin{align}
E_{Na} &= V_T \ln \frac{na_o}{na_i} \\
E_{K } &= V_T \ln\frac{k_o}{k_i}  \\
E_{Ca} &= 0.5V_T \ln \frac{ca_o}{ca_{sl}} \\
E_{Kr} &= V_T \ln\left( \frac{0.98 k_o + 0.02na_o}{0.98 k_i + 0.02a_i} \right) \\
\end{align}
$$
