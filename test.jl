using ProgressLogging
using ModelingToolkit
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using LsqFit
using CaMKIIModel
using CaMKIIModel: μM
Plots.default(lw=2)

@parameters Ca=0μM ROS=0μM
sys = get_camkii_sys(Ca; ROS, phospho_rate=0, simplify=true)

tend=100.0
prob = ODEProblem(sys, [], tend)
alg = Rodas5()
callback = TerminateSteadyState()
ca = exp10.(range(log10(1e-2μM), log10(10μM), length=1001))
prob_func = (prob, i, repeat) -> begin
    remake(prob, p=[Ca => ca[i]])
end
trajectories = length(ca)
sol0 = solve(remake(prob, p=[Ca => 0.0]), alg; callback) ## warmup
sim = solve(EnsembleProblem(prob; prob_func, safetycopy=false), alg; trajectories, callback)

extract(sim, k) = map(s->s[k][end], sim)
hil_s(x, k, n) = sign(x*k) * (abs(x)^n) / (abs(x)^n + abs(k)^n)
sqrt_s(x) = sign(x) * sqrt(abs(x))

sol0[sys.CaMKAct][end]
plot(ca .* 1000, extract(sim, sys.CaMKAct), lab="Active CaMKII", xlabel="Ca (μM)", xscale=:log10)

xdata = ca
ydata = extract(sim, sys.CaMKAct)

function model(x, p)
    A = 30μM
    B = 70μM
    kmca = p[1]
    nca = p[2]
    kfb = p[3]
    k0 = p[4]
    y = map(x) do ca
        keq = kfb * hil_s(ca, kmca, nca) + k0
        xterm = A + B + 1/keq
        camkcam4 = 0.5 * ( xterm - sqrt_s(xterm^2 - 4 * A * B))
        y = camkcam4 / B
    end
end
