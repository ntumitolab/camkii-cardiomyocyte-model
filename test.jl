using Plots

f1(x) = 0.2205 + 0.2339 * (x^0.9932) / (x^0.9932 + 0.00726^0.9932)
f2(x) = 0.2205 + 0.2339 * (x) / (x + 0.00726)

xx = logrange(1e-4, 3, length=1001)
y1 = f1.(xx)
y2 = f2.(xx)

plot(xx, [y1 y2], lab=["full" "reduced"], xscale=:log10)

plot(xx, (y2 .- y1) ./ y1 .* 100, lab="Rel err", xscale=:log10)
