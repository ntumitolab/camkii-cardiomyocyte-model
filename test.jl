using Plots

exprel(x) = x / expm1(x)

exprel2(x) = log2(1 + exp2(-x))

xx = -3:0.01:+6
y1 = exprel.(xx)
y2 = exprel2.(xx)

plot(xx, [y1 y2], lab=["exprel ln"])

plot(xx, (y2 .- y1) ./ y1 .* 100, lab="Rel err")
