using KernelDensitySJ
using Gadfly
using Colors
using StableRNGs

set_default_plot_size(30cm, 30cm)

rng = StableRNG(123456)

px = -1 .+ 2.0.*rand(rng, 60)
py = cos.(π.*px) .+ 0.2.*randn(rng, length(px))

x = range(-1.0,stop=1.0,length=1001)
bandwidths = [0.3, 0.1, 0.03]
yd = density(px, bandwidths', x)
yd ./= maximum(yd;dims=1)

colors = [colorant"blue",colorant"orange",colorant"purple"]

layers = [layer(x=x,y=yd[:,i],color=[colors[i]],Geom.line) for i=1:length(bandwidths)]
pl1 = plot(layer(x=px,y=zeros(size(px)),color=[colorant"black"],Geom.point), layers...)

ys = smooth(px, py, bandwidths', x)
layers = [layer(x=x,y=ys[:,i],color=[colors[i]],Geom.line) for i=1:length(bandwidths)]
pl2 = plot(layer(x=px,y=py,color=[colorant"black"],Geom.point), layers...)



bwopt = bwsj(px)
@show bwopt
ydopt = density(px, bwopt, x)
ydopt ./= maximum(ydopt)
pl3 = plot(layer(x=px,y=zeros(size(px)),color=[colorant"black"],Geom.point), layer(x=x,y=ydopt,Geom.line))

ysopt = smooth(px, py, bwopt, x)
pl4 = plot(layer(x=px,y=py,color=[colorant"black"],Geom.point), layer(x=x,y=ysopt,Geom.line))

gridstack([pl1 pl2; pl3 pl4])
