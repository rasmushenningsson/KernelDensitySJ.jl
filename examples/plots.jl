using KernelDensitySJ
using Gadfly
using Colors
using DataFrames
using StableRNGs

theme = Theme(background_color=colorant"white")

rng = StableRNG(123456)

px = -1 .+ 2.0.*rand(rng, 40)
py = cos.(π.*px) .+ 0.2.*randn(rng, length(px))

x = range(-1.0,stop=1.0,length=201)
bandwidths = [0.03, 0.1, 0.3]
yd = density(px, bandwidths', x)
yd ./= maximum(yd;dims=1)
ys = smooth(px, py, bandwidths', x)

df = DataFrame(x=repeat(x,3), bandwidth=repeat(string.(bandwidths),inner=length(x)), density=vec(yd), smoothed=vec(ys))

pl1 = plot(layer(df, x=:x, y=:density, color=:bandwidth, Geom.line),
           layer(x=px,y=0*py,color=[colorant"black"],Geom.point), theme)
pl2 = plot(layer(df, x=:x, y=:smoothed, color=:bandwidth, Geom.line),
           layer(x=px,y=py,color=[colorant"black"],Geom.point), theme)

bwopt = bwsj(px)
ydopt = density(px, bwopt, x)
ydopt ./= maximum(ydopt)
ysopt = smooth(px, py, bwopt, x)

dfopt = DataFrame(x=x, bandwidth=string(round(bwopt,digits=3)), density=ydopt, smoothed=ysopt)

pl3 = plot(layer(dfopt, x=:x, y=:density, color=:bandwidth, Geom.line),
           layer(x=px,y=0*py,color=[colorant"black"],Geom.point), theme)
pl4 = plot(layer(dfopt, x=:x, y=:smoothed, color=:bandwidth, Geom.line),
           layer(x=px,y=py,color=[colorant"black"],Geom.point), theme)

draw(SVG("pl1.svg", 16cm, 12cm), pl1)
draw(SVG("pl2.svg", 16cm, 12cm), pl2)
draw(SVG("pl3.svg", 16cm, 12cm), pl3)
draw(SVG("pl4.svg", 16cm, 12cm), pl4)
