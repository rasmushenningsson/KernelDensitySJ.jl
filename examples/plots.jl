using KernelDensitySJ
using Gadfly
using Colors

px = rand(10)
py = Float64.(px.>0.5)

x = range(-0.5,stop=1.5,length=201)
bandwidths = [0.02, 0.06, 0.2]
yd = density(px, bandwidths', x)
yd ./= maximum(yd;dims=1)

colors = [colorant"blue",colorant"orange",colorant"purple"]

layers = [layer(x=x,y=yd[:,i],color=[colors[i]],Geom.line) for i=1:length(bandwidths)]
display(plot(layer(x=px,y=zeros(size(px)),color=[colorant"black"],Geom.point), layers...))

ys = smooth(px, py, bandwidths', x)
layers = [layer(x=x,y=ys[:,i],color=[colors[i]],Geom.line) for i=1:length(bandwidths)]
display(plot(layer(x=px,y=py,color=[colorant"black"],Geom.point), layers...))



bwopt = bwsj(px)
ydopt = density(px, bwopt, x)
ydopt ./= maximum(ydopt)
display(plot(layer(x=px,y=zeros(size(px)),color=[colorant"black"],Geom.point), layer(x=x,y=ydopt,Geom.line)))

ysopt = smooth(px, py, bwopt, x)
display(plot(layer(x=px,y=py,color=[colorant"black"],Geom.point), layer(x=x,y=ysopt,Geom.line)))
