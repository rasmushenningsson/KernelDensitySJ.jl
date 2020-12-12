# KernelDensitySJ.jl

[![Build Status](https://travis-ci.com/rasmushenningsson/KernelDensitySJ.jl.svg?branch=master)](https://travis-ci.com/rasmushenningsson/KernelDensitySJ.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/rasmushenningsson/KernelDensitySJ.jl?svg=true)](https://ci.appveyor.com/project/rasmushenningsson/KernelDensitySJ-jl)
[![Coverage](https://codecov.io/gh/rasmushenningsson/KernelDensitySJ.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rasmushenningsson/KernelDensitySJ.jl)
[![Coverage](https://coveralls.io/repos/github/rasmushenningsson/KernelDensitySJ.jl/badge.svg?branch=master)](https://coveralls.io/github/rasmushenningsson/KernelDensitySJ.jl?branch=master)

## Overview

KernelDensitySJ.jl is a package for [Kernel Density Estimation](https://en.wikipedia.org/wiki/Kernel_density_estimation) and function estimation by [Gaussian Kernel Smoothing](https://en.wikipedia.org/wiki/Kernel_smoother#Gaussian_kernel_smoother).

More specifically, this package implements accurate, robust and reasonably fast algorithms for:

* Kernel Density Estimation in one dimension using Sheather and Jones algorithm for bandwidth selection, also known as width.SJ or bw.SJ. (The specific version of the algorithm is sometimes referred to as "solve-the-equation".) Only Gaussian kernels are supported.
* Evaluation of the density.
* Gaussian Kernel Smoothing, i.e. given a set of points (x and y values), computes a weighted average of the y values with weights coming from gaussian kernels placed at each x value.

Extra care is taken to handle outliers and extreme values for different parameters, ensuring that the results are accurate for a very wide range of inputs.

## Installation
Start the Julia REPL, press ']' to enter package mode and type
```
add KernelDensitySJ
```

## Usage Example
Let's start with some random points on the x-axis for a simple example.
```julia
using KernelDensitySJ
px = -1 .+ 2.0.*rand(40)
```
Now we can estimate the bandwidth
```julia
bw = bwsj(px)
```
and evaluate the density
```julia
x = range(-1.0,stop=1.0,length=201)
yd = density(px, bw, x)
```
![density](https://raw.githubusercontent.com/rasmushenningsson/KernelDensitySJ.jl/master/images/pl3.svg)


If we also have y values for the points
```julia
py = cos.(π.*px) .+ 0.2.*randn(40)
```
then we can make a smooth estimate of the underlying function by
```julia
ys = smooth(px, py, bw, x)
```
![smooth](https://raw.githubusercontent.com/rasmushenningsson/KernelDensitySJ.jl/master/images/pl4.svg)


The desired precision can be controlled by the `rtol` and `atol` keyword arguments.


Here is how it looks like with a couple of different values for the bandwidth:
![density2](https://raw.githubusercontent.com/rasmushenningsson/KernelDensitySJ.jl/master/images/pl1.svg)
![smooth2](https://raw.githubusercontent.com/rasmushenningsson/KernelDensitySJ.jl/master/images/pl2.svg)


## Technical Details
The bandwidth is computed by solving equation (12) in [Sheather and Jones, 1991] using root-finding methods.
To solve the equation both accurately and fast, this packages computes and successively refines lower and upper bounds of the target function for a given bandwidth until the bounds no longer includes 0, which is enough to bracket the root.

The `density` and `smooth` functions use a similar approach by refining upper and lower bounds until the desired tolerance has been achieved. It's common that the weight is very close to zero for points far away from the evaluation point and this strategy ensures that we do not spend unnecessary computing power there, while still having accurate bounds for the error.

The upper and lower bounds are computed by splitting the estimated functions into intervals where the functions are convex/concave and using [Jensen's inequality](https://en.wikipedia.org/wiki/Jensen's_inequality) for one direction and a linear approximation for the other.



## References
* Sheather, S. J., & Jones, M. C. (1991). A reliable data‐based bandwidth selection method for kernel density estimation. Journal of the Royal Statistical Society: Series B (Methodological), 53(3), 683-690.


# Other Kernel Density Estimation Packages
The following packages also implement Kernel Density Estimation in Julia:

* [KernelDensity.jl](https://github.com/JuliaStats/KernelDensity.jl)
* [KernelDensityEstimate.jl](https://github.com/JuliaRobotics/KernelDensityEstimate.jl)
