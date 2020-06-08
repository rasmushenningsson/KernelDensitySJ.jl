# KernelDensitySJ.jl

[![Build Status](https://travis-ci.com/rasmushenningsson/KernelDensitySJ.jl.svg?branch=master)](https://travis-ci.com/rasmushenningsson/KernelDensitySJ.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/rasmushenningsson/KernelDensitySJ.jl?svg=true)](https://ci.appveyor.com/project/rasmushenningsson/KernelDensitySJ-jl)
[![Coverage](https://codecov.io/gh/rasmushenningsson/KernelDensitySJ.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rasmushenningsson/KernelDensitySJ.jl)
[![Coverage](https://coveralls.io/repos/github/rasmushenningsson/KernelDensitySJ.jl/badge.svg?branch=master)](https://coveralls.io/github/rasmushenningsson/KernelDensitySJ.jl?branch=master)

Bandwidth selection using Sheather and Jones algorithm also known as width.SJ or bw.SJ. (The specific version of the algorithm is sometimes referred to as "solve-the-equation".)
Only Gaussian kernels are supported.

The aim of this package is to provide implementation of the Sheather and Jones algorithm for bandwidth selection that is:
* Accurate
* Reasonably fast


## Installation
Start the Julia REPL, press ']' to enter package mode and type
```
add https://github.com/rasmushenningsson/KernelDensitySJ.jl.git
```

## Usage Example
```julia
using KernelDensitySJ
X = randn(1000);
bwsj(X)
```

## Technical Details
The bandwidth is computed by solving equation (12) in [Sheather and Jones, 1991] using root-finding methods.
To solve the equation both accurately and fast, this packages computes and successively refines lower and upper bounds of the target function for a given bandwidth until the bounds no longer includes 0, which is enough to bracket the root.

## References
* Sheather, S. J., & Jones, M. C. (1991). A reliable data‐based bandwidth selection method for kernel density estimation. Journal of the Royal Statistical Society: Series B (Methodological), 53(3), 683-690.
