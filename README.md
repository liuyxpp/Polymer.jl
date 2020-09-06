# Polymer.jl [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://liuyxpp.github.io/Polymer.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://liuyxpp.github.io/Polymer.jl/dev) [![Build Status](https://github.com/liuyxpp/Polymer.jl/workflows/CI/badge.svg)](https://github.com/liuyxpp/Polymer.jl/actions) [![Build Status](https://travis-ci.com/liuyxpp/Polymer.jl.svg?branch=master)](https://travis-ci.com/liuyxpp/Polymer.jl) [![Build Status](https://ci.appveyor.com/api/projects/status/github/liuyxpp/Polymer.jl?svg=true)](https://ci.appveyor.com/project/liuyxpp/Polymer-jl) [![Coverage](https://codecov.io/gh/liuyxpp/Polymer.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/liuyxpp/Polymer.jl)

**Polymer.jl** defines a common interface to describe a polymer system.

*Warning: Be aware that this package is currently under active development. The interface is highly unstable and subjects to change frequently.*

## Usage

```julia
julia> using Polymer

# Create A and B monomers.
julia> sA = KuhnSegment(:A)
julia> sB = KuhnSegment(:B)

# Create free end and branch point (joint)
julia> eA = FreeEnd(:A1)
julia> eAB = BranchPoint(:AB)
julia> eB = FreeEnd(:B1)

# Create A and B blocks
julia> A = PolymerBlock(:A, sA, 0.5, eA, eAB)
julia> B = PolymerBlock(:B, sB, 0.5, eB, eAB)

# Create a AB diblock copolymer chain
julia> chainAB = BlockCopolymer(:AB, [A,B])

# Create a homopolymer chain
julia> hA = PolymerBlock(:hA, sA, 1.0, FreeEnd(), FreeEnd())
julia> chainA = BlockCopolymer(:hA, [hA])

# Create components
julia> polymerAB = Component(chainAB; ϕ=0.5)
julia> polymerA = Component(chainA; ϕ=0.5)

# Create AB/A polymer blend system.
julia> AB_A = PolymerSystem([polymerAB, polymerA]; χN_map=Dict(Set([:A, :B])=>20.0))
```

Convenient functions are also provided to create common polymer chains and systems. For example, above AB chain, A chain, AB/A polymer blend system can be simply created by a single line of code.

```julia
julia> diblock_chain() # AB chain
julia> homopolymer_chain() # A chain
julia> AB_A_system() # AB/A polymer blend
```

At present, there is no documentation. Consult the testing codes reside in `test` folder to learn more.

## Contribute

* Star the package on [github.com](https://github.com/liuyxpp/Polymer.jl).
* File an issue or make a pull request on [github.com](https://github.com/liuyxpp/Polymer.jl).
* Contact the author via email <lyx@fudan.edu.cn>.

## Links

* [Source code](https://github.com/liuyxpp/Polymer.jl)