module Polymer

using ArgCheck
using REPL: symbol_latex
using LaTeXStrings

include("utils.jl")
export unicodesymbol2string

include("parameters.jl")
export AbstractParameter, PolymerParameter
export χParam, NParam, χNParam, fParam, ϕParam, RgParam, CParam, bParam, αParam, τParam
export description, as_variable_name, as_ascii_label, as_plot_label

include("types.jl")
export PolymerSystemType, NeatPolymer, PolymerBlend, PolymerSolution
export PolymerType, Homopolymer, Copolymer, BlockCopolymer, RandomCopolymer
export PolymerArchitecture, LinearArchitecture, BranchedArchitecture, StarArchitecture, CombArchitecture, RingArchitecture
export SpaceDimension, D1, D2, D3
export ConfinementType, BulkConfinement, BoxConfinement, SlabConfinement, DiskConfinement, SphereConfinement, CylinderConfinement
export ChargedType, Neutral, SmearedCharge, DiscreteCharge
export BlockEnd, FreeEnd, BranchPoint, PolymerBlock
export AbstractSpecie, KuhnSegment, SmallMolecule
export AbstractComponent, PolymerComponent, SmallMoleculeComponent, ParticleComponent, GiantMoleculeComponent, PolymerSystem
export islinearchain, isconfined, ischarged, multicomponent, ncomponents, species, nspecies, systemtype

include("systems.jl")
export homopolymer_chain, diblock_chain

end # module
