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
export PolymerArchitecture, LinearArchitecture, BranchedArchitecture, StarArchitecture, CombArchitecture, RingArchitecture
export SpaceDimension, D1, D2, D3
export ConfinementType, BulkConfinement, BoxConfinement, SlabConfinement, DiskConfinement, SphereConfinement, CylinderConfinement
export ChargedType, Neutral, SmearedCharge, DiscreteCharge
export BlockEnd, FreeEnd, BranchPoint, PolymerBlock
export AbstractSpecie, KuhnSegment
export AbstractMolecule, SmallMolecule, AbstractPolymer, BlockCopolymer, RandomCopolymer, AlternatingCopolymer, Particle, GiantMolecule
export  AbstractComponent, Component, AbstractSystem, PolymerSystem
export islinearchain, isconfined, ischarged, isfreeblockend, multicomponent, ncomponents, specie, species, nspecies, systemtype

include("systems.jl")
export homopolymer_chain, diblock_chain, solvent, AB_system, AB_A_system, AB_S_system

# include("graph.jl")

end # module
