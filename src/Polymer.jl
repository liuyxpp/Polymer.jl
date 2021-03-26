module Polymer

using LinearAlgebra
using ArgCheck
using REPL: symbol_latex
using LaTeXStrings
using YAML

include("utils.jl")
export
    unicodesymbol2string,
    reverse_dict

include("parameters.jl")
export
    AbstractParameter,
    PolymerParameter
export
    χParam,
    NParam,
    χNParam,
    fParam,
    ϕParam,
    RgParam,
    CParam,
    bParam,
    αParam,
    τParam,
    DEFAULT_PARAMETERS
export
    description,
    as_variable_name,
    as_ascii_label,
    as_plot_label

include("chiN.jl")
export
    AbstractχNMatrix,
    χNMatrix

include("types.jl")
export
    PolymerSystemType,
    NeatPolymer, PolymerBlend,
    PolymerSolution
export
    SpaceDimension,
    D1, D2, D3
export
    ConfinementType,
    BulkConfinement,
    BoxConfinement,
    SlabConfinement,
    DiskConfinement,
    SphereConfinement,
    CylinderConfinement
export
    ChargedType,
    Neutral,
    SmearedCharge,
    DiscreteCharge
export
    BlockEnd,
    FreeEnd,
    BranchPoint,
    PolymerBlock
export
    AbstractSpecie,
    KuhnSegment
export
    AbstractMolecule,
    SmallMolecule,
    AbstractPolymer,
    BlockCopolymer,
    RandomCopolymer,
    AlternatingCopolymer,
    Particle,
    GiantMolecule
export
    AbstractComponent,
    Component,
    AbstractSystem,
    PolymerSystem
export
    isconfined,
    ischarged,
    isfreeblockend,
    multicomponent,
    ncomponents,
    specie,
    species,
    nspecies,
    nblocks,
    systemtype

include("systems.jl")
export
    homopolymer_chain,
    diblock_chain,
    solvent,
    AB_system,
    AB_A_system,
    AB_S_system

include("make.jl")
export
    ObjType
export
    load_config,
    make

end # module
