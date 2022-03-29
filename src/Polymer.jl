module Polymer

using LinearAlgebra
using ArgCheck
using REPL: symbol_latex
using LaTeXStrings
using YAML
using Setfield

include("utils.jl")
export
    unicodesymbol2string,
    reverse_dict

include("parameters.jl")
export
    AbstractParameter,
    PolymerParameter
export
    χType,
    NType,
    χNType,
    fType,
    ϕType,
    RgType,
    CType,
    bType,
    αType,
    τType,
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
    value_type,
    variable_symbol,
    as_variable_name,
    as_ascii_label,
    as_plot_label
export
    AbstractControlParameter,
    ϕControlParameter,
    αControlParameter,
    fControlParameter,
    χNControlParameter,
    bControlParameter

include("chiN.jl")
export
    AbstractχNMatrix,
    χNMatrix
    # Not exported but useful:
    # χN, χNmap

include("types.jl")
export
    PolymerSystemType,
    NeatPolymer,
    PolymerBlend,
    PolymerSolution
export
    ComponentNumberType,
    MonoSystem,
    BinarySystem,
    TernarySystem,
    MultiComponentSystem
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

include("properties.jl")
export
    isconfined,
    ischarged,
    isfreeblockend,
    multicomponent,
    ncomponents,
    specie, species, nspecies,
    nblocks,
    systemtype,
    component_number_type,
    component_label, component_labels,
    component_id, component_ids,
    molecule_id, molecule_ids,
    molecule_label, molecule_labels,
    block_id, block_ids,
    block_label, block_labels,
    block_length, block_lengths,
    block_specie, block_species,
    block_b, block_bs
    # not exported but useful functions:
    # ϕ, ϕs, α, αs, molecules, molecule,
    # blocks, block
    # getparam

include("systems.jl")
export
    homopolymer_chain,
    diblock_chain,
    solvent,
    AB_system,
    AB_A_system,
    AB_S_system,
    A_B_S_system,
    A_B_S1_S2_system

include("make.jl")
export
    ObjType
export
    load_config,
    make

include("update.jl")
# update! is not exported because of possible name confilication
# setparam! is an alias of update!

end # module
