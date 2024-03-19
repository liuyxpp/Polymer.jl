abstract type AbstractχNMatrix{T} <: AbstractMatrix{T} end

struct χNMatrix{T} <: AbstractχNMatrix{T}
    map::Dict{Set{Symbol}, T}
    mat::Matrix{T}
    imat::Matrix{T}  # the inverse of mat
end

function _species(map::Dict{Set{Symbol}, T}) where T<:Real
    sps = Symbol[]
    for key in keys(map)
        append!(sps, key)
    end
    return unique(sps) |> sort
end

function _check_consistency(map::Dict{Set{Symbol}, T}) where T<:Real
    # Now sps is a list of sorted and unique symbols.
    # which make sure following loops sp1 and sp2 are distinct.
    sps = _species(map)
    for i in 1:length(sps)
        sp1 = sps[i]
        for j in (i+1):length(sps)
            sp2 = sps[j]
            !haskey(map, Set((sp1,sp2))) && return false
        end
    end
    return true
end

function _χNmap_to_matrix(map::Dict{Set{Symbol}, T}) where T<:Real
    _check_consistency(map) || error("Full map of interactions should be specified!")
    sps = _species(map)
    nc = length(sps)
    mat = zeros(T, nc, nc)
    for i in 1:nc
        sp1 = sps[i]
        for j in (i+1):nc
            sp2 = sps[j]
            χN = map[Set((sp1,sp2))]
            mat[i, j] = χN
            mat[j, i] = χN
        end
    end
    r = rank(mat)
    (r == size(mat)[1]) || error("χN map is singular! Possible indistinguishable speices are presented in the map. Consider to combine them.")
    return mat
end

function χNMatrix(map::Dict{Set{Symbol}, T}) where T<:Real
    mat = _χNmap_to_matrix(map)
    imat = inv(mat)
    e1 = first(imat)
    mapnew = Dict(keys(map) .=> [promote(v, e1)[1] for v in values(map)])
    return χNMatrix(mapnew, promote(mat, imat)...)
end

"""
    χNMatrix(m)

Convert a dictionary of interaction pairs to an interaction matrix. The input argument `m` can be a dictionary of `Dict{Tuple{Symobl}, <:Real}`, `Dict{Tuple{String}, <:Real}`, `Dict{<:AbstractArray{Symbol}, <:Real}`, `Dict{<:AbstractArray{String}, <:Real}`.

Example:

```julia
map = Dict((:A, :B)=>10.0, (:A, :C)=>20.0, (:B, :C)=>30.0)
map = Dict(("A", "B")=>10.0, ("A", "C")=>20.0, ("B", "C")=>30.0)
map = Dict([:A, :B]=>10, [:A, :C]=>10.0, [:B, :C]=>20.0)
map = Dict(["A", "B"]=>10.0, ["A", "C"]=>20.0, ["B", "C"]=>30.0)
```
"""
function χNMatrix(m)
    mapnew = Dict(Set.([Symbol.(k) for k in keys(m)]) .=> promote(values(m)...))
    return χNMatrix(mapnew)
end

species(m::χNMatrix) = _species(m.map)

Base.size(A::χNMatrix) = size(A.mat)
Base.getindex(A::χNMatrix, i::Int) = getindex(A.mat, i)
Base.getindex(A::χNMatrix, I::Vararg{Int, N}) where N = getindex(A.mat, I...)

Base.getindex(::χNMatrix{T}, ::Symbol) where T = zero(T)

function Base.getindex(A::χNMatrix{T}, sp1::Symbol, sp2::Symbol) where T
    (sp1 == sp2) && return zero(T)
    return A.map[Set((sp1, sp2))]
end

function _setindex!(A::χNMatrix, v, sp1::Symbol, sp2::Symbol)
    A.map[Set((sp1, sp2))] = v
    # call _χNmap_to_matrix to validate the matrix
    mat = _χNmap_to_matrix(A.map)
    # If no error, set the value
    # update the inverse matrix
    A.mat .= mat
    A.imat .= inv(mat)
end

function _setindex!(A::χNMatrix, v, i, j)
    # Do not set value for j = k
    (i == j) && return nothing
    sps = _species(A.map)
    _setindex!(A, v, sps[i], sps[j])
end

function Base.setindex!(A::χNMatrix, v, i::Int)
    j, k = Tuple(CartesianIndices(A)[i])
    _setindex!(A, v, j, k)
    return nothing
end

function Base.setindex!(A::χNMatrix, v, i::Int, j::Int)
    _setindex!(A, v, i, j)
    return nothing
end

Base.setindex!(::χNMatrix, v, ::Symbol) = nothing

function Base.setindex!(A::χNMatrix, v, sp1::Symbol, sp2::Symbol)
    _setindex!(A, v, sp1, sp2)
end

Base.eltype(::χNMatrix{T}) where T = T

Base.showarg(io::IO, A::χNMatrix, toplevel) = print(io, typeof(A), " of species: ", _species(A.map))

function Base.show(io::IO, m::χNMatrix)
    println(io, "Flory-Huggins interaction parameters betwen species:")
    for (k, v) in m.map
        sp1, sp2 = k
        println(io, "    ($sp1, $sp2) => $v")
    end
end

function Base.show(io::IO, ::MIME"text/plain", m::χNMatrix)
    println(io, "Flory-Huggins interaction parameters between species:")
    for (k, v) in m.map
        sp1, sp2 = k
        println(io, "    ($sp1, $sp2) => $v")
    end
    println(io, "as a matrix:")
    show(io, "text/plain", m.mat)
end