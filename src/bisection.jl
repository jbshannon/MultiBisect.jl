struct BisectionGrid{T <: Number, N}
    signs::Array{Tuple{Bool, Bool}, N}
    evaluated::Array{Bool, N}
    toevaluate::Array{Bool, N}
    domain::NTuple{N, AbstractRange{T}}
    evaluations::Vector{T}
    evalinds::Vector{CartesianIndex{N}}
end

function BisectionGrid(grid)
    sz = size(CartesianIndices(eachindex.(grid)))
    signs = fill((true, true), sz)  # All start as unevaluated
    toevaluate = ones(Bool, sz)
    evaluated = zeros(Bool, sz)
    domain = convert.(promote_type(typeof.(grid)...), grid)
    N = length(domain)
    T = eltype(first(domain))
    evaluations = T[]
    evalinds = CartesianIndex{N}[]
    return BisectionGrid{T, N}(signs, evaluated, toevaluate, domain, evaluations, evalinds)
end

## Helper functions for Tuple{Bool, Bool} sign system

# Functions on Tuple{Bool, Bool} directly
ispositive(sign::Tuple{Bool, Bool}) = sign == (true, false)
isnegative(sign::Tuple{Bool, Bool}) = sign == (false, true)
iszerosign(sign::Tuple{Bool, Bool}) = sign == (false, false)  # Renamed to avoid conflict with Base.iszero
isevaluated(sign::Tuple{Bool, Bool}) = sign != (true, true)

# Edge detection including zero-adjacent edges
function differentsigns(signs::Array{Tuple{Bool, Bool}, N}, I1, I2) where N
    s1, s2 = signs[I1], signs[I2]
    # Different if at least one boolean differs AND both are evaluated
    return (s1 != s2) && isevaluated(s1) && isevaluated(s2)
end

function countedges(BG::BisectionGrid)
    S = BG.signs
    CI = CartesianIndices(S)
    s = 0
    for I in CI, d in forwardinds(length(I))
        inside = min(I+d, last(CI)) != I
        inside && differentsigns(S, I, I+d) && (s += 1)
    end
    return s
end

function show(io::IO, G::BisectionGrid{T, N}) where {T, N}
    println(io, typeof(G))
    println(io, "       Domain: $(G.domain)")
    println(io, "  Grid points: $(length(G.signs))")
    println(io, "  Evaluations: $(sum(G.evaluated))")
    print(  io, "        Edges: $(countedges(G))")
end

## Property retrieval API

evaluated(BG::BisectionGrid) = BG.evaluated
evaluations(BG::BisectionGrid) = sum(evaluated(BG))
domain(BG::BisectionGrid) = BG.domain
efficiency(BG::BisectionGrid) = 1 - evaluations(BG)/length(evaluated(BG))


"""
    splitsign(BG::BisectionGrid)

Split the evaluated points of the bisection grid into a vector of positive points and a vector of negative points.

# Examples
```jldoctest
julia> BG = bisect(x -> 1 - sum(abs2, x), (0.0:1.0, 0.0:1.0); iterations=2)
BisectionGrid{Float64, 2}
       Domain: (0.0:0.5:1.0, 0.0:0.5:1.0)
  Grid points: 9
  Evaluations: 9
        Edges: 4

julia> posx, negx = splitsign(BG);


julia> posx
4-element Vector{Tuple{Float64, Float64}}:
 (0.0, 0.0)
 (0.5, 0.0)
 (0.0, 0.5)
 (0.5, 0.5)

julia> negx
5-element Vector{Tuple{Float64, Float64}}:
 (1.0, 0.0)
 (1.0, 0.5)
 (0.0, 1.0)
 (0.5, 1.0)
 (1.0, 1.0)
```
"""
function splitsign(BG::BisectionGrid)
    E = evaluated(BG)
    S = BG.signs
    CI = Iterators.filter(I -> E[I], CartesianIndices(E)) # evaluated indices
    posinds = eltype(CI)[]
    neginds = eltype(CI)[]
    foreach(I -> begin
        sign = S[I]
        if ispositive(sign)
            push!(posinds, I)
        elseif isnegative(sign)
            push!(neginds, I)
        end
        # Zeros are excluded from both arrays
    end, CI)
    _getindex = Base.Fix1(domainindex, domain(BG))
    return (_getindex.(posinds), _getindex.(neginds))
end

## Expand

function filldouble!(B, A)
    C = CartesianIndices(B)
    CF = first(C)
    CL = last(C)
    C2 = 2*oneunit(CF)
    for (AC, BC) in enumerate(CF:C2:CL)
        B[BC] = A[AC]
    end
    return B
end

# find the size of the doubled array
doublesize(A) = ntuple(d -> 1+2*(size(A, d)-1), ndims(A))

"""
    expandzeros(A)

Expand the array A, inserting a zero between each element.

# Examples

```jldoctest
julia> A = reshape(1:4, 2, 2)
2×2 reshape(::UnitRange{Int64}, 2, 2) with eltype Int64:
 1  3
 2  4

julia> expandzeros(A)
3×3 Matrix{Int64}:
 1  0  3
 0  0  0
 2  0  4

julia> B = reshape(1:6, 3, 2)
3×2 reshape(::UnitRange{Int64}, 3, 2) with eltype Int64:
 1  4
 2  5
 3  6

julia> expandzeros(B)
5×3 Matrix{Int64}:
 1  0  4
 0  0  0
 2  0  5
 0  0  0
 3  0  6
```
"""
function expandzeros(A)
    if eltype(A) == Tuple{Bool, Bool}
        B = fill((false, false), doublesize(A))  # Fill with zero values for expansion
    elseif eltype(A) == Bool
        B = zeros(eltype(A), doublesize(A))  # Original behavior for Bool arrays
    else
        B = zeros(eltype(A), doublesize(A))  # Fallback for other types
    end
    filldouble!(B, A)
end

function expandones(A)
    B = fill((true, true), doublesize(A))  # Fill with unevaluated values
    filldouble!(B, A)
end

function expanddomain(X)
    ntuple(length(X)) do i
        r = X[i]
        return range(first(r), last(r); length=1+2*(length(r)-1))
    end
end

function expandinds!(evalinds, domain)
    L = length.(domain)
    N = length(domain)
    newI = ntuple(d -> 1:2:(2L[d]-1), N)
    for i in eachindex(evalinds)
        I = Tuple(evalinds[i])
        evalinds[i] = CartesianIndex(ntuple(j -> newI[j][I[j]], N))
    end
    return evalinds
end

function expand(G::BisectionGrid)
    return BisectionGrid(
        expandzeros(G.signs),
        expandzeros(G.evaluated),  # This works for Bool arrays
        expandzeros(G.toevaluate),  # This works for Bool arrays
        expanddomain(G.domain),
        G.evaluations,
        expandinds!(G.evalinds, G.domain),
    )
end

## Hypercubes

# iterate over all elements of the hypercube
hypercube(H; dist=2) = range(H, H + dist*oneunit(H))

"""
    corners(CI::CartesianIndices)

Construct an iterator over the corners of the hypercube defined by `CI`.

# Examples
```jldoctest
julia> CI = CartesianIndices((1:3, 1:4))
CartesianIndices((1:3, 1:4))

julia> collect(CI)
3×4 Matrix{CartesianIndex{2}}:
 CartesianIndex(1, 1)  CartesianIndex(1, 2)  …  CartesianIndex(1, 4)
 CartesianIndex(2, 1)  CartesianIndex(2, 2)     CartesianIndex(2, 4)
 CartesianIndex(3, 1)  CartesianIndex(3, 2)     CartesianIndex(3, 4)

julia> corners(CI)
CartesianIndices((1:2:3, 1:3:4))

julia> collect(corners(CI))
2×2 Matrix{CartesianIndex{2}}:
 CartesianIndex(1, 1)  CartesianIndex(1, 4)
 CartesianIndex(3, 1)  CartesianIndex(3, 4)
```
"""
function corners(CI::CartesianIndices)
    F = first(CI)
    L = last(CI)
    return F:L-F:L
end

# find starting indices of all hypercubes in the expanded array
"""
    findhypercubes(A; dist=2)

Find the starting vertex of all hypercubes in the array `A`.

# Examples
```jldoctest
julia> A = ones(Bool, 5, 5)
5×5 Matrix{Bool}:
 1  1  1  1  1
 1  1  1  1  1
 1  1  1  1  1
 1  1  1  1  1
 1  1  1  1  1

julia> cubes = findhypercubes(A)
CartesianIndices((1:2:3, 1:2:3))

julia> collect(cubes)
2×2 Matrix{CartesianIndex{2}}:
 CartesianIndex(1, 1)  CartesianIndex(1, 3)
 CartesianIndex(3, 1)  CartesianIndex(3, 3)

julia> cubes3 = findhypercubes(ones(Bool, 5, 5, 5))
CartesianIndices((1:2:3, 1:2:3, 1:2:3))

julia> collect(cubes3)
2×2×2 Array{CartesianIndex{3}, 3}:
[:, :, 1] =
 CartesianIndex(1, 1, 1)  CartesianIndex(1, 3, 1)
 CartesianIndex(3, 1, 1)  CartesianIndex(3, 3, 1)

[:, :, 2] =
 CartesianIndex(1, 1, 3)  CartesianIndex(1, 3, 3)
 CartesianIndex(3, 1, 3)  CartesianIndex(3, 3, 3)
```
"""
function findhypercubes(A; dist=2)
    CI = CartesianIndices(A)
    F = first(CI)
    L = last(CI)
    Δ = dist*oneunit(F)
    return F:Δ:L-Δ
end

equalcorners(A, I; dist=2) = all(i -> isevaluated(A[i]) && A[i] == A[I], corners(hypercube(I; dist)))

isadjacent(I, J; dist=2) = abs(sum(Tuple(I-J))) == dist

# if corner values are all the same, fill the hypercube with those values
function fillhypercube!(A, I, val; dist=2)
    H = hypercube(I; dist)
    for i in Iterators.filter(!in(corners(H)), H)
        A[i] = val
    end
    return A
end

## Iterations

function fillnonbracketing!(signs, toevaluate; monotonic=false)
    cubes = findhypercubes(signs)
    for H in cubes
        if equalcorners(signs, H)
            fillhypercube!(signs, H, signs[H]) # fill in sign where it is known
        else
            fillhypercube!(toevaluate, H, true) # mark points to be evaluated
        end
    end
    if monotonic # handle points on edges that have the same sign
        for H in cubes, V in corners(hypercube(H))
            H_sign = signs[H]
            V_sign = signs[V]
            if isadjacent(V, H) && isevaluated(H_sign) && isevaluated(V_sign) && H_sign == V_sign
                midpoint = (H:V)[2] # index between two vertices
                signs[midpoint] = H_sign # sign is known (same as both edges)
                toevaluate[midpoint] = false # do not evaluate
            end
        end
    end
    return signs, toevaluate
end

function fillnonbracketing!(BG::BisectionGrid; kwargs...)
    fillnonbracketing!(BG.signs, BG.toevaluate; kwargs...)
    return BG
end

function _evaluate!(f, signs, toevaluate, evaluated, domain, evaluations, evalinds, I)
    val = f(ntuple(j -> domain[j][Tuple(I)[j]], ndims(signs)))
    if val > 0
        signs[I] = (true, false)
    elseif val < 0
        signs[I] = (false, true)
    else  # val == 0
        signs[I] = (false, false)
    end
    evaluated[I] = true
    toevaluate[I] = false
    push!(evaluations, val)
    push!(evalinds, I)
    return val
end

function evaluate!(f, signs, toevaluate, evaluated, domain, evaluations, evalinds; threaded=false)
    if threaded
        @threads for I in filter(I -> toevaluate[I], CartesianIndices(signs))
            _evaluate!(f, signs, toevaluate, evaluated, domain, evaluations, evalinds, I)
        end
    else
        for I in Iterators.filter(I -> toevaluate[I], CartesianIndices(signs))
            _evaluate!(f, signs, toevaluate, evaluated, domain, evaluations, evalinds, I)
        end
    end
    return signs, evaluated
end

function evaluate!(f, BG::BisectionGrid; kwargs...)
    evaluate!(f, BG.signs, BG.toevaluate, BG.evaluated, BG.domain, BG.evaluations, BG.evalinds; kwargs...)
end

## Revert
# Only useful if you want to try a more expensive interpolation method

# function _revert(A::AbstractArray; steps=1)
#     CI = CartesianIndices(A)
#     IF = first(CI)
#     IL = last(CI)
#     Δ = 2^steps * one(IF)
#     return A[IF:Δ:IL]
# end

# function _revert(X::NTuple{N, R}; steps=1) where {N, R <: AbstractRange}
#     ntuple(length(X)) do i
#         r = X[i]
#         newlength = 1 + (length(r) - 1) ÷ (2^steps)
#         return range(first(r), last(r); length=newlength)
#     end
# end


# """
#     revert(BG::BisectionGrid; steps=1)

# Revert the grid `BG` by a given number of `steps` (e.g. go from the fourth iteration back to the third iteration).

# # Examples

# ```jldoctest
# julia> f(x) = 1 - sum(abs2, x);

# julia> grid = (0.0:1.0, 0.0:1.0);

# julia> BG3 = bisect(f, grid; iterations=3)
# BisectionGrid{Float64, 2}
#        Domain: (0.0:0.25:1.0, 0.0:0.25:1.0)
#   Grid points: 25
#   Evaluations: 22
#         Edges: 8


# julia> BG4 = bisect(f, grid; iterations=4)
# BisectionGrid{Float64, 2}
#        Domain: (0.0:0.125:1.0, 0.0:0.125:1.0)
#   Grid points: 81
#   Evaluations: 51
#         Edges: 16


# julia> revert(BG4)
# BisectionGrid{Float64, 2}
#        Domain: (0.0:0.25:1.0, 0.0:0.25:1.0)
#   Grid points: 25
#   Evaluations: 22
#         Edges: 8
# ```
# """
# function revert(BG::BisectionGrid; steps=1)
#     BisectionGrid(
#         _revert(BG.signs; steps),
#         _revert(BG.evaluated; steps),
#         _revert(BG.toevaluate; steps),
#         _revert(BG.domain; steps),
#     )
# end

## Bisection function

"""
    bisect(f, grid; threaded=false, monotonic=false, iterations=5, verbose=false)

Perform the multidimensional bisection algorithm of function `f` on an initial `grid` of
evaluation points.

# Arguments

## Positional Arguments

- `f`: function to be evaluated
- `grid`: a tuple of ranges describing the initial evaluation grid

## Keyword Arguments

- `threaded`: set to `true` for multithreaded evaluation
- `monotonic`: set to true if `f` is known to be monotonic in each dimension
- `iterations`: number of times to iterate the algorithm
- `verbose`: set to `true` to print the number of new function evaluations at each iteration

# Examples

```jldoctest
julia> grid = (0.0:1.0, 0.0:1.0);

julia> bisect(x -> 1 - sum(abs2, x), grid)
BisectionGrid{Float64, 2}
       Domain: (0.0:0.0625:1.0, 0.0:0.0625:1.0)
  Grid points: 289
  Evaluations: 112
        Edges: 32


julia> bisect(x -> x[2] - x[1], grid)
BisectionGrid{Float64, 2}
       Domain: (0.0:0.0625:1.0, 0.0:0.0625:1.0)
  Grid points: 289
  Evaluations: 112
        Edges: 32


julia> bisect(x -> x[2] - x[1], grid; monotonic=true)
BisectionGrid{Float64, 2}
       Domain: (0.0:0.0625:1.0, 0.0:0.0625:1.0)
  Grid points: 289
  Evaluations: 67
        Edges: 32
```

# Extended help

## Algorithm

The algorithm begins with an initial evaluation grid of ``N``-dimensional hypercubes (lines, squares, cubes, etc.) and proceeds by splitting it in "half" at each stage, but only in regions where the function is known to change sign. After the function is evaluated, each hypercube is split in half (a line into two lines, a square into four squares, a cube into eight cubes, etc.) and the function sign at each vertex is checked. If not all vertices share the same sign, the function must change sign somewhere within the hypercube, so the vertices of the subcubes (e.g. the midpoint of a line) will be marked for evaluation in the next stage; if the vertices all share the same sign, the function is presumed to have the same sign everywhere throughout the hypercube, so the sign is inferred at all interior points and the function is not evaluated in the hypercube again.

By passing the keyword `monotonic=true`, the algorithm further assumes that the function cannot change sign in any one direction between two points with the same sign. This means that the sign of the function should be checked at the vertices of each edge of the hypercube; if the signs are the same, the function is assumed to have the same sign at the midpoint along the edge and will not be evaluated.

## Multithreading

Set the keyword argument `threaded=true` to make evaluation of the function at each iteration multithreaded. For functions that are very cheap to evaluate, the overhead from setting up threads will actually make each iteration longer. However, this effect disappears quickly.

### Example
```julia-repl
julia> ffast(x) = 1 - sum(abs2, x);

julia> grid = (0.0:1.0, 0.0:1.0);

julia> @btime bisect(ffast, \$grid; threaded=false);
  49.720 μs (40 allocations: 3.59 KiB)

julia> @btime bisect(ffast, \$grid; threaded=true);
  146.173 μs (163 allocations: 18.23 KiB)

julia> fslow(x) = (sleep(0.1); ffast(x))
fslow (generic function with 1 method)

julia> @btime bisect(fslow, \$grid; threaded=false);
  11.397 s (605 allocations: 20.83 KiB)

julia> @btime bisect(fslow, \$grid; threaded=true);
  3.170 s (736 allocations: 35.45 KiB)
```

## Monotonicity

If the function `f` is known to be monotonic, then the number of evaluations can be reduced by inferring the behavior of `f` along edges where both vertices have the same sign.

### Example
```jldoctest
julia> bisect(x -> 1 - sum(abs2, x), (0.0:1.0, 0.0:1.0); monotonic=false)
BisectionGrid{Float64, 2}
       Domain: (0.0:0.0625:1.0, 0.0:0.0625:1.0)
  Grid points: 289
  Evaluations: 112
        Edges: 32


julia> bisect(x -> 1 - sum(abs2, x), (0.0:1.0, 0.0:1.0); monotonic=true)
BisectionGrid{Float64, 2}
       Domain: (0.0:0.0625:1.0, 0.0:0.0625:1.0)
  Grid points: 289
  Evaluations: 80
        Edges: 32
```

## Iterations

For finer grids, repeat the algorithm for more iterations.

### Example
```jldoctest
julia> bisect(x -> 1 - sum(abs2, x), (0.0:1.0, 0.0:1.0))
BisectionGrid{Float64, 2}
       Domain: (0.0:0.0625:1.0, 0.0:0.0625:1.0)
  Grid points: 289
  Evaluations: 112
        Edges: 32


julia> bisect(x -> 1 - sum(abs2, x), (0.0:1.0, 0.0:1.0); iterations=3)
BisectionGrid{Float64, 2}
       Domain: (0.0:0.25:1.0, 0.0:0.25:1.0)
  Grid points: 25
  Evaluations: 22
        Edges: 8


julia> bisect(x -> 1 - sum(abs2, x), (0.0:1.0, 0.0:1.0); iterations=7)
BisectionGrid{Float64, 2}
       Domain: (0.0:0.015625:1.0, 0.0:0.015625:1.0)
  Grid points: 4225
  Evaluations: 490
        Edges: 128
```

## Logging

For longer-running bisections, it may be useful to report the number of new evaluations as a kind of thread-safe progress meter.

### Example
```jldoctest
julia> bisect(x -> 1 - sum(abs2, x), (0.0:1.0, 0.0:1.0); verbose=true)
[ Info: Iteration 1: 4 gridpoints (4 evaluations)
[ Info: Iteration 2: 9 gridpoints (5 evaluations)
[ Info: Iteration 3: 25 gridpoints (13 evaluations)
[ Info: Iteration 4: 81 gridpoints (29 evaluations)
[ Info: Iteration 5: 289 gridpoints (61 evaluations)
BisectionGrid{Float64, 2}
       Domain: (0.0:0.0625:1.0, 0.0:0.0625:1.0)
  Grid points: 289
  Evaluations: 112
        Edges: 32


julia> bisect(x -> 1 - sum(abs2, x), (0.0:1.0, 0.0:1.0); verbose=true, iterations=7)
[ Info: Iteration 1: 4 gridpoints (4 evaluations)
[ Info: Iteration 2: 9 gridpoints (5 evaluations)
[ Info: Iteration 3: 25 gridpoints (13 evaluations)
[ Info: Iteration 4: 81 gridpoints (29 evaluations)
[ Info: Iteration 5: 289 gridpoints (61 evaluations)
[ Info: Iteration 6: 1089 gridpoints (125 evaluations)
[ Info: Iteration 7: 4225 gridpoints (253 evaluations)
BisectionGrid{Float64, 2}
       Domain: (0.0:0.015625:1.0, 0.0:0.015625:1.0)
  Grid points: 4225
  Evaluations: 490
        Edges: 128
```
"""
function bisect(f, grid; threaded=false, monotonic=false, iterations=5, verbose=false)
    BG = BisectionGrid(grid)
    for i in 1:iterations
        if i > 1
            BG = expand(BG)
            fillnonbracketing!(BG; monotonic)
        end
        verbose && @info("Iteration $i: $(length(BG.signs)) gridpoints ($(sum(BG.toevaluate)) evaluations)")
        evaluate!(f, BG; threaded)
    end
    return BG
end
