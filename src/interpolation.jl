# get CartesianIndex steps forward in each dimension
"""
    forwardinds(d)

Get CartesianIndex steps forward in each of `d` dimensions.

# Examples

```jldoctest
julia> forwardinds(1)
(CartesianIndex(1,),)

julia> forwardinds(2)
(CartesianIndex(1, 0), CartesianIndex(0, 1))

julia> forwardinds(3)
(CartesianIndex(1, 0, 0), CartesianIndex(0, 1, 0), CartesianIndex(0, 0, 1))
```
"""
forwardinds(d) = ntuple(j -> CartesianIndex(ntuple(==(j), d)), d)


"""
    domainindex(X, I::CartesianIndex)

Like `getindex`, except `X` is a tuple of ranges instead of an array.

# Examples
```jldoctest
julia> X = (1.0:5.0, 0.5:2.5);

julia> domainindex(X, CartesianIndex(1, 1))
(1.0, 0.5)

julia> domainindex(X, CartesianIndex(2, 1))
(2.0, 0.5)

julia> domainindex(X, CartesianIndex(1, 2))
(1.0, 1.5)

julia> domainindex(X, last(CartesianIndices(eachindex.(X))))
(5.0, 2.5)
```
"""
domainindex(X, I::CartesianIndex{N}) where N = ntuple(i -> getindex(X[i], I[i]), N)


"""
    edges(BG::BisectionGrid)

Find all hypercube edges where the sign of the function changes. Returns a `Vector` of 2-tuples of points in the domain.
"""
function edges(BG::BisectionGrid{T, N}) where {T, N}
    S = BG.signs
    _getindex = Base.Fix1(domainindex, domain(BG))
    CI = CartesianIndices(S)
    edgevec = NTuple{2, NTuple{N, T}}[] # vector of edges
    for I in CI, d in forwardinds(N)
        inside = min(I+d, last(CI)) != I # the forward step does not go outside the array
        inside && (S[I] != S[I+d]) && push!(edgevec, _getindex.((I, I+d)))
    end
    return edgevec
end


"""
    edgedim(edge)

Return the dimension of the edge (i.e. the only direction in which the vertices vary).

# Examples
```jldoctest
julia> edgedim(((0.5, 1.0), (1.0, 1.0)))
1

julia> edgedim(((0.5, 1.0), (0.5, 1.5)))
2

julia> edgedim(((1.0, 1.0, 1.0), (1.0, 1.0, 2.0)))
3
```
"""
edgedim(edge) = findfirst(!=(0), edge[2] .- edge[1])


"""
    edgebounds(edge)

Return the bounds of the edge in its varying dimension.

# Examples

```jldoctest
julia> edgebounds(((0.5, 1.0), (1.0, 1.0)))
(0.5, 1.0)

julia> edgebounds(((0.5, 1.0), (0.5, 1.5)))
(1.0, 1.5)

julia> edgebounds(((1.0, 1.0, 1.0), (1.0, 1.0, 2.0)))
(1.0, 2.0)
```
"""
edgebounds(edge) = getindex.(edge, edgedim(edge))


"""
    edgetuple(edge)

Return a function converting a value along the edge dimension to a point.

# Examples

```jldoctest
julia> tup = edgetuple(((0.5, 1.0), (1.0, 1.0)));

julia> tup(0.75)
(0.75, 1.0)

julia> tup(0.9)
(0.9, 1.0)

julia> edgetuple(((0.5, 1.0), (0.5, 1.5)))(1.2)
(0.5, 1.2)
```
"""
edgetuple(edge) = x -> ntuple(j -> j == edgedim(edge) ? x : edge[1][j], length(edge[1]))

"""
    edgeroot(f, edge, M; kwargs...)

Find the root along `edge` with a call to `Roots.find_zero`. Both `M` and `kwargs` are
passed directly to `find_zero`.

# Examples
```jldoctest
julia> f(x) = 1 - sum(abs2, x); # unit circle

julia> BG = bisect(f, (0.0:1.0, 0.0:1.0); iterations=3);

julia> E = edges(BG)
8-element Vector{Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}}:
 ((0.75, 0.0), (1.0, 0.0))
 ((0.75, 0.25), (1.0, 0.25))
 ((0.75, 0.5), (1.0, 0.5))
 ((0.75, 0.5), (0.75, 0.75))
 ((0.0, 0.75), (0.0, 1.0))
 ((0.25, 0.75), (0.25, 1.0))
 ((0.5, 0.75), (0.75, 0.75))
 ((0.5, 0.75), (0.5, 1.0))

julia> tracks = Roots.Tracks();

julia> root = edgeroot(f, E[3], Roots.ITP(); tracks)
(0.8660254037844387, 0.5)

julia> f(root)
0.0

julia> tracks
Results of univariate zero finding:

* Converged to: 0.8660254037844387
* Algorithm: Roots.ITP{Float64, Int64}(0.2, 2, 1)
* iterations: 7
* function evaluations ≈ 9
* stopped as f(x_n) = 0

Trace:
(a₀, b₀) = ( 0.75, 1 )
(a₁, b₁) = ( 0.75, 0.86964285714285705 )
(a₂, b₂) = ( 0.86290337975046705, 0.86964285714285705 )
(a₃, b₃) = ( 0.86290337975046705, 0.8660279692952968 )
(a₄, b₄) = ( 0.86602344653979335, 0.8660279692952968 )
(a₅, b₅) = ( 0.86602344653979335, 0.86602540378563075 )
(a₆, b₆) = ( 0.86602540378367243, 0.86602540378563075 )
(a₇, b₇) = ( 0.86602540378443871, 0.86602540378563075 )
```
"""
function edgeroot(f, edge, M; kwargs...)
    iszero(f(edge[1])) && return edge[1]
    iszero(f(edge[2])) && return edge[2]
    root = find_zero(f ∘ edgetuple(edge), edgebounds(edge), M; kwargs...)
    return edgetuple(edge)(root)
end

function edgeroot(f, M; kwargs...)
    return edge -> edgeroot(f, edge, M; kwargs...)
end


"""
    marchingsquares(edge)

Return the midpoint of the edge, in a manner similar to the [marching squares algorithm](https://en.wikipedia.org/wiki/Marching_squares). Does not require any function evaluations.

# Examples

```jldoctest
julia> f(x) = 1 - sum(abs2, x); # unit circle

julia> BG = bisect(f, (0.0:1.0, 0.0:1.0); iterations=3);


julia> E = edges(BG)
8-element Vector{Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}}:
 ((0.75, 0.0), (1.0, 0.0))
 ((0.75, 0.25), (1.0, 0.25))
 ((0.75, 0.5), (1.0, 0.5))
 ((0.75, 0.5), (0.75, 0.75))
 ((0.0, 0.75), (0.0, 1.0))
 ((0.25, 0.75), (0.25, 1.0))
 ((0.5, 0.75), (0.75, 0.75))
 ((0.5, 0.75), (0.5, 1.0))

julia> marchingsquares.(E)
8-element Vector{Tuple{Float64, Float64}}:
 (0.875, 0.0)
 (0.875, 0.25)
 (0.875, 0.5)
 (0.75, 0.625)
 (0.0, 0.875)
 (0.25, 0.875)
 (0.625, 0.75)
 (0.5, 0.875)
```
"""
function marchingsquares(edge)
    d = edgedim(edge)
    midpoint = (edge[1][d] + edge[2][d])/2
    return edgetuple(edge)(midpoint)
end

"""
    linearroot(f, edge)

Construct a linear interpolation of `f` through `edge` and return its root. Requires two
function evaluations.

# Examples

```jldoctest
julia> f(x) = 1 - sum(abs2, x); # unit circle

julia> BG = bisect(f, (0.0:1.0, 0.0:1.0); iterations=3);

julia> E = edges(BG)
8-element Vector{Tuple{Tuple{Float64, Float64}, Tuple{Float64, Float64}}}:
 ((0.75, 0.0), (1.0, 0.0))
 ((0.75, 0.25), (1.0, 0.25))
 ((0.75, 0.5), (1.0, 0.5))
 ((0.75, 0.5), (0.75, 0.75))
 ((0.0, 0.75), (0.0, 1.0))
 ((0.25, 0.75), (0.25, 1.0))
 ((0.5, 0.75), (0.75, 0.75))
 ((0.5, 0.75), (0.5, 1.0))

julia> linearroot.(f, E)
8-element Vector{Tuple{Float64, Float64}}:
 (1.0, 0.0)
 (0.9642857142857143, 0.25)
 (0.8571428571428571, 0.5)
 (0.75, 0.65)
 (0.0, 1.0)
 (0.25, 0.9642857142857143)
 (0.65, 0.75)
 (0.5, 0.8571428571428571)
```
"""
function linearroot(f, edge)
    x = edgebounds(edge)
    y = f.(edge)
    b = (y[2] - y[1])/(x[2] - x[1])
    a = y[1] - b*x[1]
    return edgetuple(edge)(-a/b)
end

linearroot(f) = edge -> linearroot(f, edge)

"""
    interpolate(rootfinder, BG::BisectionGrid{T, N}; threaded=false) where {T, N}

Compute the roots along all bracketing edges of `BG`.

The function `rootfinder` can be any method that receives an edge and returns a root, such
as [`edgeroot`](@ref), [`linearroot`](@ref), or [`marchingsquares`](@ref).

# Examples

```jldoctest
julia> f(x) = 1 - sum(abs2, x); # unit circle

julia> BG = bisect(f, (0.0:1.0, 0.0:1.0); iterations=3);

julia> interpolate(marchingsquares, BG)
8-element Vector{Tuple{Float64, Float64}}:
 (0.875, 0.0)
 (0.875, 0.25)
 (0.875, 0.5)
 (0.75, 0.625)
 (0.0, 0.875)
 (0.25, 0.875)
 (0.625, 0.75)
 (0.5, 0.875)

julia> interpolate(edge -> linearroot(f, edge), BG)
8-element Vector{Tuple{Float64, Float64}}:
 (1.0, 0.0)
 (0.9642857142857143, 0.25)
 (0.8571428571428571, 0.5)
 (0.75, 0.65)
 (0.0, 1.0)
 (0.25, 0.9642857142857143)
 (0.65, 0.75)
 (0.5, 0.8571428571428571)

julia> interpolate(edge -> edgeroot(f, edge, Roots.ITP()), BG)
8-element Vector{Tuple{Float64, Float64}}:
 (1.0, 0.0)
 (0.9682458365518543, 0.25)
 (0.8660254037844387, 0.5)
 (0.75, 0.6614378277661477)
 (0.0, 1.0)
 (0.25, 0.9682458365518543)
 (0.6614378277661477, 0.75)
 (0.5, 0.8660254037844387)
```
"""
function interpolate(rootfinder, BG::BisectionGrid; threaded=false)
    if threaded
        spawned = map(edge -> @spawn(rootfinder(edge)), edges(BG))
        return fetch.(spawned)
    else
        return rootfinder.(edges(BG))
    end
end

"""
    interpolate(BG::BisectionGrid{T, N}) where {T, N}

Compute the roots along all bracketing edges of `BG` with a linear approximation using precomputed function evaluations. Equivalent to `interpolate(linearroot(f), BG)`, but does not require re-evaluating the function at each edge.

# Examples

```jldoctest
julia> f(x) = 1 - sum(abs2, x); # unit circle

julia> BG = bisect(f, (0.0:1.0, 0.0:1.0); iterations=3);

julia> interpolate(linearroot(f), BG)
8-element Vector{Tuple{Float64, Float64}}:
 (1.0, 0.0)
 (0.9642857142857143, 0.25)
 (0.8571428571428571, 0.5)
 (0.75, 0.65)
 (0.0, 1.0)
 (0.25, 0.9642857142857143)
 (0.65, 0.75)
 (0.5, 0.8571428571428571)

julia> interpolate(BG)
8-element Vector{Tuple{Float64, Float64}}:
 (1.0, 0.0)
 (0.9642857142857143, 0.25)
 (0.8571428571428571, 0.5)
 (0.75, 0.65)
 (0.0, 1.0)
 (0.25, 0.9642857142857143)
 (0.65, 0.75)
 (0.5, 0.8571428571428571)

"""
function interpolate(BG::BisectionGrid)
    # Map from points in the domain to evaluations
    _getindex = Base.Fix1(domainindex, domain(BG))
    IV = zip(BG.evalinds, BG.evaluations)
    f = Dict(_getindex(I) => val for (I, val) in IV)

    # Find roots using linear interpolation
    roots = map(edges(BG)) do edge
        x = edgebounds(edge)
        b = (f[edge[2]] - f[edge[1]])/(x[2] - x[1])
        a = f[edge[1]] - b*x[1]
        return edgetuple(edge)(-a/b)
    end

    return roots
end
