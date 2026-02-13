module MultiBisect

using Base.Threads: @threads, @spawn, fetch
using Logging
using Roots

import Base.show

export BisectionGrid
export domain, evaluated, evaluations, efficiency, splitsign
export bisect
export edges, interpolate, edgeroot, linearroot, marchingsquares
export ispositivesign, isnegativesign, iszerosign, isevaluated, differentsigns

include("bisection.jl")
include("interpolation.jl")

end
