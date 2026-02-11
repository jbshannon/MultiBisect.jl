# AGENTS.md - MultiBisect.jl Development Guide

This guide provides essential information for agentic coding agents working on the MultiBisect.jl package.

## Quick Reference

### Essential Commands
```bash
# Install and precompile dependencies
julia --project -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"

# Run all tests
julia --project -e "using Pkg; Pkg.test()"

# Run specific test (add @testset blocks first)
julia --project -e "using Pkg; Pkg.test(test_args=\"test_name\")"

# Build package after changes
julia --project -e "using Pkg; Pkg.build()"

# Generate documentation
julia --project=docs docs/make.jl"

# Check Julia version compatibility
julia --project -e "using Pkg; Pkg.status()"
```

### Development Workflow
1. Always run `Pkg.test()` after significant changes
2. Use `julia --project` for package environment
3. Include doctests in all new functions
4. Update `Project.toml` version when adding breaking changes

## Code Style Guidelines

### Module Structure
```julia
# Main module (src/MultiBisect.jl)
module MultiBisect

# Imports - grouped by purpose
using Base.Threads: @threads, @spawn, fetch
using Logging
using Roots

# Selective imports for extending functionality
import Base.show

# Export list - alphabetical order
export BisectionGrid, MultiBisectMethod, bisect, evaluate, findroots

# Include submodules
include("bisection.jl")
include("interpolation.jl")

end # module MultiBisect
```

### Type Definitions
```julia
# Parametric types with constraints
struct BisectionGrid{T <: Number, N}
    data::Array{T, N}
    ranges::NTuple{N, AbstractRange{T}}
end

# Abstract types for extensibility
abstract type AbstractBisectionMethod end

# Concrete implementations
struct MultiBisectMethod <: AbstractBisectionMethod
    # fields...
end
```

### Function Signatures
```julia
# Public functions with full type annotations
function bisect(
    f::Function, 
    grid::BisectionGrid{T, N};
    threaded::Bool=false, 
    monotonic::Bool=false,
    kwargs...
) where {T <: Number, N}
    # implementation
end

# Private functions (prefixed with _)
function _evaluate!(result::AbstractArray, f::Function, grid::BisectionGrid)
    # mutation function uses !
end
```

### Naming Conventions
- **Types**: PascalCase (`BisectionGrid`, `MultiBisectMethod`)
- **Functions**: snake_case (`splitsign`, `fillnonbracketing!`)
- **Variables**: 
  - camelCase for coordinates/domain points (`posx`, `negx`)
  - snake_case for general variables (`grid_data`, `method_config`)
- **Constants**: UPPER_SNAKE_CASE
- **Private functions**: prefix with `_`

### Import Organization
```julia
# Standard library imports first
using Base.Threads: @threads, @spawn, fetch
using Logging

# External packages next
using Roots

# Selective imports for extending
import Base.show
import Base: getindex, setindex!
```

### Error Handling
```julia
# Use Julia's built-in error system
function risky_operation(x)
    x < 0 && throw(ArgumentError("x must be non-negative"))
    # implementation
end

# Boundary checks for array operations
function safe_index(arr, i)
    i < 1 || i > length(arr) && 
        throw(BoundsError(arr, i))
    return arr[i]
end

# Zero checks for numerical operations
iszero(result) && return zero(result)
```

### Documentation Standards
```julia
"""
    bisect(f, grid; threaded=false, monotonic=false)

Find roots of function `f` on the given `grid` using multi-dimensional bisection.

# Arguments
- `f::Function`: Function to find roots of
- `grid::BisectionGrid`: Grid defining the search domain

# Keywords
- `threaded::Bool=false`: Enable multithreading
- `monotonic::Bool=false`: Assume monotonic function

# Returns
- `Vector{T}`: Found root locations

# Examples
```jldoctest
julia> grid = BisectionGrid((-1.0:0.1:1.0, -1.0:0.1:1.0));
julia> roots = bisect(x -> x^2 + y^2 - 1, grid);
```

# Extended help
The algorithm implements a multi-dimensional extension of the classic bisection method...
"""
```

### Testing Patterns
```julia
# In test/runtests.jl
using Test
using MultiBisect

@testset "MultiBisect.jl" begin
    @testset "BisectionGrid" begin
        grid = BisectionGrid((-1.0:0.1:1.0,))
        @test eltype(grid) == Float64
        @test ndims(grid) == 1
    end
    
    @testset "Root finding" begin
        f(x) = x^2 - 1
        grid = BisectionGrid((-2.0:0.1:2.0,))
        roots = bisect(f, grid)
        @test length(roots) == 2
    end
end
```

## Performance Guidelines

### Memory Management
- Use mutation functions (`!`) for large arrays
- Preallocate arrays when possible
- Leverage `CartesianIndices` for efficient multi-dimensional access

### Threading
- Use `@threads` for embarrassingly parallel loops
- Use `@spawn`/`fetch` for task-based parallelism
- Always test both threaded and non-threaded paths

### Type Stability
- Ensure function return types are predictable
- Use `where` clauses for parametric functions
- Avoid `Any` type annotations in hot paths

## Version Compatibility

### Supported Julia Versions
- Minimum: 1.6
- Recommended: 1.9+
- Tested: 1.6, 1.9, nightly (CI)

### Dependencies
- `Roots`: v2+ for root finding utilities
- `Logging`: 1.6+ (standard library)
- Development: `Test` (standard library)
- Documentation: `Documenter` (docs only)

## Code Quality Checklist

Before submitting changes:
- [ ] All tests pass (`Pkg.test()`)
- [ ] New functions have doctests
- [ ] Code follows naming conventions
- [ ] Types are properly annotated
- [ ] Error handling is appropriate
- [ ] Documentation is updated
- [ ] Version number updated if breaking changes