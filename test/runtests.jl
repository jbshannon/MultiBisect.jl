using MultiBisect
using Test
using MultiBisect: forwardinds, domainindex, edgedim, edgebounds

@testset "MultiBisect.jl" begin
    @testset "BisectionGrid" begin
        # Test 1D grid construction
        grid1d = (0.0:1.0,)
        BG1 = BisectionGrid(grid1d)
        @test BG1 isa BisectionGrid{Float64, 1}
        @test length(BG1.signs) == 2
        @test size(BG1.signs) == (2,)
        @test size(BG1.evaluated) == (2,)
        @test size(BG1.toevaluate) == (2,)
        @test domain(BG1) == (0.0:1.0:1.0,)
        
        # Test 2D grid construction
        grid2d = (0.0:1.0, 0.0:0.5:1.0)
        BG2 = BisectionGrid(grid2d)
        @test BG2 isa BisectionGrid{Float64, 2}
        @test size(BG2.signs) == (2, 3)
        @test domain(BG2) == (0.0:1.0:1.0, 0.0:0.5:1.0)
        
        # Test property functions
        @test evaluations(BG2) == 0  # No evaluations initially
        @test efficiency(BG2) == 1.0  # 100% efficient initially (no evals)
        @test evaluated(BG2) == BG2.evaluated
    end
    
    @testset "Bisection Algorithm" begin
        # Test unit circle example from README
        f(x) = 1 - sum(abs2, x)  # unit circle: x² + y² = 1
        grid = (0.0:1.0, 0.0:1.0)
        
        # Basic bisect test
        BG = bisect(f, grid; iterations=3)
        @test BG isa BisectionGrid{Float64, 2}
        @test evaluations(BG) > 0
        @test evaluations(BG) < length(BG.signs)  # Should be efficient
        
        # Test monotonic option
        BG_mono = bisect(f, grid; iterations=3, monotonic=true)
        @test evaluations(BG_mono) <= evaluations(BG)  # Should use fewer evaluations
        
        # Test different iteration counts
        BG2 = bisect(f, grid; iterations=2)
        BG3 = bisect(f, grid; iterations=3)
        @test length(BG3.signs) > length(BG2.signs)  # More iterations = finer grid
        
        # Test verbose mode doesn't error
        BG_verbose = bisect(f, grid; iterations=2, verbose=true)
        @test evaluations(BG_verbose) > 0
    end
    
    @testset "Utility Functions" begin
        f(x) = 1 - sum(abs2, x)
        BG = bisect(f, (0.0:1.0, 0.0:1.0); iterations=3)
        
        # Test splitsign
        posx, negx = splitsign(BG)
        # Count zero values directly from signs array
        zero_count = count(sign -> sign == (false, false), BG.signs[BG.evaluated])
        pos_count = length(posx)
        neg_count = length(negx)
        @test pos_count + neg_count + zero_count == evaluations(BG)
        @test all(x -> f(x) > 0, posx)  # Strictly positive now
        @test all(x -> f(x) < 0, negx)  # Strictly negative now
        
        # Test edges
        edge_list = edges(BG)
        @test length(edge_list) > 0
        @test all(edge -> length(edge) == 2, edge_list)  # Each edge has 2 endpoints
        @test all(edge -> length(edge[1]) == 2, edge_list)  # 2D points
        
        # Test efficiency calculation
        eff = efficiency(BG)
        @test 0 <= eff <= 1
    end
    
    @testset "Interpolation Methods" begin
        f(x) = 1 - sum(abs2, x)
        BG = bisect(f, (0.0:1.0, 0.0:1.0); iterations=3)
        edge_list = edges(BG)
        
        # Test marchingsquares
        ms_roots = interpolate(marchingsquares, BG)
        @test length(ms_roots) == length(edge_list)
        @test all(root -> length(root) == 2, ms_roots)
        
        # Test linearroot
        linear_roots = interpolate(linearroot(f), BG)
        @test length(linear_roots) == length(edge_list)
        @test all(root -> length(root) == 2, linear_roots)
        
        # Test edgeroot with Roots
        using Roots
        root_roots = interpolate(edge -> edgeroot(f, edge, Roots.ITP()), BG)
        @test length(root_roots) == length(edge_list)
        @test all(root -> length(root) == 2, root_roots)
        
        # Test interpolate without rootfinder (uses saved evaluations)
        saved_roots = interpolate(BG)
        @test length(saved_roots) == length(edge_list)
        @test all(root -> length(root) == 2, saved_roots)
        
        # Test that different methods give different results
        # Use elementwise comparison for tuples
        differences = [!isapprox(ms_roots[i][j], linear_roots[i][j]) for i in 1:length(ms_roots) for j in 1:length(ms_roots[i])]
        @test any(differences)
    end
    
    @testset "Higher Dimensions" begin
        f(x) = 1 - sum(abs2, x)  # unit hypersphere
        
        # Test 3D
        grid3d = (0.0:1.0, 0.0:1.0, 0.0:1.0)
        BG3 = bisect(f, grid3d; iterations=2)
        @test BG3 isa BisectionGrid{Float64, 3}
        @test evaluations(BG3) > 0
        
        # Test interpolation in 3D
        roots3 = interpolate(marchingsquares, BG3)
        @test all(root -> length(root) == 3, roots3)
        
        # Test 5D example from README
        grid5d = ntuple(i -> (0.0:1.0), 5)
        BG5 = bisect(f, grid5d; iterations=1)  # Fewer iterations for performance
        @test BG5 isa BisectionGrid{Float64, 5}
        @test evaluations(BG5) > 0
        
        roots5 = interpolate(marchingsquares, BG5)
        @test all(root -> length(root) == 5, roots5)
    end
    
    @testset "Edge Cases" begin
        # Test linear function
        f_linear(x) = x[1] - x[2]  # f(x,y) = x - y
        BG = bisect(f_linear, (0.0:1.0, 0.0:1.0); iterations=3)
        @test evaluations(BG) > 0
        
        # Test function with no sign changes
        f_positive(x) = 1 + sum(abs2, x)  # Always positive
        BG_pos = bisect(f_positive, (0.0:1.0, 0.0:1.0); iterations=2)
        @test evaluations(BG_pos) > 0
        
        # Test with different number types (converted to floats by BisectionGrid)
        BG_int = bisect(x -> 1 - sum(abs2, x), (0:1.0, 0:1.0); iterations=2)
        @test BG_int isa BisectionGrid{Float64, 2}
    end
    
    @testset "Helper Functions" begin
        # Test forwardinds
        @test forwardinds(1) == (CartesianIndex(1),)
        @test forwardinds(2) == (CartesianIndex(1, 0), CartesianIndex(0, 1))
        @test forwardinds(3) == (CartesianIndex(1, 0, 0), CartesianIndex(0, 1, 0), CartesianIndex(0, 0, 1))
        
        # Test domainindex
        X = (1.0:3.0, 0.5:1.5)
        @test domainindex(X, CartesianIndex(1, 1)) == (1.0, 0.5)
        @test domainindex(X, CartesianIndex(2, 1)) == (2.0, 0.5)
        @test domainindex(X, CartesianIndex(1, 2)) == (1.0, 1.5)
        
        # Test edgedim
        edge1 = ((0.5, 1.0), (1.0, 1.0))
        edge2 = ((0.5, 1.0), (0.5, 1.5))
        @test edgedim(edge1) == 1
        @test edgedim(edge2) == 2
        
        # Test edgebounds
        @test edgebounds(edge1) == (0.5, 1.0)
        @test edgebounds(edge2) == (1.0, 1.5)
    end
end
