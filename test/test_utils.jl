@testset "reverse_cumsum" begin
    @test CoalescentTrees.reverse_cumsum(Int[]) == Int[]
    @test CoalescentTrees.reverse_cumsum([1]) == [1]
    @test CoalescentTrees.reverse_cumsum([1, 2, 3]) == [6, 5, 3]
    @test CoalescentTrees.reverse_cumsum([0, 0, 0]) == [0, 0, 0]
    @test CoalescentTrees.reverse_cumsum([5, -2, 4]) == [7, 2, 4]
end

@testset "pop_random!" begin
    v = [10, 20, 30]
    popped = CoalescentTrees.pop_random!(v)
    @test popped in (10, 20, 30)
    @test length(v) == 2
    @test sort([popped; v...]) == [10, 20, 30]

    v = [1]
    @test CoalescentTrees.pop_random!(v) == 1
    @test isempty(v)

    v = Int[]
    @test_throws ArgumentError CoalescentTrees.pop_random!(v)
end