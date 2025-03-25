@testset "Unbounded Coalescent Tree test" begin
    # --- Valid input test ---
    Ne = 1.0
    sampled_sequences = [2, 1]
    sequence_times = [1.0, 3.0]
    tree = sample_tree(sampled_sequences, sequence_times, Ne)

    @test tree isa Tree
    @test length(tree) == 2 * sum(sampled_sequences) - 1
    @test all(isfinite.(tree.time))
    @test issorted(tree.time)
    @test all(tree.left[i] ≥ i for i in eachnode(tree) if tree.left[i] != 0)
    @test all(tree.right[i] ≥ i for i in eachnode(tree) if tree.right[i] != 0)

    # --- Assertion checks ---
    @test_throws AssertionError sample_tree(sampled_sequences, [3.], Ne)
    @test_throws AssertionError sample_tree([-1, 1], sequence_times, Ne)
    @test_throws AssertionError sample_tree(sampled_sequences, reverse(sequence_times), 0.0)
    @test_throws AssertionError sample_tree(sampled_sequences, sequence_times, -1.0)

    # --- Edge case: single lineage ---
    sampled_sequences = [1]
    sequence_times = [0.0]
    tree = sample_tree(sampled_sequences, sequence_times, Ne)
    @test length(tree) == 1
    @test tree.time[1] == 0.0
    @test is_leaf(tree, 1)

end