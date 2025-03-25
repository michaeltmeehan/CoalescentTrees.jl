using Test
using CoalescentTrees


@testset "Tree Constructor Tests" begin
    tree = Tree()
    @test length(tree) == 0
    @test in_tree(tree, 1) == false
    @test tree.time == []
    @test tree.left == []
    @test tree.right == []
    @test get_leaves(tree) == []
    @test_throws AssertionError get_children(tree, 1)
    @test_throws AssertionError get_parent(tree, 1)
    @test_throws AssertionError get_ancestors(tree, 1)
    @test_throws AssertionError get_descendants(tree, 1)

    time = [-0.33, -0.17, -0.15, 0., 0.5, 0.6, 1., 1.5, 2.]
    left = [2, 5, 8, 0, 0, 7, 0, 0, 0]
    right = [4, 3, 6, 0, 0, 9, 0, 0, 0]
    tree = Tree(time, left, right)
    @test length(tree) == 9
    @test in_tree(tree, 1) == true
    @test tree.time == time
    @test tree.left == left
    @test tree.right == right
end


@testset "Tree Equality and Hashing" begin
    time = [-0.33, -0.17, -0.15, 0., 0.5, 0.6, 1., 1.5, 2.]
    left = [2, 5, 8, 0, 0, 7, 0, 0, 0]
    right = [4, 3, 6, 0, 0, 9, 0, 0, 0]
    tree1 = Tree(time, left, right)
    tree2 = Tree(time, left, right)
    tree3 = Tree(time, left, [right[end]; right[1:end-1]])

    @test tree1 == tree2
    @test hash(tree1) == hash(tree2)
    @test tree1 != tree3
end


@testset "Tree Indexing and Slicing" begin
    time = [-0.33, -0.17, -0.15, 0., 0.5, 0.6, 1., 1.5, 2.]
    left = [2, 5, 8, 0, 0, 7, 0, 0, 0]
    right = [4, 3, 6, 0, 0, 9, 0, 0, 0]
    tree = Tree(time, left, right)

    node = tree[2]
    @test node.time == time[2]
    @test node.left == left[2]
    @test node.right == right[2]

    subtree = tree[1:3]
    @test isa(subtree, Tree)
    @test length(subtree) == 3
    @test subtree.time == time[1:3]
    @test issorted(subtree.time)
    @test all(subtree.left[i] ≥ i for i in eachnode(subtree) if subtree.left[i] != 0)
    @test all(subtree.right[i] ≥ i for i in eachnode(subtree) if subtree.right[i] != 0)

end


@testset "Tree Structure Utilities" begin
    time = [-0.33, -0.17, -0.15, 0., 0.5, 0.6, 1., 1.5, 2.]
    left = [2, 5, 8, 0, 0, 7, 0, 0, 0]
    right = [4, 3, 6, 0, 0, 9, 0, 0, 0]
    tree = Tree(time, left, right)

    @test is_leaf(tree, 1) == false
    @test is_leaf(tree, 2) == false
    @test is_leaf(tree, 4) == true
    @test is_internal(tree, 1) == true
    @test is_internal(tree, 2) == true
    @test is_internal(tree, 4) == false
    @test is_binary(tree, 1) == true
    @test is_binary(tree, 2) == true
    @test is_binary(tree, 4) == false
    @test get_leaves(tree) == [4, 5, 7, 8, 9]
    @test get_children(tree, 1) == [2, 4]
    @test get_children(tree, 2) == [5, 3]
    @test get_children(tree, 4) == []
    @test isnothing(get_parent(tree, 1))
    @test get_parent(tree, 2) == 1
    @test get_parent(tree, 3) == 2
    @test get_parent(tree, 4) == 1
    @test get_parent(tree, 5) == 2
    @test get_parent(tree, 6) == 3
    @test get_parent(tree, 7) == 6
    @test get_parent(tree, 8) == 3
    @test get_parent(tree, 9) == 6
    @test get_ancestors(tree, 1) == []
    @test get_ancestors(tree, 2) == [1]
    @test get_ancestors(tree, 3) == [2, 1]
    @test get_ancestors(tree, 4) == [1]
    @test get_ancestors(tree, 5) == [2, 1]
    @test get_ancestors(tree, 6) == [3, 2, 1]

    subtree = get_subtree(tree, 2)
    @test subtree isa Tree
    @test subtree.time[1] == tree.time[2]
    @test length(subtree) == length(get_descendants(tree, 2)) + 1

end


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


@testset "Coalescent Probability" begin
    Ne = 100.0

    # --- Valid input tests ---

    p = coalescent_probability(5, 1, 1.0, Ne)
    @test isapprox(p, 7.206712e-8; atol=1e-14)

    p = coalescent_probability(5, 2, 1.0, Ne)
    @test isapprox(p, 2.854016e-5; atol=1e-5)

    p = coalescent_probability(3, 3, 0.0, Ne)
    @test isapprox(p, 1.0; atol=1e-10)

    p = coalescent_probability(3, 3, 0.5, Ne)
    @test isapprox(p, 0.9851119; atol=1e-6)

    p = coalescent_probability(5, 1, 1000.0, Ne)
    @test isapprox(p, 0.9999092; atol=1e-6)

    p = coalescent_probability(5, 1, 1e-10, Ne)
    @test isapprox(p, 0.0; atol=1e-6)

    # --- Invalid input tests ---

    @test_throws AssertionError coalescent_probability(3, 4, 1.0, Ne)  # n_small > n_big
    @test_throws AssertionError coalescent_probability(0, 1, 1.0, Ne)  # n_big < 1
    @test_throws AssertionError coalescent_probability(1, 0, 1.0, Ne)  # n_small < 1
    @test_throws AssertionError coalescent_probability(5, 2, -1.0, Ne) # Δt < 0
    @test_throws AssertionError coalescent_probability(5, 2, 1.0, 0.0) # Ne <= 0

end