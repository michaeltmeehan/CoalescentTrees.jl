@testset "Tree Constructor Tests" begin
    tree = CoalescentTrees.Tree()
    @test length(tree) == 0
    @test CoalescentTrees.in_tree(tree, 1) == false
    @test tree.time == []
    @test tree.left == []
    @test tree.right == []
    @test CoalescentTrees.get_leaves(tree) == []
    @test_throws AssertionError CoalescentTrees.get_children(tree, 1)
    @test_throws AssertionError CoalescentTrees.get_parent(tree, 1)
    @test_throws AssertionError CoalescentTrees.get_ancestors(tree, 1)
    @test_throws AssertionError CoalescentTrees.get_descendants(tree, 1)

    time = [-0.33, -0.17, -0.15, 0., 0.5, 0.6, 1., 1.5, 2.]
    left = [2, 5, 8, 0, 0, 7, 0, 0, 0]
    right = [4, 3, 6, 0, 0, 9, 0, 0, 0]
    tree = CoalescentTrees.Tree(time, left, right)
    @test length(tree) == 9
    @test CoalescentTrees.in_tree(tree, 1) == true
    @test tree.time == time
    @test tree.left == left
    @test tree.right == right
end


@testset "Tree Equality and Hashing" begin
    time = [-0.33, -0.17, -0.15, 0., 0.5, 0.6, 1., 1.5, 2.]
    left = [2, 5, 8, 0, 0, 7, 0, 0, 0]
    right = [4, 3, 6, 0, 0, 9, 0, 0, 0]
    tree1 = CoalescentTrees.Tree(time, left, right)
    tree2 = CoalescentTrees.Tree(time, left, right)
    tree3 = CoalescentTrees.Tree(time, left, [right[end]; right[1:end-1]])

    @test tree1 == tree2
    @test hash(tree1) == hash(tree2)
    @test tree1 != tree3
end


@testset "Tree Indexing and Slicing" begin
    time = [-0.33, -0.17, -0.15, 0., 0.5, 0.6, 1., 1.5, 2.]
    left = [2, 5, 8, 0, 0, 7, 0, 0, 0]
    right = [4, 3, 6, 0, 0, 9, 0, 0, 0]
    tree = CoalescentTrees.Tree(time, left, right)

    node = tree[2]
    @test node.time == time[2]
    @test node.left == left[2]
    @test node.right == right[2]

    subtree = tree[1:3]
    @test isa(subtree, Tree)
    @test length(subtree) == 3
    @test subtree.time == time[1:3]
    @test issorted(subtree.time)
    @test all(subtree.left[i] ≥ i for i in CoalescentTrees.eachnode(subtree) if subtree.left[i] != 0)
    @test all(subtree.right[i] ≥ i for i in CoalescentTrees.eachnode(subtree) if subtree.right[i] != 0)

end


@testset "Tree Structure Utilities" begin
    time = [-0.33, -0.17, -0.15, 0., 0.5, 0.6, 1., 1.5, 2.]
    left = [2, 5, 8, 0, 0, 7, 0, 0, 0]
    right = [4, 3, 6, 0, 0, 9, 0, 0, 0]
    tree = CoalescentTrees.Tree(time, left, right)

    @test CoalescentTrees.is_leaf(tree, 1) == false
    @test CoalescentTrees.is_leaf(tree, 2) == false
    @test CoalescentTrees.is_leaf(tree, 4) == true
    @test CoalescentTrees.is_internal(tree, 1) == true
    @test CoalescentTrees.is_internal(tree, 2) == true
    @test CoalescentTrees.is_internal(tree, 4) == false
    @test CoalescentTrees.is_binary(tree, 1) == true
    @test CoalescentTrees.is_binary(tree, 2) == true
    @test CoalescentTrees.is_binary(tree, 4) == false
    @test CoalescentTrees.get_leaves(tree) == [4, 5, 7, 8, 9]
    @test CoalescentTrees.get_children(tree, 1) == [2, 4]
    @test CoalescentTrees.get_children(tree, 2) == [5, 3]
    @test CoalescentTrees.get_children(tree, 4) == []
    @test isnothing(CoalescentTrees.get_parent(tree, 1))
    @test CoalescentTrees.get_parent(tree, 2) == 1
    @test CoalescentTrees.get_parent(tree, 3) == 2
    @test CoalescentTrees.get_parent(tree, 4) == 1
    @test CoalescentTrees.get_parent(tree, 5) == 2
    @test CoalescentTrees.get_parent(tree, 6) == 3
    @test CoalescentTrees.get_parent(tree, 7) == 6
    @test CoalescentTrees.get_parent(tree, 8) == 3
    @test CoalescentTrees.get_parent(tree, 9) == 6
    @test CoalescentTrees.get_ancestors(tree, 1) == []
    @test CoalescentTrees.get_ancestors(tree, 2) == [1]
    @test CoalescentTrees.get_ancestors(tree, 3) == [2, 1]
    @test CoalescentTrees.get_ancestors(tree, 4) == [1]
    @test CoalescentTrees.get_ancestors(tree, 5) == [2, 1]
    @test CoalescentTrees.get_ancestors(tree, 6) == [3, 2, 1]

    subtree = CoalescentTrees.get_subtree(tree, 2)
    @test subtree isa Tree
    @test subtree.time[1] == tree.time[2]
    @test length(subtree) == length(CoalescentTrees.get_descendants(tree, 2)) + 1

end