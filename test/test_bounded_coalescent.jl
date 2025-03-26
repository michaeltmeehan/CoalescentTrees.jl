
@testset "Coalescent Probability" begin
    Ne = 100.0

    # --- Valid input tests ---

    p = CoalescentTrees.coalescent_probability(5, 1, 1.0, Ne)
    @test isapprox(p, 7.206712e-8; atol=1e-14)

    p = CoalescentTrees.coalescent_probability(5, 2, 1.0, Ne)
    @test isapprox(p, 2.854016e-5; atol=1e-5)

    p = CoalescentTrees.coalescent_probability(3, 3, 0.0, Ne)
    @test isapprox(p, 1.0; atol=1e-10)

    p = CoalescentTrees.coalescent_probability(3, 3, 0.5, Ne)
    @test isapprox(p, 0.9851119; atol=1e-6)

    p = CoalescentTrees.coalescent_probability(5, 1, 1000.0, Ne)
    @test isapprox(p, 0.9999092; atol=1e-6)

    p = CoalescentTrees.coalescent_probability(5, 1, 1e-10, Ne)
    @test isapprox(p, 0.0; atol=1e-6)

    # --- Invalid input tests ---

    @test_throws AssertionError CoalescentTrees.coalescent_probability(3, 4, 1.0, Ne)  # n_small > n_big
    @test_throws AssertionError CoalescentTrees.coalescent_probability(0, 1, 1.0, Ne)  # n_big < 1
    @test_throws AssertionError CoalescentTrees.coalescent_probability(1, 0, 1.0, Ne)  # n_small < 1
    @test_throws AssertionError CoalescentTrees.coalescent_probability(5, 2, -1.0, Ne) # Δt < 0
    @test_throws AssertionError CoalescentTrees.coalescent_probability(5, 2, 1.0, 0.0) # Ne <= 0

end


@testset "Forward and Backwards Probability Calculator" begin
    sampled_sequences = [1, 3, 2]
    sequence_times = [1.0, 2.0, 3.0]
    bound_time = 0.0
    Ne = 10.0

    forward_probs = CoalescentTrees.calc_forward_probs(sampled_sequences, sequence_times, bound_time, Ne)

    @test size(forward_probs) == (sum(sampled_sequences), length(sequence_times) + 1)
    @test all(0. .≤ forward_probs .≤ 1.)
    @test all(sum(forward_probs, dims=1) .≈ 1)

    lineages, backwards_probs = CoalescentTrees.sample_lineages_backward(forward_probs, sampled_sequences, sequence_times, bound_time, Ne)

    @test size(backwards_probs) == size(forward_probs)
    @test all(lineages .≥ 1)
    @test lineages[1] == 1
    @test all(0. .≤ backwards_probs .≤ 1.)
    @test all(sum(backwards_probs, dims=1) .≈ 1)

    max_lineages = [sum(sampled_sequences); CoalescentTrees.reverse_cumsum(sampled_sequences)]
    @test lineages[end] == max_lineages[end]
    @test all(lineages .≤ max_lineages)
    @test all(lineages .≥ [0; sampled_sequences])

end


@testset "Partition Intervals" begin
    sampled_sequences = [2, 3]
    sequence_times = [2.0, 3.0]
    bound_time = 0.0
    Ne = 10.0

    # Use forward-backward to get input lineages
    fp = CoalescentTrees.calc_forward_probs(sampled_sequences, sequence_times, bound_time, Ne)
    lineages, _ = CoalescentTrees.sample_lineages_backward(fp, sampled_sequences, sequence_times, bound_time, Ne)

    times, new_lineages, new_samples = CoalescentTrees.partition_intervals(lineages, sampled_sequences, sequence_times, bound_time, Ne)

    @test issorted(times)
    @test length(times) == length(new_lineages) == length(new_samples)

    # Check that every interval has ≤ 1 coalescent event
    for i in 1:(length(times)-1)
        n_start = new_lineages[i] - new_samples[i]
        n_end   = new_lineages[i+1]
        @test n_end - n_start ≤ 1
    end

    # Check the lineage counts at the start and end of the intervals
    for i in 1:(length(times) - 1)
        @test new_lineages[i+1] ≥ new_lineages[i] - new_samples[i]
    end

    # Check that sample points were preserved
    @test sum(new_samples) == sum(sampled_sequences)
end


using Test, Statistics

@testset "Sample Coalescent Time" begin
    t_lower = 0.0
    t_upper = 1.0
    Ne = 10.0

    for n_small in [2, 5, 10]
        times = [CoalescentTrees.sample_coalescent_time(t_lower, t_upper, n_small, Ne) for _ in 1:10_000]
        @test all(t_lower .≤ times .≤ t_upper)
    end

    # histogram([CoalescentTrees.sample_coalescent_time(0.0, 1.0, 5, 10.0) for _ in 1:10_000], bins=50)

    # Degenerate interval
    t_fixed = CoalescentTrees.sample_coalescent_time(2.5, 2.5, 5, Ne)
    @test t_fixed == 2.5
end


@testset "Construct Event List" begin
    sampled_sequences = [1, 2]
    sequence_times = [1.0, 2.0]
    bound_time = 0.0
    Ne = 10.0

    # Get lineages from forward-backward
    fp = CoalescentTrees.calc_forward_probs(sampled_sequences, sequence_times, bound_time, Ne)
    lineages, _ = CoalescentTrees.sample_lineages_backward(fp, sampled_sequences, sequence_times, bound_time, Ne)
    interval_times, lineages, samples = CoalescentTrees.partition_intervals(lineages, sampled_sequences, sequence_times, bound_time, Ne)

    events = CoalescentTrees.construct_event_list(interval_times, lineages, samples, Ne)

    # Check event structure
    @test all(e -> e.type in (0, 1, 2), events)
    @test any(e -> e.type == 2, events)  # Should include at least one coalescent
    @test count(e -> e.type == 1, events) == 1  # Exactly one root
    @test count(e -> e.type == 0, events) == sum(sampled_sequences)

    # Check times are within bounds
    @test all(e -> e.time ≥ minimum(interval_times) && e.time ≤ maximum(interval_times), events)

    # Optional: ensure no duplicate sampling times unless samples > 1
    sampling_times = [e.time for e in events if e.type == 0]
    @test length(sampling_times) == sum(samples)
end


@testset "Sample Topology" begin
    sampled_sequences = [1, 3, 2]
    sequence_times = [1.0, 2.0, 3.0]
    bound_time = 0.0
    Ne = 10.0

    fp = CoalescentTrees.calc_forward_probs(sampled_sequences, sequence_times, bound_time, Ne)
    lineages, _ = CoalescentTrees.sample_lineages_backward(fp, sampled_sequences, sequence_times, bound_time, Ne)
    interval_times, lineages, samples = CoalescentTrees.partition_intervals(lineages, sampled_sequences, sequence_times, bound_time, Ne)
    events = CoalescentTrees.construct_event_list(interval_times, lineages, samples, Ne)
    tree = CoalescentTrees.sample_topology(events)

    n_leaves = length(get_leaves(tree))
    @test n_leaves == sum(sampled_sequences)
    @test !CoalescentTrees.is_binary(tree, 1)
    @test length(tree.time) == length(tree.left) == length(tree.right) == 2 * n_leaves

    @test all(minimum(interval_times) .≤ tree.time .≤ maximum(interval_times))
    # Internal nodes should reference valid indices or 0
    @test all(tree.left[i] ≥ i for i in CoalescentTrees.eachnode(tree) if tree.left[i] != 0)
    @test all(tree.right[i] ≥ i for i in CoalescentTrees.eachnode(tree) if tree.right[i] != 0)

    # Check every node has one and only one parent (except root)
    for i in CoalescentTrees.eachnode(tree)
        if i != 1
            parent = CoalescentTrees.get_parent(tree, i)
            @test !isnothing(parent)
        end
    end
    # No duplicate children
    @test all(tree[node].left .!= tree[node].right for node in CoalescentTrees.eachnode(tree) if !CoalescentTrees.is_leaf(tree, node))
end