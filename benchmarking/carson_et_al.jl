using CoalescentTrees
using Plots

sampled_sequences = fill(1, 5)
sequence_times = [0.0, 0.5, 1.0, 1.5, 2.0]
bound_time = -0.5
Ne = 1.0

# Empirical estimate for bound probability
n_trees = 500_000
accepted = 0
for _ in 1:n_trees
    tree = sample_tree(sampled_sequences, sequence_times, Ne)
    if tree.time[1] ≥ bound_time
        accepted += 1
    end
end

println("Empirical probability: ", accepted / n_trees) # Compare with published value of 0.234



# Distribution of coalescent times
function get_coalescent_times(tree::Tree)
    coalescent_times = Float64[]
    for i in eachnode(tree)
        if is_binary(tree, i)
            push!(coalescent_times, tree.time[i])
        end
    end
    return coalescent_times
end

n_trees = 500_000
bounded_coalescent_times = []
unbounded_coalescent_times = []
for _ in 1:n_trees
    bounded_tree = sample_tree(sampled_sequences, sequence_times, bound_time, Ne)
    push!(bounded_coalescent_times, get_coalescent_times(bounded_tree))
    unbounded_tree = sample_tree(sampled_sequences, sequence_times, Ne)
    if unbounded_tree.time[1] ≥ bound_time
        push!(unbounded_coalescent_times, get_coalescent_times(unbounded_tree))
    end
end

plots = []
for i in 1:4
    p = histogram(getindex.(bounded_coalescent_times, i), bins=50, normalize=:pdf, alpha=0.5, label="Bounded")
    histogram!(getindex.(unbounded_coalescent_times, i), bins=50, normalize=:pdf, alpha=0.5, label="Unbounded")
    push!(plots, p)
end
plot(plots..., layout=(2,2))