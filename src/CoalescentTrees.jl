module CoalescentTrees

using DataFrames
using Plots
using RecipesBase
using StatsBase

include("utils.jl")
include("tree_representation.jl")
include("standard_coalescent.jl")
include("bounded_coalescent.jl")

export pop_random!, reverse_cumsum
export sample_tree, Tree, eachnode, get_leaves, is_leaf, get_parent, get_children, is_internal, is_binary, in_tree, get_ancestors
export get_descendants, get_subtree
export coalescent_probability, calc_forward_probs, calc_backward_probs, partition_intervals, sample_coalescence_time, sample_event_times, sample_topology, reverse_cumsum

end # module CoalescentTrees
