module CoalescentTrees

using DataFrames
using Plots
using RecipesBase
using StatsBase

include("tree_representation.jl")
include("standard_coalescent.jl")
include("bounded_coalescent.jl")

export sample_tree, Tree, eachnode, get_leaves, is_leaf, get_parent, get_children, is_internal, is_binary, in_tree, get_ancestors
export get_descendants, get_subtree
export coalescent_probability

end # module CoalescentTrees
