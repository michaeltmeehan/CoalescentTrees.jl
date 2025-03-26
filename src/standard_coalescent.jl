###############################################################################
# Coalescent Tree Sampler: Standard (Unbounded) Coalescent Model
#
# This script implements a simulation of genealogical trees under the standard
# (Kingman) coalescent process with a constant effective population size (Ne).
# Unlike the bounded coalescent model, this version assumes that the coalescent
# process continues indefinitely backward in time, without conditioning on a
# fixed bound or root time.
#
# -----------------------------------------------------------------------------
# Algorithm Overview:
#
# 1. Traverse backward in time across sampling time points.
# 2. At each time point:
#    a. Add newly sampled lineages to the active pool.
#    b. Simulate exponential waiting times between coalescent events:
#       - The rate is determined by the number of active lineages: λ = k(k−1)/(2Ne).
#       - Continue simulating coalescent events while the next event occurs before
#         the previous sampling time.
# 3. For each coalescence:
#    - Sample two lineages at random to coalesce.
#    - Create an internal node with its time and descendant pointers.
#
# -----------------------------------------------------------------------------
# Inputs:
# - `sampled_sequences::Vector{Int}`: number of sequences sampled at each time.
# - `sequence_times::Vector{Float64}`: corresponding sampling times (ascending order).
# - `Ne::Float64`: effective population size.
#
# Output:
# - `Tree`: a binary tree with fields:
#     - `time::Vector{Float64}`: coalescent/sampling times for all nodes.
#     - `left::Vector{Int}`, `right::Vector{Int}`: indices of child nodes.
#
# Notes:
# - Output tree includes both leaf and internal nodes, with total size 2n − 1.
# - Sampling events and coalescences are processed backward in time.
# - The root is the final coalescent event joining the last two lineages.
#
###############################################################################


import Distributions: Exponential


function sample_tree(sampled_sequences::Vector{Int}, sequence_times::Vector{Float64}, Ne::Float64)
    # Assertions to validate input arguments
    @assert length(sampled_sequences) == length(sequence_times) "sampled_sequences and sequence_times must have the same length"
    @assert all(x -> x >= 0, sampled_sequences) "sampled_sequences must contain non-negative integers"
    @assert issorted(sequence_times) "sequence_times must be sorted in ascending order"
    @assert Ne > 0 "Effective population size (Ne) must be positive"

    total_sequences = sum(sampled_sequences)
    total_nodes = 2 * total_sequences - 1
    time = zeros(Float64, total_nodes)
    left = zeros(Int, total_nodes)
    right = zeros(Int, total_nodes)

    active_nodes = Int[]
    time_steps = [-Inf; sequence_times]
    new_sequences = [0; sampled_sequences]
    node = total_nodes
    for t_idx in length(time_steps):-1:2
        current_time = time_steps[t_idx]
        # Add new sampled lineages to tree
        for _ in 1:new_sequences[t_idx]
            time[node] = current_time
            push!(active_nodes, node)
            node -= 1
        end

        # Coalesce lineages
        next_time = time_steps[t_idx-1]
        n_active = length(active_nodes)
        time_to_coalescence = rand(Exponential(2 * Ne / (n_active * (n_active - 1))))
        while current_time - time_to_coalescence > next_time # && active_nodes > 1
            current_time -= time_to_coalescence
            time[node] = current_time
            left_child, right_child = pop_random!(active_nodes), pop_random!(active_nodes)
            left[node], right[node] = left_child, right_child
            push!(active_nodes, node)
            node -= 1
            n_active -= 1
            n_active < 2 && break
            time_to_coalescence = rand(Exponential(2 * Ne / (n_active * (n_active - 1))))
        end
    end
    return Tree(time, left, right)
end