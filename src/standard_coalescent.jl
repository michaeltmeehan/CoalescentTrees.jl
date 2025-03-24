import Distributions: Exponential


function sample_tree(sampled_sequences::Vector{Int}, sequence_times::Vector{Float64}, Ne::Float64)
    # Assertions to validate input arguments
    @assert length(sampled_sequences) == length(sequence_times) "sampled_sequences and sequence_times must have the same length"
    @assert all(x -> x >= 0, sampled_sequences) "sampled_sequences must contain non-negative integers"
    @assert issorted(sequence_times) "sequence_times must be sorted in descending order"
    @assert Ne > 0 "Effective population size (Ne) must be positive"

    total_sequences = sum(sampled_sequences)
    total_nodes = 2 * total_sequences - 1
    time = zeros(Float64, total_nodes)
    left = zeros(Int, total_nodes)
    right = zeros(Int, total_nodes)

    active_nodes = Int[]
    time_steps = [-Inf; sequence_times]
    new_sequences = [0; sampled_sequences]
    node_idx = total_nodes
    for t_idx in length(time_steps):-1:2
        current_time = time_steps[t_idx]
        # Add new sampled lineages to tree
        for _ in 1:new_sequences[t_idx]
            time[node_idx] = current_time
            insert!(active_nodes, rand(1:length(active_nodes)+1), node_idx)
            node_idx -= 1
        end

        # Coalesce lineages
        next_time = time_steps[t_idx-1]
        n_active = length(active_nodes)
        time_to_coalescence = rand(Exponential(2 * Ne / (n_active * (n_active - 1))))
        while current_time - time_to_coalescence > next_time # && active_nodes > 1
            current_time -= time_to_coalescence
            time[node_idx] = current_time
            left_child, right_child = pop!(active_nodes), pop!(active_nodes)
            left[node_idx], right[node_idx] = left_child, right_child
            insert!(active_nodes, rand(1:length(active_nodes)+1), node_idx)
            node_idx -= 1
            n_active -= 1
            n_active < 2 && break
            time_to_coalescence = rand(Exponential(2 * Ne / (n_active * (n_active - 1))))
        end
    end
    return Tree(time, left, right)
end