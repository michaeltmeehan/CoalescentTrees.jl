import Distributions: Exponential


function sample_tree(sampled_sequences::Vector{Int}, sequence_times::Vector{Float64}, Ne::Float64)
    total_sequences = sum(sampled_sequences)
    total_nodes = 2 * total_sequences - 1
    time = zeros(Float64, total_nodes)
    left = zeros(Int, total_nodes)
    right = zeros(Int, total_nodes)

    active_nodes = Int[]
    time_steps = [-Inf; sequence_times]
    new_sequences = [0; sampled_sequences]
    node_idx = 1
    for t_idx in length(time_steps):-1:2
        current_time = time_steps[t_idx]
        # Add new sampled lineages to tree
        for _ in 1:new_sequences[t_idx]
            time[node_idx] = current_time
            # insert_random!(active_nodes, node_idx)
            insert!(active_nodes, rand(1:length(active_nodes)+1), node_idx)
            node_idx += 1
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
            # insert_random!(active_nodes, node_idx)
            insert!(active_nodes, rand(1:length(active_nodes)+1), node_idx)
            node_idx += 1
            n_active -= 1
            time_to_coalescence = rand(Exponential(2 * Ne / (n_active * (n_active - 1))))
        end
    end
    return time, left, right
end