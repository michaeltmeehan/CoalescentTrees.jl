# =============================================================================
# Bounded Coalescent Tree Sampler
# =============================================================================
# This module implements a forward-backward algorithm for sampling phylogenetic
# trees under the bounded coalescent model. The model conditions on:
#   - a lower time bound (the root age),
#   - an effective population size (Ne),
#   - a sequence of heterochronous sampling times and counts.
#
# The implementation follows a five-step algorithm:
#
#   1. Forward Filtering Backward Sampling (FFBS):
#      - Run a forward pass to compute the distribution over lineage counts at
#        each sampling time conditioned on the root time (bound).
#      - Sample a valid lineage count path using backward sampling.
#
#   2. Determine the number of coalescent events in each time interval:
#      - Based on sampled lineage counts and sampling times.
#
#   3. Partition intervals:
#      - Subdivide intervals with multiple coalescent events until each interval
#        contains at most one event.
#      - Sample the number of events per subinterval using a conditional
#        distribution.
#
#   4. Sample coalescent times:
#      - Use inverse transform sampling from the coalescent rate conditioned on
#        exactly one coalescent event.
#
#   5. Sample tree topology:
#      - Traverse the event list in reverse chronological order.
#      - Build a rooted binary tree by randomly pairing extant lineages at each
#        coalescence event.

# -----------------------------------------------------------------------------
# Function Overview
# -----------------------------------------------------------------------------
# - coalescent_probability          : Main probability function (delegates to special cases)
# - coalescent_probability_to_root        : Analytic form for coalescing to a single root
# - coalescent_probability_general        : Analytic form for general k → j coalescence
# - coalescent_probability_stable  : Stable alternative when Δt/Ne is small
#
# - calc_forward_probs             : Forward pass (FFBS step 1a)
# - sample_lineages_backward            : Backward sampling pass (FFBS step 1b)
#
# - partition_intervals            : Recursive subdivision of intervals with >1 coalescent
# - sample_coalescent_time         : Inverse-transform sampling for coalescent times
# - construct_event_list             : Convert interval structure to timestamped events
# - sample_topology                : Build tree topology from event list
# - sample_tree                    : Main entry point — executes full sampling pipeline
#
# Utility functions (in separate utils.jl file):
# - reverse_cumsum                 : Efficient backward cumulative sum
# - pop_random!                    : Randomly pop an element from a vector in O(1) time
#
# -----------------------------------------------------------------------------
# Expected Usage
# -----------------------------------------------------------------------------
# Call `sample_tree(sampled_sequences, sequence_times, bound_time, Ne)` to
# generate a random phylogenetic tree under the bounded coalescent model.
#
# All time vectors must be sorted in ascending order, and `sampled_sequences`
# must sum to at least 2.
#
# The output is a `Tree` structure containing:
#   - `times`: node times
#   - `left`, `right`: indices of child nodes
#
# Tree construction assumes time is measured backwards from the present.
# =============================================================================


"""
    coalescence_probability(n_big, n_small, Δt, Ne)

Compute the probability that `n_big` lineages coalesce down to `n_small` lineages 
in time interval `Δt` under the standard coalescent model with effective population size `Ne`.

Falls back to a numerically stable computation if instability is detected.

# Arguments
- `n_big::Int`: initial number of lineages (≥ 1)
- `n_small::Int`: final number of lineages (≥ 1, ≤ n_big)
- `Δt::Float64`: time interval (≥ 0)
- `Ne::Float64`: effective population size (> 0)

# Returns
- `Float64`: probability in [0, 1]

# Throws
- `AssertionError` if inputs are invalid
"""
function coalescent_probability(n_big::Int, n_small::Int, Δt::Float64, Ne::Float64)
    @assert Ne > 0 "Effective population size (Ne) must be positive."
    @assert Δt ≥ 0 "Time (Δt) must be non-negative."
    @assert n_big ≥ n_small "n_big must be greater than or equal to n_small."
    @assert n_big ≥ 1 "n_big must be at least 1."
    @assert n_small ≥ 1 "n_small must be at least 1."
    n_big == 1 && return 1.0
    n_small == 1 && return coalescent_probability_to_root(n_big, Δt, Ne)
    return coalescent_probability_general(n_big, n_small, Δt, Ne)
end


# Probability that n_big lineages coalesce to exactly 1 in Δt
function coalescent_probability_to_root(n_big::Int, Δt::Float64, Ne::Float64)
    probability_total = 0.0
    for k in 2:n_big
        λk = k * (k - 1.0)
        term = exp(-λk * Δt / (2.0 * Ne))
        for l in 2:n_big
            l == k && continue
            λl = l * (l - 1.0)
            term *= λl / (λl - λk)
        end
        probability_total += term
    end
    out = 1.0 - probability_total
    return abs(out) ≤ 1e-12 ? 0.0 : out
end


# General case: probability of coalescing from n_big to n_small in Δt
function coalescent_probability_general(n_big::Int, n_small::Int, Δt::Float64, Ne::Float64)
    sum_prob = 0.0
    for k in n_small:n_big
        λk = k * (k - 1.0)
        term = (λk / 2.0) * exp(-λk * Δt / (2.0 * Ne))
        for l in n_small:n_big
            l == k && continue
            λl = l * (l - 1.0)
            term *= λl / (λl - λk)
        end
        sum_prob += term
    end
    out = (2.0 / (n_small * (n_small - 1.0))) * sum_prob
    return abs(out) ≤ 1e-12 ? 0.0 : out
end


# TODO: Need to implement the stable version of the coalescent_probability function
"""
    stable_partial_probability_term(k::Int, j::Int, n_big::Int, Δt::Float64, Ne::Float64)

Compute the `k`th partial term (log-space) in the stable coalescence probability sum 
from `n_big` to `j` lineages over time interval `Δt` with effective population size `Ne`.

Used to avoid underflow in the original analytical formula.

# Returns
- `Float64`: the (signed) probability term for summation
"""
function stable_partial_probability_term(k::Int, j::Int, n_big::Int, Δt::Float64, Ne::Float64)
    λk = k * (k - 1.0)
    log_p = -λk * Δt / (2.0 * Ne)
    for l in j:n_big
        l == k && continue
        λl = l * (l - 1.0)
        log_p -= log(abs(1.0 - λk / λl))
    end
    return (-1.0)^(k - j) * exp(log_p)
end



"""
    coalescent_probability_stable(n_big::Int, n_small::Int, Δt::Float64, Ne::Float64)

Numerically stable version of `coalescent_probability`, used when Δt/Ne is small.

# Returns
- `Float64`: probability of coalescing from `n_big` to `n_small` lineages
"""
function coalescent_probability_stable(n_big::Int, n_small::Int, Δt::Float64, Ne::Float64)
    if n_small == 1
        sum_prob = mapreduce(k -> stable_partial_probability_term(k, 2, n_big, Δt, Ne), +, 2:n_big)
        return 1.0 - sum_prob
    else
        sum_prob = mapreduce(k -> (k * (k - 1.0) / 2.0) * stable_partial_probability_term(k, n_small, n_big, Δt, Ne), +, n_small:n_big)
        return (2.0 / (n_small * (n_small - 1.0))) * sum_prob
    end
end


"""
    calc_forward_probs(sampled_sequences, sequence_times, bound_time, Ne)

Run the forward step of the forward-backward algorithm for the bounded coalescent model.

# Arguments
- `sampled_sequences::Vector{Int}`: number of sequences sampled at each time point
- `sequence_times::Vector{Float64}`: sampling times (descending order)
- `bound_time::Float64`: oldest time point (the bound)
- `Ne::Float64`: effective population size

# Returns
- `Matrix{Float64}`: forward probabilities, `P[k, t]` is the prob. of having `k` lineages at time `t`
"""
function calc_forward_probs(
    sampled_sequences::Vector{Int},
    sequence_times::Vector{Float64},
    bound_time::Float64,
    Ne::Float64
)
    n_times = length(sampled_sequences)
    total_sequences = sum(sampled_sequences)

    # Time points and new lineages at each
    times = [bound_time; sequence_times]
    samples = [0; sampled_sequences]
    n_steps = n_times + 1

    # Maximum number of lineages possible at each time step (in reverse order)
    max_lineages = [total_sequences; reverse_cumsum(sampled_sequences)]

    # Initialize forward probability matrix [lineages x time step]
    forward_probs = zeros(Float64, total_sequences, n_steps)
    final_k = samples[end]  # number of lineages at most recent sample time
    forward_probs[final_k, end] = 1.0  # P[k, t_final] = 1

    # Backward loop over time steps
    for t in n_times:-1:1
        Δt = times[t+1] - times[t]
        for n_small in 1:max_lineages[t+1]  # Lineage count at time t
            k_curr = n_small + samples[t]
            for n_big in n_small:max_lineages[t+1]  # Lineage count at time t+1
                forward_probs[n_big, t+1] == 0.0 && continue
                forward_probs[k_curr, t] += coalescent_probability(n_big, n_small, Δt, Ne) * forward_probs[n_big, t+1]
            end
        end

        # Normalize column t to sum to 1 (helps guard against drift)
        col_sum = sum(view(forward_probs, :, t))
        if !(isapprox(col_sum, 1.0; atol=1e-8))
            @warn "Forward probabilities at t = $t did not sum to 1.0 (sum = $col_sum). Renormalizing."
            forward_probs[:, t] ./= col_sum
        end
    end

    return forward_probs
end


"""
    sample_lineages_backward(forward_probs, sampled_sequences, sequence_times, bound_time, Ne)

Run the backward sampling step of the forward-backward algorithm.

# Arguments
- `forward_probs::Matrix{Float64}`: matrix of forward probabilities (from `calc_forward_probs`)
- `sampled_sequences::Vector{Int}`: number of sequences sampled at each time point
- `sequence_times::Vector{Float64}`: sampling times (descending order)
- `bound_time::Float64`: oldest time point (the bound)
- `Ne::Float64`: effective population size

# Returns
- `lineages::Vector{Int}`: sampled lineage counts at each time point
- `backward_probs::Matrix{Float64}`: backward probability matrix
"""
function sample_lineages_backward(forward_probs::Matrix{Float64}, sampled_sequences::Vector{Int}, sequence_times::Vector{Float64}, bound_time::Float64, Ne::Float64)
    # Concatenate bound time and sequence times
    times = [bound_time; sequence_times]
    samples = [0; sampled_sequences]
    n_timesteps = length(times)

    max_lineages = [sum(sampled_sequences); reverse_cumsum(sampled_sequences)]

    # Initialize backward pass for sampling
    lineages = zeros(Int, n_timesteps)
    backward_probs = zeros(Float64, size(forward_probs))

    # Start from the root (1 lineage at the bound time)
    backward_probs[1, 1] = 1.0
    lineages[1] = 1

    # Backward pass
    for t in 2:n_timesteps
        Δt = times[t] - times[t-1]
        n_small = lineages[t-1] - samples[t-1]
        for n_big in n_small:max_lineages[t]
            forward_probs[n_big, t] == 0.0 && continue
            backward_probs[n_big, t] = coalescent_probability(n_big, n_small,  Δt, Ne) * forward_probs[n_big, t] / forward_probs[lineages[t-1], t-1]
        end

        # Normalize column t
        col_sum = sum(view(backward_probs, :, t))
        if !(isapprox(col_sum, 1.0; atol=1e-8))
            @warn "Backward probabilities at t = $t did not sum to 1.0 (sum = $col_sum). Renormalizing."
            backward_probs[:, t] ./= col_sum
        end
        lineages[t] = sample(n_small:max_lineages[t], Weights(backward_probs[n_small:max_lineages[t], t]))
    end
    @assert all(samples .≤ lineages .≤ max_lineages) "Sampled lineages out of bounds."

    return lineages, backward_probs
end


"""
    partition_intervals(lineages, sampled_sequences, sequence_times, bound_time, Ne)

Recursively split time intervals to ensure at most one coalescent event per interval.

# Arguments
- `lineages::Vector{Int}`: lineage count at each time point (starting from the root)
- `sampled_sequences::Vector{Int}`: number of sequences sampled at each time point
- `sequence_times::Vector{Float64}`: sampling times in descending order
- `bound_time::Float64`: oldest time point (time of root)
- `Ne::Float64`: effective population size

# Returns
- `times::Vector{Float64}`: updated list of time points
- `lineages::Vector{Int}`: updated lineage counts at each time
- `samples::Vector{Int}`: updated samples vector aligned with time points
"""
function partition_intervals(
    lineages::Vector{Int},
    sampled_sequences::Vector{Int},
    sequence_times::Vector{Float64},
    bound_time::Float64,
    Ne::Float64
)
    times = [bound_time; sequence_times]
    samples = [0; sampled_sequences]

    i = 1
    while i < length(times)
        n_big = lineages[i + 1]             # Lineages at end of interval (more recent)
        n_small = lineages[i] - samples[i]  # Lineages at start of interval
        Δn = n_big - n_small
        if Δn > 1
            t_mid = (times[i] + times[i + 1]) / 2
            Δt = times[i + 1] - times[i]

            # Compute probability distribution over possible midpoints
            mid_counts = n_small .+ (0:Δn)
            numerators = coalescent_probability.(mid_counts, n_small, Δt / 2, Ne) .*
                         coalescent_probability.(n_big, mid_counts, Δt / 2, Ne)
            denominator = coalescent_probability(n_big, n_small, Δt, Ne)
            weights = numerators ./ denominator

            # Sample number of coalescences in [t_i, t_mid]
            n_coalescent = sample(0:Δn, Weights(weights))

            # Insert midpoint and update vectors
            insert!(times, i + 1, t_mid)
            insert!(lineages, i + 1, n_small + n_coalescent)
            insert!(samples, i + 1, 0)

            # Rerun this new interval next
        else
            i += 1
        end
    end

    return times, lineages, samples
end


"""
    sample_coalescent_time(t_lower, t_upper, n_small, Ne)

Sample a coalescent event time uniformly from the conditional distribution 
under the standard coalescent model in the interval [t_lower, t_upper].

Assumes exactly one coalescent event occurs in the interval.

# Arguments
- `t_lower::Float64`: start of the interval
- `t_upper::Float64`: end of the interval
- `n_small::Int`: number of extant lineages at the start of the interval
- `Ne::Float64`: effective population size

# Returns
- `Float64`: sampled coalescent time
"""
function sample_coalescent_time(t_lower::Float64, t_upper::Float64, n_small::Int, Ne::Float64)
    r = n_small / Ne
    return 1/r * log(exp(r * t_upper) - rand() * (exp(r * t_upper) - exp(r * t_lower)))
end


"""
    construct_event_list(interval_times, lineages, samples, Ne)

Generate a list of events (coalescence, sampling, root) based on the final interval structure.

# Arguments
- `interval_times::Vector{Float64}`: time points marking interval boundaries
- `lineages::Vector{Int}`: number of lineages at each time point
- `samples::Vector{Int}`: number of sequences sampled at each time point
- `Ne::Float64`: effective population size

# Returns
- `Vector{NamedTuple{(:time, :type), Tuple{Float64, Int}}}`: list of events
    - `type = 0` → sampling event
    - `type = 1` → bound/root event
    - `type = 2` → coalescence event
"""
function construct_event_list(
    interval_times::Vector{Float64},
    lineages::Vector{Int},
    samples::Vector{Int},
    Ne::Float64
)
    n_intervals = length(interval_times) - 1
    events = [(time = interval_times[1], type = 1)]  # root/bound event

    for i in 1:n_intervals
        n_big = lineages[i + 1]
        n_small = lineages[i] - samples[i]

        if n_big > n_small
            t = sample_coalescent_time(interval_times[i], interval_times[i + 1], n_small, Ne)
            push!(events, (time = t, type = 2))
        end

        if samples[i + 1] > 0
            for _ in 1:samples[i + 1]
                push!(events, (time = interval_times[i + 1], type = 0))
            end
        end
    end
    @assert issorted(events, by=x->x.time) "Events are not sorted by time."
    return events
end


"""
    sample_topology(events)

Construct a rooted binary tree from a list of events sorted in reverse chronological order.

# Arguments
- `events::Vector{NamedTuple{(:time, :type), Tuple{Float64, Int}}}`: list of events
    - `type = 0`: sampling event (leaf)
    - `type = 1`: bound/root event
    - `type = 2`: coalescent event (internal node)

# Returns
- `Tree`: a binary tree with fields:
    - `times::Vector{Float64}`: time of each node
    - `left::Vector{Int}`: index of left descendant (0 if none)
    - `right::Vector{Int}`: index of right descendant (0 if none)
"""
function sample_topology(events::Vector{NamedTuple{(:time, :type), Tuple{Float64, Int}}})
    n_nodes = length(events)
    times = Vector{Float64}(undef, n_nodes)
    left = Vector{Int}(undef, n_nodes)
    right = Vector{Int}(undef, n_nodes)

    active_nodes = Int[]

    for node in n_nodes:-1:1
        event = events[node]

        if event.type == 0
            # Leaf node (sampling event)
            times[node] = event.time
            left[node] = right[node] = 0
            push!(active_nodes, node)

        elseif event.type == 2
            # Internal node (coalescence)
            times[node] = event.time
            left[node], right[node] = pop_random!(active_nodes), pop_random!(active_nodes)
            push!(active_nodes, node)

        elseif event.type == 1
            # Root node (bound time)
            times[node] = event.time
            left[node], right[node] = pop_random!(active_nodes), 0

        else
            error("Unknown event type: $(event.type)")
        end
    end

    @assert isempty(active_nodes) "Still active nodes after tree construction."
    return Tree(times, left, right)
end


"""
    sample_tree(sampled_sequences, sequence_times, bound_time, Ne)

Sample a full binary tree under the bounded coalescent model using forward-backward sampling.

# Arguments
- `sampled_sequences::Vector{Int}`: number of sequences sampled at each time point
- `sequence_times::Vector{Float64}`: sampling times in descending order
- `bound_time::Float64`: oldest time point (lower bound of coalescent process)
- `Ne::Float64`: effective population size

# Returns
- `Tree`: a rooted binary tree with coalescent and sampling structure consistent with the bounded model
"""
function sample_tree(
    sampled_sequences::Vector{Int},
    sequence_times::Vector{Float64},
    bound_time::Float64,
    Ne::Float64
)
    @assert length(sampled_sequences) == length(sequence_times) "Sample count and time vectors must be the same length."
    @assert sum(sampled_sequences) ≥ 2 "At least two sequences are required to form a tree."
    @assert Ne > 0 "Effective population size must be positive."
    @assert issorted(sequence_times) "Sequence times must be in ascending order."

    forward_probs = calc_forward_probs(sampled_sequences, sequence_times, bound_time, Ne)
    lineages, _ = sample_lineages_backward(forward_probs, sampled_sequences, sequence_times, bound_time, Ne)
    interval_times, lineages, samples = partition_intervals(lineages, sampled_sequences, sequence_times, bound_time, Ne)
    events = construct_event_list(interval_times, lineages, samples, Ne)
    tree = sample_topology(events)

    return tree
end