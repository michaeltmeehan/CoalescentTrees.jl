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

    if n_small == 1
        prob_total = 0.0
        for k in 2:n_big
            prob_increment = exp(-k * (k - 1.0) * Δt / (2.0 * Ne))
            # prob_increment = 1.0
            for l in 2:n_big
                l == k && continue
                prob_increment *= l * (l - 1.0) / (l * (l - 1.0) - k * (k - 1.0))
                # prob_increment *= l * (l - 1.0) / (l * (l - 1.0) - k * (k - 1.0))
            end
            # prob_increment *= 1.0 - exp(-k * (k - 1.0) * Δt / (2.0 * Ne))
            prob_total += prob_increment
        end
        out = 1.0 - prob_total
        # out = prob_total
    else
        sum_prob = 0.0
        for k in n_small:n_big
            partial_prob = k * (k - 1.0) / 2.0 * exp(-k * (k - 1.0) * Δt / (2.0 * Ne))
            for l in n_small:n_big
                l == k && continue
                partial_prob *= l * (l - 1.0) / (l * (l - 1.0) - k * (k - 1.0))
            end
            sum_prob += partial_prob
        end
        out = 2.0 / (n_small * (n_small - 1.0)) * sum_prob
    end
    abs(out) ≤ 1e-12 && return 0.0
    # if !(0.0 ≤ out ≤ 1.0)
    #     println("Warning: Coalescence probability is not in [0, 1].")
    #     println("Trying stable version...")
    #     out = coalescence_probability_stable(n_big, n_small, Δt, Ne)
    # end
    # @assert 0. ≤ out ≤ 1.0 "Coalescence probability must be in [0, 1]."
    return out
end

# TODO: Need to implement the stable version of the coalescence_probability function
# Values become negative when Δt / Ne is small
function stable_partial_probability(k, j, n_big, Δt, Ne)
    log_p = -k * (k - 1.0) * Δt / (2.0 * Ne)
    for l in j:n_big
        l == k && continue
        log_p -= log(abs(1. - k * (k - 1.0) / (l * (l - 1.0))))
    end
    return (-1.)^(k - j) * exp(log_p)
end


function coalescent_probability_stable(n_big::Int, n_small::Int, Δt::Float64, Ne::Float64)
    # (n_big == 1 && n_small = 1) && return 1.0
    if n_small == 1
        sum_prob = mapreduce(k -> stable_partial_probability(k, 2, n_big, Δt, Ne), +, 2:n_big)
        return 1.0 - sum_prob
    else
        sum_prob = mapreduce(k -> k * (k - 1.0) / 2.0 * stable_partial_probability(k, n_small, n_big, Δt, Ne), +, n_small:n_big)
        return (2.0 / (n_small * (n_small - 1.0))) * sum_prob
    end
end


"""
    calc_forward_probs(sampled_sequences, sequence_times, bound_time, Ne)

Run the forward step of the forward-backward algorithm for the bounded coalescent model.

# Arguments
- `sampled_sequences::Vector{Int}`: number of sequences sampled at each time point
- `sequence_times::Vector{Float64}`: sampling times (descending order)
- `bound_time::Float64`: lower time bound (oldest point in time)
- `Ne::Float64`: effective population size

# Returns
- A matrix `P[k, t]` giving the probability of having `k` lineages at time `t`
"""
function calc_forward_probs(sampled_sequences::Vector{Int}, sequence_times::Vector{Float64}, bound_time::Float64, Ne::Float64)
    n = length(sampled_sequences)
    total_sequences = sum(sampled_sequences)
    @assert n == length(sequence_times) "Number of sampled sequences and sequence times must match."
    @assert total_sequences > 1 "At least two sequences are required."
    @assert Ne > 0 "Effective population size (Ne) must be positive."

    # Concatenate bound time and sequence times
    times = [bound_time; sequence_times]
    new_lineages = [0; sampled_sequences]
    n_timesteps = n + 1

    # Initialize forward pass
    forward_probabilities = zeros(Float64, total_sequences, n_timesteps)
    forward_probabilities[new_lineages[end], end] = 1.0
    max_lineages = reverse(cumsum(sampled_sequences))

    # Forward pass
    for t in n:-1:1
        Δt = times[t+1] - times[t]
        for j in 1:max_lineages[t]
            for i in j:max_lineages[t]
                # Number of lineages at time t: j lineages formed from i at t+1 plus new lineages
                forward_probabilities[j+new_lineages[t], t] += coalescent_probability(i, j, Δt, Ne) * forward_probabilities[i, t+1]
            end
        end
    end
    @assert (all(forward_probabilities .≥ 0.)  && all(forward_probabilities .≤ 1.) && all(sum(forward_probabilities, dims=1) .≈ 1.0)) "Failed on forward pass"
    return forward_probabilities
end


function calc_backward_probs(forward_probs::Matrix{Float64}, sampled_sequences::Vector{Int}, sequence_times::Vector{Float64}, bound_time::Float64, Ne::Float64)
    n = length(sampled_sequences)
    total_sequences = sum(sampled_sequences)
    @assert n == length(sequence_times) "Number of sampled sequences and sequence times must match."
    @assert total_sequences > 1 "At least two sequences are required."
    @assert Ne > 0 "Effective population size (Ne) must be positive."

    # Concatenate bound time and sequence times
    times = [bound_time; sequence_times]
    new_lineages = [0; sampled_sequences]
    n_timesteps = n + 1

    max_lineages = [reverse(cumsum(sampled_sequences)); 0]

    # Initialize backward pass for sampling
    drawn_lineages = zeros(Int, n_timesteps)
    backward_probabilities = zeros(Float64, total_sequences, n_timesteps)
    backward_probabilities[1, 1] = 1.0
    drawn_lineages[1] = 1

    # Backward pass
    for t in 2:n_timesteps
        Δt = times[t] - times[t-1]
        for j in (drawn_lineages[t-1] - 2*new_lineages[t-1]):max_lineages[t]
            backward_probabilities[j+new_lineages[t], t] = coalescent_probability(j + new_lineages[t], drawn_lineages[t-1] - new_lineages[t-1],  Δt, Ne) * 
                                                                forward_probs[j + new_lineages[t], t] / forward_probs[drawn_lineages[t-1], t-1]
        end
        drawn_lineages[t] = sample(1:total_sequences, Weights(backward_probabilities[:, t]))
    end
    @assert (all(backward_probabilities .≥ 0.)  && all(backward_probabilities .≤ 1.) && all(sum(backward_probabilities, dims=1) .≈ 1.0)) "Failed on backward pass"
    return drawn_lineages, backward_probabilities
end


function partition_intervals(extant_lineages::Vector{Int}, sampled_sequences::Vector{Int}, sequence_times::Vector{Float64}, bound_time::Float64, Ne::Float64)
    times = [bound_time; sequence_times]
    samples = vcat([0], sampled_sequences)
    lineages = copy(extant_lineages)

    i = 1
    while i < length(times)
        n_big = lineages[i+1]   # Number of lineages at the end of the interval (times[i])
        n_small = lineages[i] - samples[i]  # Number of lineages at the beginning of the interval (times[i+1])
        Δn = n_big - n_small
        if Δn > 1
            # Compute midpoint of the interval
            t_mid = (times[i] + times[i+1]) / 2

            # Sample number of coalescent events in each half
            Δt = times[i+1] - times[i]

            # Calculate probability for lineage counts at midpoint
            prob_a_mid = coalescent_probability.(n_small .+ (0:Δn), n_small, Δt / 2, Ne) .* 
                            coalescent_probability.(n_big, n_small .+ (0:Δn), Δt / 2, Ne) / 
                                coalescent_probability(n_big, n_small, Δt, Ne)

            # Draw number of coalescent events in in the interval [t_i, t_mid]
            n_coalescent = sample(0:Δn, Weights(prob_a_mid))

            # Insert midpoint and update counts
            insert!(times, i+1, t_mid)
            insert!(lineages, i+1, n_small + n_coalescent)
            insert!(samples, i+1, 0)

            # Don't increment `i` since the new subintervals may still need splitting
        else
            i += 1  # Move to the next interval only if it does not need further splitting
        end
    end

    return times, lineages, samples
end


function sample_coalescence_time(t_lower::Float64, t_upper::Float64, n_small::Int, Ne::Float64)
    r = n_small / Ne
    return 1/r * log(exp(r * t_upper) - rand() * (exp(r * t_upper) - exp(r * t_lower)))
end


function sample_coalescent_times(interval_times::Vector{Float64}, lineages::Vector{Int}, samples::Vector{Int}, Ne::Float64)
    n_intervals = length(interval_times) - 1
    coalescence_times = Vector{Float64}()
    for i in 1:n_intervals
        n_big = lineages[i+1]   # Number of lineages at the end of the interval (times[i])
        n_small = lineages[i] - samples[i]  # Number of lineages at the beginning of the interval (times[i+1])
        if n_big > n_small
            t = sample_coalescence_time(interval_times[i], interval_times[i+1], n_small, Ne)
            push!(coalescence_times, t)
        end
    end
    return coalescence_times
end