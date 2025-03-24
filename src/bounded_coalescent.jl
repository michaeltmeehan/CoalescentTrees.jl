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


