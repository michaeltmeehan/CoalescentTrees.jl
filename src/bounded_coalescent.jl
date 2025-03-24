function coalescence_probability(n_big::Int, n_small::Int, Δt::Float64, Ne::Float64)
    @assert Ne > 0 "Effective population size (Ne) must be positive."
    @assert Δt ≥ 0 "Time (Δt) must be non-negative."
    @assert n_big ≥ n_small "n_big must be greater than or equal to n_small."
    @assert n_big ≥ 1 "n_big must be at least 1."
    @assert n_small ≥ 1 "n_small must be at least 1."

    if n_small == 1
        sum_prob = 0.0
        for k in 2:n_big
            partial_prob = exp(-k * (k - 1.0) * Δt / (2.0 * Ne))
            for l in 2:n_big
                l == k && continue
                partial_prob *= l * (l - 1.0) / (l * (l - 1.0) - k * (k - 1.0))
            end
            sum_prob += partial_prob
        end
        out = 1.0 - sum_prob
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
    if !(0.0 ≤ out ≤ 1.0)
        println("Warning: Coalescence probability is not in [0, 1].")
        println("Trying stable version...")
        out = coalescence_probability_stable(n_big, n_small, Δt, Ne)
    end
    @assert 0. ≤ out ≤ 1.0 "Coalescence probability must be in [0, 1]."
    return out
end


function stable_partial_probability(k, j, n_big, Δt, Ne)
    log_p = -k * (k - 1.0) * Δt / (2.0 * Ne)
    for l in j:n_big
        l == k && continue
        log_p -= log(abs(1. - k * (k - 1.0) / (l * (l - 1.0))))
    end
    return (-1.)^(k - j) * exp(log_p)
end


function coalescence_probability_stable(n_big::Int, n_small::Int, Δt::Float64, Ne::Float64)
    if n_small == 1
        sum_prob = mapreduce(k -> stable_partial_probability(k, 2, n_big, Δt, Ne), +, 2:n_big)
        return 1.0 - sum_prob
    else
        sum_prob = mapreduce(k -> k * (k - 1.0) / 2.0 * stable_partial_probability(k, n_small, n_big, Δt, Ne), +, n_small:n_big)
        return (2.0 / (n_small * (n_small - 1.0))) * sum_prob
    end
end