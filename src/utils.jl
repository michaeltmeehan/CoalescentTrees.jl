"""
    pop_random!(v::Vector)

Remove and return a random element from vector `v`, swapping with the last element
for efficiency. Operates in constant time.

# Arguments
- `v::Vector`: input vector (must be non-empty)

# Returns
- `el`: a randomly chosen element from `v`

# Throws
- `ArgumentError` if the input vector is empty
"""
function pop_random!(v::Vector)
    isempty(v) && throw(ArgumentError("Cannot pop from an empty vector."))

    i = rand(1:length(v))
    v[i], v[end] = v[end], v[i]   # swap with last element
    return pop!(v)
end


"""
    reverse_cumsum(x::AbstractVector{T}) where T

Compute the reverse cumulative sum of vector `x`, such that the output at index `i`
is the sum of elements from `x[i]` to the end.

Equivalent to `reverse(cumsum(reverse(x)))`, but more efficient.

# Arguments
- `x::AbstractVector{T}`: input vector of numeric type `T`

# Returns
- `Vector{T}`: reverse cumulative sum of `x`, same length and type as `x`

# Example
```julia
reverse_cumsum([1, 2, 3])  # returns [6, 5, 3]
"""
function reverse_cumsum(x::AbstractVector{T}) where T
    n = length(x)
    out = similar(x, T)
    s = zero(T)
    for i in n:-1:1
        s += x[i]
        out[i] = s
    end
    return out
end