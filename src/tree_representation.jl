import Base: show


mutable struct Tree
    time::Vector{Float64}
    left::Vector{Int}
    right::Vector{Int}
    data::Dict{Int, Any}
end


function show(io::IO, tree::Tree)
    n_nodes = length(tree.time)
    n_tips = sum(tree.left .== 0 .&& tree.right .== 0)
    has_data = !isempty(tree.data)

    println(io, "Tree with $n_nodes nodes ($n_tips tips)")
    n_nodes > 0 && println(io, "Root height: $(tree.time[1])")
    println(io, has_data ? "Metadata attached to $(length(tree.data)) nodes" : "No node metadata attached")
end


Tree(time::Vector{Float64}, left::Vector{Int}, right::Vector{Int}) = Tree(time, left, right, Dict{Int, Any}())

Tree() = Tree(Float64[], Int[], Int[], Dict{Int, Any}())

Tree(triple::Tuple{Vector{Float64}, Vector{Int}, Vector{Int}}) = Tree(triple[1], triple[2], triple[3])

Base.getindex(tree::Tree, i::Int) = (
    time = tree.time[i],
    left = tree.left[i],
    right = tree.right[i],
    data = get(tree.data, i, nothing)
)


function Base.getindex(tree::Tree, inds::AbstractVector{<:Integer})
    !issorted(inds) && sort!(inds)
    index_map = Dict(old => new for (new, old) in enumerate(inds))

    new_time  = tree.time[inds]
    new_left  = [get(index_map, tree.left[i], 0) for i in inds]
    new_right = [get(index_map, tree.right[i], 0) for i in inds]

    new_data = Dict(index_map[i] => tree.data[i] for i in inds if haskey(tree.data, i))

    return Tree(new_time, new_left, new_right, new_data)
end

Base.isequal(tree1::Tree, tree2::Tree) = tree1.time == tree2.time && tree1.left == tree2.left && tree1.right == tree2.right

Base.:(==)(tree1::Tree, tree2::Tree) = isequal(tree1, tree2)

Base.hash(tree::Tree, h::UInt) = hash((tree.time, tree.left, tree.right), h)

Base.length(tree::Tree) = length(tree.time)

Base.eachindex(tree::Tree) = eachindex(tree.time)

eachnode(tree::Tree; start::Int=1, stop::Int=length(tree)) = start:stop

in_tree(tree::Tree, node::Int) = 1 ≤ node ≤ length(tree)

is_leaf(tree::Tree, node::Int) = tree.left[node] == 0 && tree.right[node] == 0

function get_leaves(tree::Tree)
    leaves = Int[]
    for node in eachnode(tree)
        is_leaf(tree, node) && push!(leaves, node)
    end
    return leaves
end

function is_internal(tree::Tree, node::Int)
    @assert in_tree(tree, node) "Node $node is not in the tree"
    return !is_leaf(tree, node)
end


function is_binary(tree::Tree, node::Int)
    @assert in_tree(tree, node) "Node $node is not in the tree"
    return tree.left[node] != 0 && tree.right[node] != 0
end


function get_children(tree::Tree, node::Int)
    @assert in_tree(tree, node) "Node $node is not in the tree"
    return [child for child in [tree.left[node], tree.right[node]] if child != 0]
end


function get_parent(tree::Tree, node::Int)
    @assert in_tree(tree, node) "Node $node is not in the tree"
    for parent in eachnode(tree, stop=node-1)
        if tree.left[parent] == node || tree.right[parent] == node
            return parent
        end
    end
    return nothing
end


function get_ancestors(tree::Tree, node::Int)
    @assert in_tree(tree, node) "Node $node is not in the tree"
    ancestors = Int[]
    current_node = node
    for nd in reverse(eachnode(tree, stop=node-1))
        (tree.left[nd] == current_node || tree.right[nd] == current_node) && begin
            push!(ancestors, nd)
            current_node = nd
        end
    end
    return ancestors
end


function is_ancestor(tree::Tree, node::Int, ancestor::Int)
    @assert in_tree(tree, node) "Node $node is not in the tree"
    @assert in_tree(tree, ancestor) "Node $ancestor is not in the tree"
    return ancestor in get_ancestors(tree, node)
end


function get_descendants(tree::Tree, node::Int)
    @assert in_tree(tree, node) "Node $node is not in the tree"
    unvisited = [node]
    descendants = Int[]
    while !isempty(unvisited)
        current_node = pop!(unvisited)
        children = get_children(tree, current_node)
        append!(descendants, children)
        append!(unvisited, children)
    end
    return descendants
end


function get_subtree(tree::Tree, node::Int)
    @assert in_tree(tree, node) "Node $node is not in the tree"
    subtree_nodes = get_descendants(tree, node)
    return tree[[node; subtree_nodes]]
end


@recipe function f(tree::Tree; d=1.0)
    # Suppress legend and grid
    legend := false
    grid := false
    yticks := false

    # Initialize x and y coordinates
    x = tree.time
    y = zeros(Float64, length(x))
    y_counter = d

    # Stack for depth-first traversal
    stack = [(1, false)]

    # Traverse tree
    while !isempty(stack)
        node, children_visited = stack[end]
        left, right = tree.left[node], tree.right[node]
        if children_visited
            pop!(stack)
            if left == 0 && right == 0
                y[node] = y_counter
                y_counter += d
                @series begin
                    seriestype := :scatter
                    color := :blue
                    [x[node]], [y[node]]
                end
            elseif left == 0
                y[node] = y[right]
                @series begin
                    seriestype := :path
                    color := :black
                    [x[node], x[right]], [y[node], y[right]]
                end
            elseif right == 0
                y[node] = y[left]
                @series begin
                    seriestype := :path
                    color := :black
                    [x[node], x[left]], [y[node], y[left]]
                end
            else
                y[node] = (y[left] + y[right]) / 2
                @series begin
                    seriestype := :path
                    color := :black
                    [x[node], x[node]], [y[right], y[left]]
                end
                @series begin
                    seriestype := :path
                    color := :black
                    [x[node], x[left]], [y[left], y[left]]
                end
                @series begin
                    seriestype := :path
                    color := :black
                    [x[node], x[right]], [y[right], y[right]]
                end
            end
        else
            stack[end] = (node, true)
            left != 0 && push!(stack, (left, false))
            right != 0 && push!(stack, (right, false))
        end
    end
end