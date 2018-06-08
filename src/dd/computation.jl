include_dependency("core.jl")

"""
Find a longest path on DD.

# Input
- `dd::DecisionDiagram`: The DD.
- `func::Array{Function}`: An array of weight functions (each component represents the separable function for that variable)

# Output
- The longest path (composed of arc labels)
- The length of the longest path
"""
function longest_path(dd::DecisionDiagram, func::Array{Function})

    n = nvars(dd)
    @assert(n == length(func))       #makes sure the size of the fractional point matches the dimension of variables represented by dd

    node_num = nnodes(g)
    node_obj_vals = Array{Float64}(node_num)   #stores the maximum objective value computed from the source to each node
    node_obj_paths = Array{Int}(node_num)   #at each node, stores the index of the previous node of longest path from source to that node
    node_obj_labels = Array{Int}(node_num)  #stores the label value of the last arc of the longest path from source to each node

    #setting the root objective value equal to zero
    node_obj_vals[root(dd)] = 0

    #computing the longest path from source to a node, by searching through all nodes and their incoming arcs, layer by layer
    for i=2:n+1, node in dd.layers[i]
        temp_v = -Inf
        temp_p = 0
        temp_l = 0
        for (parent, label) in inneighbors(dd, node)
            val = node_obj_vals[parent] + func[i-1](label)    #computes the head node value based on the tail and arc label
            if val > temp_v
                temp_v = val
                temp_p = parent
                temp_l = label
            end
        end
        node_obj_vals[node] = temp_v
        node_obj_paths[node] = temp_p
        node_obj_labels[node] = temp_l
    end

    terminal_node = terminal(dd)
    lpv = node_obj_vals[terminal_node]   #the length of the longest path

    lp = Array{Int}(0)
    p = terminal_node
    #computing the arc labels on the longest path through a backtracking
    for i=n:-1:1
        # push!(lp, get_node_obj_label(dd, p))
        push!(lp, node_obj_labels[p])
        # p = get_node_obj_path(dd, p)
        p = node_obj_paths[p]
    end

    return reverse(lp), lpv
end



"""
Second method for function Longest_Path, where weight functions are linear for all variables.

# Input
- `dd::DecisionDiagram`: The DD.
- `cf::Array{T} where T<:Float64`: An array of coefficients of the variables in the linear objective function.

function longest_path(dd::DecisionDiagram, cf::Array{T} where T<:Float64)
    func = Array{Function}(length(cf))
    for (i, c) in enumerate(cf)                 #creates an array of separable functions for each variable
        func[i] = (x::Int64) -> c*x
    end
    return longest_path(dd, func)
end
"""


function longest_path(dd::DecisionDiagram, cf::Array{T} where T<:Float64)

    n = nvars(dd)
    @assert(n == length(cf))       #makes sure the size of the fractional point matches the dimension of variables represented by dd

    g = dd.graph
    node_num = nv(g)
    arc_num = ne(g)

    #getting the index of the source and the terminal node
    source_ind = dd.layers[1][1]            #gets the index value of the source node

    node_obj_vals = Array{Float64}(node_num)   #stores the maximum objective value computed from the source to each node
    node_obj_paths = Array{Int}(node_num)   #at each node, stores the index of the previous node of longest path from source to that node
    node_obj_labels = Array{Int}(node_num)  #stores the label value of the last arc of the longest path from source to each node

    #setting the source objective value equal to zero
    # set_node_obj_val!(dd, source_ind, 0)
    node_obj_vals[source_ind] = 0

    #computing the longest path from source to a node, by searching through all nodes and their incoming arcs, layer by layer
    for i=2:n+1, j in dd.layers[i]
        temp_v = -Inf
        temp_p = 0
        temp_l = 0
        for k in in_neighbors(g, j), l in get_arc_label(dd, k, j)
            val = node_obj_vals[k] + cf[i-1]*l    #computes the head node value based on the tail and arc label
            if val > temp_v
                temp_v = val
                temp_p = k
                temp_l = l
            end
        end
        # set_node_obj_val!(dd, j, temp_v)
        # set_node_obj_path!(dd, j, temp_p)
        # set_node_obj_label!(dd, j, temp_l)
        node_obj_vals[j] = temp_v
        node_obj_paths[j] = temp_p
        node_obj_labels[j] = temp_l
    end

    terminal_ind = dd.layers[end][1]
    # lpv = get_node_obj_val(dd, terminal_ind)   #the length of the longest path
    lpv = node_obj_vals[terminal_ind]   #the length of the longest path

    lp = Array{Int}(0)
    p = terminal_ind
    #computing the arc labels on the longest path through a backtracking
    for i=n:-1:1
        # push!(lp, get_node_obj_label(dd, p))
        push!(lp, node_obj_labels[p])
        # p = get_node_obj_path(dd, p)
        p = node_obj_paths[p]
    end

    return reverse(lp), lpv

end
