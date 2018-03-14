include("dd.jl")

"""
Function: Finds a longest path on DD.
Input: The DD, and an array of functions (each component represents the separable function for that variable)
Output: The longest path (composed of arc labels), and its length
"""
function longest_path!(dd::DecisionDiagram, func::Array{Function})

    n = size(dd.layers, 1) - 1     #the number of arc layers of dd
    @assert(n == length(func))       #makes sure the size of the fractional point matches the dimension of variables represented by dd

    g = dd.graph
    node_num = nv(g)
    arc_num = ne(g)

    #getting the index of the source and the terminal node
    source_ind = dd.layers[1][1]            #gets the index value of the source node

    node_obj_vals = Array{Real}(node_num)
    node_obj_paths = Array{Int}(node_num)
    node_obj_labels = Array{Int}(node_num)

    #setting the source objective value at zero
    # set_node_obj_val!(dd, source_ind, 0)
    node_obj_vals[source_ind] = 0

    #computing the longest path from source to a node, by searching through all nodes and their incoming arcs layer by layer
    for i=2:n+1, j in dd.layers[i]
        temp_v = -Inf
        temp_p = 0
        temp_l = 0
        for k in in_neighbors(g, j), l in get_arc_label(dd, k, j)
            # val = get_node_obj_val(dd, k) + func[i-1](l)    #computes the head node value based on the tail and arc label
            val = node_obj_vals[k] + func[i-1](l)    #computes the head node value based on the tail and arc label
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

    lp = []
    p = terminal_ind
    #computing the arc labels on the longest path through a backward search
    for i=n:-1:1
        # push!(lp, get_node_obj_label(dd, p))
        push!(lp, node_obj_labels[p])
        # p = get_node_obj_path(dd, p)
        p = node_obj_paths[p]
    end

    return reverse(lp), lpv

end



"""
Funtion: Second method for function Longest_Path!, where the separable functions are linear for all variables.
Input: The DD, and an array of coefficient of the variables in the linear objective function
Output: The longest path (composed of arc labels), and its length
"""
function longest_path!(dd::DecisionDiagram, cf::Array{T} where T<:Real)
    func = Array{Function}(length(cf))
    for (i, c) in enumerate(cf)                 #creates an array of separable functions for each variable
        func[i] = (x::Int64) -> c*x
    end
    return longest_path!(dd, func)
end
