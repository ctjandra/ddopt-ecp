# This file includes the structure of Branch-and-Bound tree.

include_dependency("../dd/core.jl")
include_dependency("../dd/ordering.jl")
include_dependency("../dd/mathprog/dp_specification.jl")
include_dependency("../dd/mathprog/function_analysis.jl")


"""
Type to store node information in b&b
"""
struct NodeInfo
    parent_id::Int                  # id of the parent node in the tree (array)
    children_id::Tuple{Int, Int}      # id of the children nodes in the tree
    depth::Int                         # depth of the node in the tree
end


"""
Type to store model information in b&b
"""
struct ModelInfo
    m::JuMP.Model                       # the model after adding cuts in the node
    #    best_primal_sol::Array{Float64, 1}          # the best feasible solution at this node
    #    obj_primal::Float64                         # objective lower (primal) bound obtained at the node
    #    best_dual_sol::Array{Float64, 1}           # the optimal solution of the relaxation at this node
    #    obj_dual::Float64                         # objective upper (dual) bound obtained at the node
end



"""
Type to store DD information in b&b
"""
struct DDInfo
    dd::Array{DecisionDiagram,1}        # array of dds of all constraints in the model
    initial::Array{Float64, 1}              # array of intial states for constraints
    transition::Array{Array{UVFunc, 1}, 1}     # 1st dim: constraints, 2nd dim: variables, contenct: transtiion function
    con_function::Array{Function, 1}            # stores the function representing variable part of constraints
    domain::Array{Tuple{Float64, Float64, Symbol}, 1}   # array of variable domains
    order::Array{Ordering, 1}                   # ordering of variables in dd, for constraints
    lb_ub::Array{Array{Tuple{Float64, Float64}, 1},1}   # 1st dim: constraint, 2nd dim: variable, content: (lowerbound, upperbound) of the univariate function over the variable range defined in specs.
end


"""
Type to store Branching information in b&b
"""
struct BranchInfo
    #    branch_policy::?                   # which of the two policies is used in this branch form the parent node: variable range (arc label) branch, constraint value (node state) branch, and on which variables or constraints?
    #    refine_policy::Array{Bool, 1}        # do we reconstruct the dds in this branch (having more width) or we just leave it as it is after branching
    #    modified_node_layers::Array{Array{Array{Bool, 1}, 1}, 1}        # 1st dim: constraints, 2nd dim and 3rd dim: dd node layers, content: whether the node is removed (0) or not from the dd of the parent node
    #    modified_arc_layers::Array{Array{Array{Bool, 1}, 1}, 1}     # similar to the modified_node_layers, but for arcs
end



"""
Type to store branch-and-bound nodes
"""
struct BBNode
    node_info::NodeInfo
    model_info::ModelInfo
    dd_info::DDInfo
    branch_info::BranchInfo
end




"""
Type to store branch-and-bound tree
"""
struct BBTree
    original_model::JuMP.Model      # the main model of the problem
    LP_model::JuMP.Model            # the LP model used in the ECP
    tree::Array{BBNode, 1}          # an array that stores branch-and-bound nodes (works as a tree)
end




"""
Create a BB tree from a model

# Input
- `m1`: Original model
- `m2`: LP relaxation of the model used in the ECP
"""
function create_BB_tree(m1::JuMP.Model, m2::JuMP.Model)::BBTree
    bb_tree = BBTree(m1, m2, Array{BBNode, 1}())
    return bb_tree
end



"""
create root node in branch-and-bound

# Input
- `bb_tree`: the newly constructed object of type BBTree, which contains a model and its tree array field must be empty.
- `prop`: Keyword that contains an initial assessment for properties of the univariate functions (optional)
"""
function create_BB_root!(bb_tree::BBTree; prop::Array{Array{Array{UVFuncProp, 1}, 1}, 1} = [])
    @assert(length(bb_tree.tree) == 0)          # ensure that the tree array field must be empty

    n::Int = MathProgBase.numvar(bb_tree.original_model)      # number of variables in model
    c_n::Int = JuMP.numnlconstr(bb_tree.original_model) + MathProgBase.numquadconstr(bb_tree.original_model)      # number of nonlinear constraints in model

    node_info = NodeInfo(0, (0,0), 0)
    branch_info = BranchInfo()
    model_info = ModelInfo(copy(bb_tree.LP_model))
    dd_info = DDInfo(Array{DecisionDiagram,1}(  c_n), Array{Float64, 1}(  c_n), Array{Array{UVFunc, 1}, 1}(  c_n), Array{Function, 1}(  c_n), Array{Tuple{Float64, Float64, Symbol}, 1}(  n), Array{Ordering, 1}(  c_n), Array{Array{Tuple{Float64, Float64}, 1},1}(  c_n))

    # Determine constraint decomposition, ordering and variable domains
    eval = evaluate_separable_constraints(bb_tree.original_model)
    create_constraint_domain_function!(dd_info.domain, bb_tree.original_model)
    for i=1:c_n
        dd_info.order[i] = NoOrdering()
    end

    # Create specs for decision diagram
    for i::Int=1:c_n
        dd_info.initial[i] = create_constraint_initial_state(eval[i])
        dd_info.con_function[i] = create_constraint_single_function(eval[i])

        dd_info.transition[i] = Array{UVFunc, 1}(n)
        if isassigned(prop, i)
            create_constraint_transition_function!(dd_info.transition[i], eval[i], prp = prop[i])            # if functional property is given as input
        else
            create_constraint_transition_function!(dd_info.transition[i], eval[i])                 # no functional properties is given as input
        end

        dd_info.lb_ub[i] = Array{Tuple{Float64, Float64}, 1}(  n)
        for j::Int=1:n
            if isassigned(dd_info.transition[i], j)
                dd_info.lb_ub[i][j] = (uvf_lb(dd_info.transition[i][j].f, dd_info.transition[i][j].prop, dd_info.domain[j]), uvf_ub(dd_info.transition[i][j].f, dd_info.transition[i][j].prop, dd_info.domain[j]))
            end
        end

        println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        # Construct decision diagram
        dd_info.dd[i] = construct_DD(dd_info.transition[i], dd_info.initial[i], dd_info.domain, dd_info.lb_ub[i], ordering=dd_info.order[i], reduced_arc=true)
        println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
    end

    push!(bb_tree.tree, BBNode(node_info, model_info, dd_info, branch_info))

end
