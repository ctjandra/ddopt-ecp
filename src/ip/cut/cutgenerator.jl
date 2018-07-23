include_dependency("../../dd/construction.jl")
include_dependency("../../dd/modification.jl")
include_dependency("../../dd/computation.jl")

using JuMP, Clp, CPLEX

"""
Generate CGLP cuts.

# Input
- `dd::DecisionDiagram`: DD.
- `fractional_point::Array{Float64}`: The point to be separated (as a vector).
- `tilim::Float64`: Time limit in seconds (optional).

# Output
- `obj`: The optimal value of the CGLP.
- `coefficients`: The coefficient vector of the cut.
- `rhs`: The right-hand-side of the cut.
- `status`: The status of the CGLP.
"""
function CGLP(dd::DecisionDiagram, fractional_point::Array{Float64, 1}; tilim::Float64 = 100000)::Tuple{Float64, Array{Float64, 1}, Float64, Symbol}

    n::Int = nvars(dd)     #the number of arc layers of dd (number of variables)
    @assert(n == length(fractional_point))       #makes sure the size of the fractional point matches the dimension of variables represented by dd

    node_num::Int = nnodes(dd)
    arc_num::Int = narcs(dd)

    C_layersize::Array{Int, 1} = Array{Int, 1}(0)     # compute the cumulative number of nodes in previous layers of the indexed layer
    push!(C_layersize, 0)
    for i::Int = 2::n+1
        push!(C_layersize, C_layersize[i-1] + get_node_layer_size(dd, i-1))
    end

    # Preparing an optimization model
    m = Model(solver = CplexSolver(CPX_PARAM_TILIM = tilim))

    # Defining variables
    # NOTE: Jump does not support Array{Array} of variables in the definition
    # There are two options:
    # 1: Use anonymous variable definition which has a few restrictions (cannot get index for a copy model)
    # 2: To define "flattend out" theta variables, by numbering them manually
    @variable(m, theta_plus[1:node_num] >= 0)
    @variable(m, theta_minus[1:node_num] >= 0)
    @variable(m, gamma_plus[1:n] >= 0)
    @variable(m, gamma_minus[1:n] >= 0)

    # Fixing the theta of the source node at zero
    setupperbound(theta_plus[1], 0)
    setupperbound(theta_minus[1], 0)

    # Defining objective function
    @objective(m, Max, sum((gamma_plus[i] - gamma_minus[i])*fractional_point[i] for i in 1:n) - theta_plus[node_num] + theta_minus[node_num])

    # Adding projection cone constraints without accessing the inner graph information
    for i::Int=1:n
        for j::int=1:get_arc_layer_size(dd, i)
            @constraint(m, theta_plus[get_arc_tail(dd, i, j) + C_layersize[i]] - theta_minus[get_arc_tail(dd, i, j) + C_layersize[i]] - theta_plus[get_arc_head(dd, i, j) + C_layersize[i+1]] + theta_minus[get_arc_head(dd, i, j) + C_layersize[i+1]] + (gamma_plus[i] - gamma_minus[i])*get_arc_label(dd, i, j) <= 0)
        end
    end


    # Adding normalization constraints
    @constraint(m, sum(gamma_plus[i] + gamma_minus[i] for i in 1:n) <= 1)

    # Solving the model
    status::Symbol = solve(m)

    # Getting the output
    obj::Float64 = getobjectivevalue(m)  #gets the optimal value

    println("--------------")
    println("status: ", status)
    println("CGLP obj: ", obj)

    coefficients::Array{FLoat64, 1} = zeros(n)   #coefficient of the CGLP inequality
    for i::Int=1:n
        coefficients[i] = getvalue(gamma_plus[i]) - getvalue(gamma_minus[i])
    end
    rhs::Float64 = getvalue(theta_plus[node_num]) - getvalue(theta_minus[node_num])  #gets the rhs value for the CGLP inequality

    return obj, coefficients, rhs, status
end



"""
Generate cuts via the subgradient method.

# Input
- `dd::DecisionDiagram`: DD.
- `fractional_point::Array{Float64}`: The point to be separated (as a vector).
- `step_rule::Int64`: Step size rule (optional).
- `starting_point::Array{Float64}`: Starting point (optional).
- `tilim::Float64`: Time limit in seconds (optional).

# Output
- `best_obj`: The optimal value of the subgradient method at termination.
- `best_coefficients`: The coefficient vector of the best inequality at termination.
- `best_rhs`: The right-hand-side of the best inequality at termination.
"""
function subgradient(dd::DecisionDiagram, fractional_point::Array{Float64, 1}; step_rule::Int64 = 0, starting_point::Array{Float64, 1} = zeros(length(fractional_point)), tilim::Float64 = 100000.0)::Tuple{Float64, Array{Float64, 1}, Float64}

    n::Int = nvars(dd)     #the number of arc layers of dd
    @assert(n == length(fractional_point))       #makes sure the size of the fractional point matches the dimension of variables represented by dd


    # Choosing the step rule
    # 0. constant step size: \rho_t = c for some c>0
    # 1. constraint step length: \rho_t = c/||g_t||, where c>0 and g_t is the subgradient at step t
    # 2. diminishing step rule: \rho_t = c/sqrt(t), where c>0
    # 3. square summable diminishing step rule: \rho_t = c/(a+t), where c>0 and a>=0
    if step_rule < 0 || step_rule > 3
        step_rule = 0
        info("Mismatch in the step size rule. The default rule has been used.")
    end

    c0::Int = 1
    a0::Int = 0
    if step_rule == 0
        rho = (t::Int64,norm_g_t::Float64) -> c0
    elseif step_rule == 1
        rho = (t::Int64,norm_g_t::Float64) -> c0/norm_g_t
    elseif step_rule == 2
        rho = (t::Int64,norm_g_t::Float64) -> c0/sqrt(t)
    elseif step_rule == 3
        rho = (t::Int64,norm_g_t::Float64) -> c0/(a0 + t)
    end


    # initialization
    best_obj::Float64 = 0            # best objecive found so far
    best_coefficients::Array{Float64, 1} = zeros(n)      # coefficient vector corresponding to the best objective
    best_rhs::Float64 = 0            # right-hand-side value corresponding to the best objective
    delta::Float64 = 0               # the violation value (objective value)
    subgrad::Array{Float64, 1} = ones(n)       # subgradient vector
    iter::Int = 1                # iteration number
    stop_criterion::Array{Int, 1} = [0, 1]                # stopping criterion
    if norm(starting_point) > 0         # normalizing the starting point and making a copy
        gamma::Array{Float64, 1} = starting_point/norm(starting_point)
    else
        gamma = copy(starting_point)
    end

    # stopping criteria for the entire algorithm
    iter_max::Int = 20      # the maximum number of iterations to execute the algorithm
    tolerance::Float64 = 0.0005    # the tolerance for relative error
    tol_num::Int = 2         # the number of previous objective values (including the current one) wrt which the tolerance is computed
    obj_history::Array{Float64, 1} = Array{Float64, 1}(tol_num)    # stores the previous (improved) objective values for tolerance check

    # function to define what stopping criterion must be used:
    # 0: time
    # 1: iteration number
    # 2: concavity error: obtained from the first order concavity inequality over the unit ball normalization constraint: at iteration t, the objective error is <= norm of subgradient - violation value
    # 3: objective improvement: the relative difference between k consequtive improvements
    function stop_rule(st::Array{Int64})::Bool   # multiple stopping criteria can be passed as an array
        p = true
        for v::Int in st
            if v == 0         # compares time elapsed with given time limit
                clock_end::Float64 = time()
                time_elapsed::Float64 = float(clock_end - clock_start)
                p = p && (time_elapsed <= tilim)
            end
            if v == 1       # uses the number of iteration criterion
                p = p && (iter <= iter_max)
            end
            if v == 2       # uses the concavity error criterion
                p = p && ((norm(subgrad) - delta)/max(norm(subgrad), 1e-9) > tolerance)
            end
            if v == 3       # uses the objective improvement criterion (ONLY considers positive violation values)
                if obj_history[end] != best_obj   # if the objective value is updated, the obj_history vector must be updated as well and the criterion is checked
                    push!(obj_history, best_obj)
                    deleteat!(obj_history, 1)
                    p = p && ((best_obj - obj_history[1])/best_obj > tolerance)
                end
            end
        end
        return p    # is false when the stopping criterion is met
    end

    # initializing time parameter
    clock_start::Float64 = time()

    # Main algorithm
    while stop_rule(stop_criterion)

        # computing the longest path with the objective function defined by the current direction
        lp1::Array{Float64, 1}, lpv1::Float64 = @time longest_path(dd, gamma)

        # computing the violatation value of the current inequality at the given fractional point
        delta = gamma'*fractional_point - lpv1

        if delta > best_obj                 # since best_obj >= 0 by construction, this means that a cut is detected
            if best_obj <= 0                # i.e., the current inequality is the first that cuts off the fractional point
                iter_max = iter + 10                   # reduce the number of iterations after detecting a cut
                stop_criterion = [0, 1]             # changes the stop rule to the one that considers both iteration number and concavity error
            end
            best_obj = delta
            best_coefficients = gamma
            best_rhs = lpv1
        end

        # computing the subgradient at current iteation
        subgrad = fractional_point - lp1

        # updating the coefficient vector
        phi::Array{Float64, 1} = gamma + rho(iter,norm(subgrad))*subgrad

        # projecting the updated vector onto the unit ball normalization cosntraint
        #if norm(phi) > 1
        #    gamma = phi./norm(phi)
        #else
            gamma = phi
        #end

        iter += 1

    end

    return best_obj, best_coefficients, best_rhs

end
