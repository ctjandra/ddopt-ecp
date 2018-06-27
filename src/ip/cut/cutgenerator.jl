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
function CGLP(dd::DecisionDiagram, fractional_point::Array{Float64}; tilim::Float64 = 100000)

    n = nvars(dd)     #the number of arc layers of dd (number of variables)
    @assert(n == length(fractional_point))       #makes sure the size of the fractional point matches the dimension of variables represented by dd

    node_num = nnodes(dd)
    arc_num = narcs(dd)

    # Preparing an optimization model
    m = Model(solver = CplexSolver(CPX_PARAM_TILIM = tilim))

    # Defining variables
    @variable(m, theta_plus[1:node_num] >= 0)
    @variable(m, theta_minus[1:node_num] >= 0)
    @variable(m, gamma_plus[1:n] >= 0)
    @variable(m, gamma_minus[1:n] >= 0)

    # Fixing the theta of the source node at zero
    source_ind = root(dd)            #gets the index value of the source node
    setupperbound(theta_plus[source_ind], 0)
    setupperbound(theta_minus[source_ind], 0)

    # Getting the index of the terminal node
    terminal_ind = terminal(dd)

    # Defining objective function
    @objective(m, Max, sum((gamma_plus[i] - gamma_minus[i])*fractional_point[i] for i in 1:n) - theta_plus[terminal_ind] + theta_minus[terminal_ind])

    # Adding projection cone constraints via accessing the inner graph information
    #g = dd.graph
    #for e in edges(g)
    #    t_e = src(e)    #gets the tail of e
    #    h_e = dst(e)    #gets the head of e
    #    lyr = get_arc_layer(dd, e)    #gets the layer of e
    #    lbl = get_arc_labels(dd, e)  #gets a vector of all labels for the arc (in case of parallel arcs)
    #    for l in lbl
    #        @constraint(m, theta_plus[t_e] - theta_minus[t_e] - theta_plus[h_e] + theta_minus[h_e] + (gamma_plus[lyr] - gamma_minus[lyr])*l <= 0)
    #    end
    #end

    # Adding projection cone constraints without accessing the inner graph information
    for i=2:n+1, node in dd.layers[i]
        for (parent, label) in inneighbors(dd, node)
            @constraint(m, theta_plus[parent] - theta_minus[parent] - theta_plus[node] + theta_minus[node] + (gamma_plus[i-1] - gamma_minus[i-1])*label <= 0)
        end
    end


    # Adding normalization constraints
    @constraint(m, sum(gamma_plus[i] + gamma_minus[i] for i in 1:n) <= 1)

    # Solving the model
    status = solve(m)

    # Getting the output
    obj = getobjectivevalue(m)  #gets the optimal value

    println("--------------")
    println("status: ", status)
    println("CGLP obj: ", obj)

    coefficients = zeros(n)   #coefficient of the CGLP inequality
    for i=1:n
        coefficients[i] = getvalue(gamma_plus[i]) - getvalue(gamma_minus[i])
    end
    rhs = getvalue(theta_plus[terminal_ind]) - getvalue(theta_minus[terminal_ind])  #gets the rhs value for the CGLP inequality

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
function subgradient(dd::DecisionDiagram, fractional_point::Array{Float64}; step_rule::Int64 = 0, starting_point::Array{Float64} = zeros(length(fractional_point)), tilim::Float64 = 100000)

    n = nvars(dd)     #the number of arc layers of dd
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

    if step_rule == 0
        c0 = 1
        rho = (t::Int64,norm_g_t::Float64) -> c0
    elseif step_rule == 1
        c0 = 1
        rho = (t::Int64,norm_g_t::Float64) -> c0/norm_g_t
    elseif step_rule == 2
        c0 = 1
        rho = (t::Int64,norm_g_t::Float64) -> c0/sqrt(t)
    elseif step_rule == 3
        c0 = 1
        a0 = 0
        rho = (t::Int64,norm_g_t::Float64) -> c0/(a0 + t)
    end


    # initialization
    best_obj = 0            # best objecive found so far
    best_coefficients = zeros(n)      # coefficient vector corresponding to the best objective
    best_rhs = 0            # right-hand-side value corresponding to the best objective
    delta = 0               # the violation value (objective value)
    subgrad = ones(n)       # subgradient vector
    iter = 1                # iteration number
    stop_criterion = [0, 1]                # stopping criterion
    if norm(starting_point) > 0         # normalizing the starting point and making a copy
        gamma = starting_point/norm(starting_point)
    else
        gamma = copy(starting_point)
    end

    # stopping criteria for the entire algorithm
    iter_max = 20      # the maximum number of iterations to execute the algorithm
    tolerance = 0.05    # the tolerance for relative error
    tol_num = 2         # the number of past objective values (including the current one) wrt which the tolerance is computed
    obj_history = Array{Float64}(tol_num)    # stores the previous (improved) objective values for tolerance check

    # function to define what stopping criterion must be used:
    # 0: time
    # 1: iteration number
    # 2: concavity error: obtained from the first order concavity inequality over the unit ball normalization constraint: at iteration t, the objective error is <= norm of subgradient - violation value
    # 3: objective improvement: the relative difference between k consequtive improvements
    function stop_rule(st::Array{Int64})   # multiple stopping criteria can be passed as an array
        p = true
        for v in st
            if v == 0         # compares time elapsed with given time limit
                clock_end = time()
                time_elapsed = float(clock_end - clock_start)
                p = p && time_elapsed <= tilim
            end
            if v == 1       # uses the number of iteration criterion
                p = p && iter <= iter_max
            end
            if v == 2       # uses the concavity error criterion
                p = p && (norm(subgrad) - delta)/max(norm(subgrad), 1e-9) > tolerance
            end
            if v == 3       # uses the objective improvement criterion (ONLY considers positive violation values)
                if obj_history[end] != best_obj   # if the objective value is updated, the obj_history vector must be updated as well and the criterion is checked
                    push!(obj_history, best_obj)
                    deleteat!(obj_history, 1)
                    p = p && (best_obj - obj_history[1])/best_obj > tolerance
                end
            end
        end
        return p    # is false when the stopping criterion is met
    end

    # initializing time parameter
    clock_start = time()

    # Main algorithm
    while stop_rule(stop_criterion)

        # computing the longest path with the objective function defined by the current direction
        lp1, lpv1 = longest_path(dd, gamma)

        # computing the violatation value of the current inequality at the given fractional point
        delta = gamma'*fractional_point - lpv1

        if delta > best_obj                 # since best_obj >= 0 by construction, this means that a cut is detected
            if best_obj <= 0                # i.e., the current inequality is the first that cuts off the fractional point
                iter_max = iter + 5                   # reduce the number of iterations after detecting a cut
                stop_criterion= [0, 1, 2]             # changes the stop rule to the one that considers both iteration number and concavity error
            end
            best_obj = delta
            best_coefficients = gamma
            best_rhs = lpv1
        end

        # computing the subgradient at current iteation
        subgrad = fractional_point - lp1

        # updating the coefficient vector
        phi = gamma + rho(iter,norm(subgrad))*subgrad

        # projecting the updated vector onto the unit ball normalization cosntraint
        if norm(phi) > 1
            gamma = phi/norm(phi)
        else
            gamma = phi
        end

        iter += 1

    end

    return best_obj, best_coefficients, best_rhs

end
