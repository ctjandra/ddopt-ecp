include("../dd/dd.jl")
include("../dd/functions.jl")

#Pkg.add("Clp")     #An open source LP solver
#Pkg.add("JuMP")    #A modeling interface
#Pkg.add("CPLEX")

using JuMP, Clp, CPLEX

"""
CGLP: function to generate CGLP cuts
Input: DD, and the point to be separated (as a vector)
Output: the optimal value of the CGLP, the coefficient vector and the right-hand-side of the resulting inequality, and the status of the CGLP
"""
function CGLP(dd::DecisionDiagram, fp::Array{Float64,1})

    n = size(dd.layers, 1) - 1     #the number of arc layers of dd
    @assert(n == length(fp))       #makes sure the size of the fractional point matches the dimension of variables represented by dd

    g = dd.graph
    node_num = nv(g)
    arc_num = ne(g)

    #Preparing an optimization model
    m = Model(solver = CplexSolver())

    #deining variables
    @variable(m, theta_plus[vertices(g)] >= 0)
    @variable(m, theta_minus[vertices(g)] >= 0)
    @variable(m, gamma_plus[1:n] >= 0)
    @variable(m, gamma_minus[1:n] >= 0)

    #fixing the theta of the source node at zero
    source_ind = dd.layers[1][1]            #gets the index value of the source node
    setupperbound(theta_plus[source_ind], 0)
    setupperbound(theta_minus[source_ind], 0)

    #getting the index of the terminal node
    terminal_ind = dd.layers[end][1]

    #defining objective function
    @objective(m, Max, sum((gamma_plus[i] - gamma_minus[i])*fp[i] for i in 1:n) - theta_plus[terminal_ind] + theta_minus[terminal_ind])

    #adding projection cone constraints
    for e in edges(g)
        t_e = src(e)    #gets the tail of e
        h_e = dst(e)    #gets the head of e
        lyr = get_arc_layer(dd, e)  #gets the layer of e
        lbl = get_arc_label(dd, e)  #gets a vector of possible labels for the arc (in case of parallel arcs)
        for l in lbl
            @constraint(m, theta_plus[t_e] - theta_minus[t_e] - theta_plus[h_e] + theta_minus[h_e] + (gamma_plus[lyr] - gamma_minus[lyr])*l <= 0)
        end
    end

    #adding normalization constraints
    @constraint(m, sum(gamma_plus[i] + gamma_minus[i] for i in 1:n) <= 1)

    #solving the model
    status = solve(m)

    #getting the output
    obj = getobjectivevalue(m)  #gets the optimal value

    cf = zeros(n)   #coefficient of the CGLP inequality
    for i=1:n
        cf[i] = getvalue(gamma_plus[i]) - getvalue(gamma_minus[i])
    end
    rhs = getvalue(theta_plus[terminal_ind]) - getvalue(theta_minus[terminal_ind])  #gets the rhs value for the CGLP inequality

    return obj, cf, rhs, status
end



"""
Subgradient: function to generate cuts via the subgradient method
Input: DD, the point to be separated (as a vector); step size rule (optional), starting point (optional)
Output: the optimal value of the CGLP, the coefficient vector and the right-hand-side of the resulting inequality
"""
function subgradient!(dd::DecisionDiagram, fp::Array{Float64}; step_rule::Int64 = 0, sp::Array{Float64} = zeros(length(fp)))

    n = size(dd.layers, 1) - 1     #the number of arc layers of dd
    @assert(n == length(fp))       #makes sure the size of the fractional point matches the dimension of variables represented by dd


    #Choosing the step rule
    #0. constant step size: \rho_t = c for some c>0
    #1. constraint step length: \rho_t = c/||g_t||, where c>0 and g_t is the subgradient at step t
    #2. diminishing step rule: \rho_t = c/sqrt(t), where c>0
    #3. square summable diminishing step rule: \rho_t = c/(a+t), where c>0 and a>=0
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


    #initialization
    best_obj = 0            #best objecive found so far
    best_cf = zeros(n)      #coefficient vector corresponding to the best objective
    best_rhs = 0            #right-hand-side value corresponding to the best objective
    iter = 1                #iteration number
    st = 0                  #stopping criterion
    if norm(sp) > 0         #normalizing the starting point and making a copy
        gamma = sp/norm(sp)
    else
        gamma = copy(sp)
    end

    #stopping criteria
    iter_max = 100      #the maximum number of iterations to execute the algorithm
    tolerance = 0.05    #the tolerance for relative objective improvement
    tol_num = 2         #the number of past objective values wrt which the tolerance is computed. This number includes the current objective value
    obj_history = [0.0 for i=1:tol_num]    #the vector that stores the previous (improved) objective values for tolerance check

    #function to define what stopping criterion must be used: 0: only iteration number, 1: only objective improvement, 2: both
    function stop_rule(st::Int64)
        p = true
        if st == 0 || st == 2       #uses the number of iteration criterion
            p = iter <= iter_max
        end
        if st == 1 || st == 2       #uses the objective improvement criterion
            if endof(obj_history) != best_obj   #if the objective value is updated, the obj_history vector must be updated as well and the criterion is checked
                push!(obj_history, best_obj)
                deleteat!(obj_history, 1)
                p = p && (best_obj - obj_history[1])/best_obj > tolerance
            end
        end
        return p    #is false when the stopping criterion is met
    end


    #Main algorithm
    while stop_rule(st)

        #computing the longest path with the objective function defined by the current direction
        lp1, lpv1 = longest_path!(dd, gamma)

        #computing the violatation value of the current inequality at the given fractional point
        delta = gamma'*fp - lpv1

        if delta > best_obj
            best_obj = delta
            best_cf = gamma
            best_rhs = lpv1
            st = 2          #changes the stop rule to the one that considers both iteration number and objective improvement
        end

        #computing the subgradient at current iteation
        subgrad = fp - lp1

        #updating the coefficient vector
        phi = gamma + rho(iter,norm(subgrad))*subgrad

        #projecting the updated vector onto the unit ball normalization cosntraint
        if norm(phi) > 1
            gamma = phi/norm(phi)
        else
            gamma = phi
        end

        iter += 1

    end

    return best_obj, best_cf, best_rhs

end
