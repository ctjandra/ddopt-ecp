include("../dd/dd.jl")

#Pkg.add("Clp")     #An open source LP solver
#Pkg.add("JuMP")    #A modeling interface
#Pkg.add("CPLEX")

using JuMP, Clp, CPLEX

function CGLP(dd::DecisionDiagram, fp::NTuple{N, Number} where N)

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
        lbl = get_arc_label(dd, e)  #gets a vector of possible labels for the arc (in case of parallel edges)
        for l in lbl
            @constraint(m, ProjCone[e,l], theta_plus[t_e] - theta_minus[t_e] - theta_plus[h_e] + theta_minus[h_e] + (gamma_plus[lyr] - gamma_minus[lyr])*l <= 0)
        end
    end

    #adding normalization constraints
    @constraint(m, NormCon, sum(gamma_plus[i] + gamma_minus[i] for i in 1:n) <= 1)

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
