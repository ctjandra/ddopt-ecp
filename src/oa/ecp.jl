include("../dd/dd.jl")
include("../cut/cutgenerator.jl")

"""
ECP: function to perform ECP with the DD cuts

Input:
n: Number of variables
m: The current model to which the cuts will be added
dd_list: An array of decision diagrams corresponding to inequalities
con_list: An array of the functions defining the inequalities
rhs_list: An array of the rhs value sof inequalities

Output:
m_output: The model after addition of cuts
x_output: Optimal solution
obj_output: Optimal value
status_output: Status of the output model
"""

function ECP(n::Int64, m::JuMP.Model, dd_list::Array{DecisionDiagram}, con_list::Array{Function}, rhs_list::Array{T} where T<:Real)

    #stopping criteria for the ECP algorithm
    iter_max_ecp = 100      #the maximum number of iterations to execute the algorithm
    tolerance_ecp = 0.05    #the tolerance for relative error
    tol_num_ecp = 5         #the number of past objective values wrt which the tolerance is computed. This number includes the current objective value
    obj_history_ecp = zeros(tol_num_ecp)    #the vector that stores the previous (improved) objective values for tolerance check

    #function to define what stopping criterion must be used:
    #0: iteration number
    #1: objective improvement: the relative difference between k consequtive improvements
    #2: error on the optimality gap (TO DO: need a primal heuristic method)
    function stop_rule_ecp(stp::Array{Int64})   #multiple stopping criteria can be passed as an array
        p = true
        for v in stp
            if v == 0       #uses the number of iteration criterion
                p = p && iter_ecp <= iter_max_ecp
            end
            if v == 1       #uses the objective improvement criterion
                if cut_added   #if a cut has been added to the model, the obj_history vector must be updated and the criterion is checked
                    push!(obj_history_ecp, obj_val)
                    deleteat!(obj_history_ecp, 1)
                    p = p && (obj_val - obj_history_ecp[1])/(abs(obj_val) + 1e-5) > tolerance_ecp
                end
            end
        end
        return p    #is false when the stopping criterion is met
    end

    #*******************************************************************************
    m2 = copy(m)
    x = getindex(m2,:x)     #to be able to use the same variables for the new model


    #initialization
    cut_technique = 1       #1:subgradient, 2:CGLP
    max_cut_num = 1     #maximum number of cuts (corresponding to violated inequalities) to be added at each iteration
    subgrad_start_point = zeros(con_num,n)  #stores the starting point for each constraint to be used in the subgradient method
    st_rule = [1, 2]        #determines the stop rule
    cut_added = false       #indicates if a cut (any type) has been added to the model
    con_num = length(rhs_list)          #the number of constraint
    violation_val = Array{Float64}(con_num)     #stores the violation value for each constraint at the current point
    iter_ecp = 1                                #iteration counter

    #getting the initial solution of the model
    status = solve(m2)
    obj_val = getobjectivevalue(m2)
    x_val = getvalue(x)
    m_output = m2
    x_output = x_val
    obj_output = obj_val
    status_output = :NOToptimal

    #************************************************************
    #performing the ecp method until the stop criteria is met
    while stop_rule_ecp(st_rule)

        #computing the violation value of the current point at each constraint
        for i=1:con_num
            violation_val[i] = con_list[i](x_val) - rhs_list[i]
        end

        sort_ind = sortperm(violation_val)     #stores the indices of the ascending sorted elements of violation_val

        if violation_val[sort_ind[end]] <= 0
            println("Optimal solution found!")
            status_output = :Optimal
            break
        end

        cut_num = 0       #number of cuts added so far
        CGLP_cut = false    #CGLP cut has been added?
        subgrad_cut = false     #subgradient cut has been added
        cut_added = false      #any cut has been added

        for i=con_num:-1:1

            if (added_cut < cut_num) || (violation_val[sort_ind[i]] <= 0)   #stop the loop if the inequality is satisfied or the number of allowable cuts is fulfilled
                break
            end

            if cut_technique == 1       #uses subgradient method to generate cut
                obj1, cut_coef1, cut_rhs1 = subgradient(dd_list[sort_ind[i]], x_val; sp = subgrad_start_point[sort_ind[i],:])

                if obj1 > 0             #if the inequality cuts off the point
                    cut_num += 1
                    subgrad_cut = true
                    aff = AffExpr(x, cut_coef1, 0)      #stores the affine expression of the inequality
                    @constraint(m2, aff <= cut_rhs1)
                    subgrad_start_point[sort_ind[i],:] = cut_coef1      #uses the current cut's normal vector as the initial point for the next subgrad iteration
                end
            end


            if cut_technique == 2       #uses CGLP to generate cut
                obj1, cut_coef1, cut_rhs1, status1 = CGLP(dd_list[sort_ind[i]], x_val)

                if (status1 == :Optimal) && (obj1 > 0)
                    cut_num += 1
                    CGLP_cut = true
                    aff = AffExpr(x, cut_coef1, 0)      #stores the affine expression of the inequality
                    @constraint(m2, aff <= cut_rhs1)
                end
            end

        end

        cut_added = CGLP_cut || subgrad_cut

        #getting the initial solution of the model if a cut has been added
        if cut_added
            status = solve(m2)
            obj_val = getobjectivevalue(m2)
            x_val = getvalue(x)
            m_output = m2
            x_output = x_val
            obj_output = obj_val
            status_output = :NOToptimal
        end

        iter_ecp += 1

    end

    return m_output, x_output, obj_output, status_output

end
