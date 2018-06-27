include("../../dd/construction.jl")
include("../cut/cutgenerator.jl")

"""
Perform ECP using DD cuts.

# Input
- `n::Int64`: Number of variables.
- `m::JuMP.Model`: The current model to which the cuts will be added.
- `dd_list::Array{DecisionDiagram}`: An array of decision diagrams corresponding to inequalities of the current model.
- `con_list::Array{SeparableFunctionEvaluation}`: An array of constraint specifications with one-to-one correspondence to the DD array.
- `time_limit::Float64`: Time limit to run the algorithm (optional).


# Output
- `m_output`: The model after the addition of cuts.
- `x_output`: Optimal solution.
- `obj_output`: Optimal value.
- `status_output`: Status of the output model.
"""

function ECP(n::Int64, m::JuMP.Model, dd_list::Array{DecisionDiagram}, con_list::Array{SeparableFunctionEvaluation}; time_limit = 10000)

    @assert (length(dd_list) == length(con_list))

    # stopping criteria for the ECP algorithm
    iter_max_ecp = 100      # the maximum number of iterations to execute the algorithm
    tolerance_ecp = 0.05    # the tolerance for relative error
    tol_num_ecp = 5         # the number of past objective values (including the current one) wrt which the tolerance is computed.
    obj_history_ecp = [-Inf for i=1:tol_num_ecp]    # the vector that stores the previous (improved) objective values for tolerance check


    # function to define what stopping criterion must be used:
    # 0: time limit
    # 1: iteration number
    # 2: objective improvement: the relative difference between k consequtive improvements
    # 3: error on the optimality gap (TO DO: need a primal heuristic method)
    function stop_rule_ecp(stp::Array{Int64})   # multiple stopping criteria can be passed as an array
        p = true
        for v in stp
            if v ==  0      # time limit criterion
                clock_end = time()
                time_elapsed = float(clock_end - clock_start)
                p = p && time_elapsed <= time_limit
                #println("v==0: ", p)
            end
            if v == 1       # the number of iteration criterion
                p = p && iter_ecp <= iter_max_ecp
                #println("v==1: ", p)
            end
            if v == 2       # uses the objective improvement criterion
                if cut_added   # if a cut has been added to the model, the `obj_history_ecp` vector must be updated and the criterion is checked
                    push!(obj_history_ecp, obj_val)
                    deleteat!(obj_history_ecp, 1)
                    p = p && (obj_val - obj_history_ecp[1])/(abs(obj_val) + 1e-5) > tolerance_ecp
                end
            end
        end
        return p    # is false when the stopping criterion is met
    end

    #*******************************************************************************
    m2 = copy(m)
    x = getindex(m2,:x)     # allows to use the same variables of the orignal model for the new model


    # initializing cut gnerating paratmeters
    subgrad_flag = true       # if subgradient method is used to generate cuts
    CGLP_flag = false             # if CGLP method is used to generate cuts
    use_CGLP = true          # determines whether CGLP will be used at all in the ECP
    max_cut_num = 10     # maximum number of cuts (one cut for each violated inequality) to be added at each iteration
    st_rule = [0, 1, 2]        # sets the stopping rule for the algorithm
    cut_added = false       # indicates whether a cut (of any type) has been added to the model
    con_num = length(con_list)          # the number of constraint
    subgrad_start_point = zeros(con_num,n)  # stores the starting point for each constraint to be used in the subgradient method
    violation_val = Array{Float64}(con_num)     # stores the violation value for each constraint at the current optimal point
    normalized_violation_val = Array{Float64}(con_num)     # stores the violation value after being normalized wrt the rhs
    iter_ecp = 1                                # iteration counter

    # set parameters to switch between solving LP and IP at iterations
    lp_flag = true          # if the LP method is used to solve the current model
    ip_flag = false         # if the IP method is used to solve the current model
    lp_to_ip_iter = 100       # number of iterations to solve as LP before switching to IP
    lp_to_ip_tol = 0.005     # relative constraint violation threshold in the LP mode before swithcing to IP

    # initializing time parameters
    clock_start = time()
    CplexSolver(CPX_PARAM_TILIM = time_limit)        #sets time limit parameter for the solver

    # get the initial solution of the model
    status = solve(m2)
    obj_val = getobjectivevalue(m2)
    x_val = getvalue(x)
    m_output = m2
    x_output = x_val
    obj_output = obj_val
    status_output = :Not_Optimal

    #************************************************************
    #performing the ecp method until the stopping criteria is met

    #println("@@@@@@@@@@@@@@@@@@@@", stop_rule_ecp(st_rule))

    while stop_rule_ecp(st_rule)

        #println("@@@@@@@@@@@@@@@@@@@@")

        # if the previous iteration has not added any cuts to the model:
        if cut_added == false && iter_ecp > 1
            # 1. if the algorithm is allowed to switch to CGLP
            if use_CGLP == true
                # 1.1. if CGLP method has not been activated: activate it
                if CGLP_flag == false
                    CGLP_flag == true
                    subgrad_flag == false
                # 1.2. if the CGLP has been solved to optimality: no cuts can be obtained from the current OA, B&B is needed to improve further
                elseif status_CGLP == :Optimal
                    println("OA (NOT global) optimal solution found!")
                    status_output = :OA_Optimal
                    break
                # 1.3. if the CGLP was not solved to optimality: cannot decide on the optimality of the current OA
                else
                    println("CGLP not optimal!")
                    status_output = :CGLP_Limit
                    break
                end
            # 2. the subgradient algorithm has not been able to detect a cut: cannot decide on the optimiality of the current OA
            else
                println("Subgradient method has not added cuts!")
                status_output = :Subgrad_Limit
                break
            end
        end


        #computing the violation value of the current point at each constraint
        for i=1:con_num
            violation_val[i] = con_list[i].single_function(x_val) + con_list[i].constant
            normalized_violation_val[i] = violation_val[i]/max(1e-9, abs(con_list[i].constant))
        end

        sort_ind = sortperm(normalized_violation_val)     #stores the indices of the ascending sorted elements of violation_val

        # if the current point satisfies all constraints, it is global optimal
        if violation_val[sort_ind[end]] <= 0
            println("Global optimal solution found!")
            status_output = :Global_Optimal
            break
        end

        # check whether to solve IP or LP based on the violation tolerance and the iteration number
        if lp_flag == true
            if (normalized_violation_val[sort_ind[end]] <= lp_to_ip_tol) || (iter_ecp > lp_to_ip_iter)
                lp_flag == false
                ip_flag == true
                for i=1:n
                    setcategory(x[i], :Int)        # makes variable type integer
                end
            end
        end

        cut_num = 0       # number of cuts added so far
        CGLP_cut = false    # CGLP cut has been added
        subgrad_cut = false     # subgradient cut has been added
        cut_added = false      # any cut has been added
        status_CGLP = :Optimal # it is assumed that CGLP will solve to optimality by default

        for i=con_num:-1:1

            clock_end = time()
            time_elapsed = float(clock_end - clock_start)
            c_index = sort_ind[i]

            # stop the loop if the inequality is satisfied or the number of allowable cuts is fulfilled or time exceeds the limit
            if (max_cut_num < cut_num) || (violation_val[c_index] <= 0) || (time_elapsed > time_limit)
                break
            end

            time_remain = max(0, time_limit - time_elapsed)     # compute the remaining time

            if subgrad_flag == true       # uses subgradient method to generate cut
                obj1, cut_coef1, cut_rhs1 = subgradient(dd_list[c_index], x_val, step_rule = 0, starting_point = subgrad_start_point[c_index,:], tilim = time_remain)

                #println("\nsubgrad alg used\n")

                if obj1 > 0             # if the inequality cuts off the point
                    cut_num += 1
                    subgrad_cut = true
                    aff = AffExpr(x, cut_coef1, 0)      # stores the affine expression of the inequality
                    @constraint(m2, aff <= cut_rhs1)
                    subgrad_start_point[c_index,:] = cut_coef1      # uses the current cut's normal vector as the initial point for the next subgrad iteration
                end
            end

            if (CGLP_flag == true) && (use_CGLP == true)       # uses CGLP to generate cut

                println("SOLVING CGLP NOW")

                obj1, cut_coef1, cut_rhs1, status1 = CGLP(dd_list[c_index], x_val, tilim = time_remain)

                if (status1 == :Optimal) && (obj1 > 0)
                    cut_num += 1
                    CGLP_cut = true
                    aff = AffExpr(x, cut_coef1, 0)      #stores the affine expression of the inequality
                    @constraint(m2, aff <= cut_rhs1)
                elseif status1 != :Optimal
                    status_CGLP == status1
                    info("CGLP terminated, unbounded or infeasible!")
                end

            end

        end

        cut_added = CGLP_cut || subgrad_cut

        #getting the initial solution of the model if a cut has been added
        if cut_added

            clock_end = time()
            time_elapsed = float(clock_end - clock_start)
            time_remain = max(0, time_limit - time_elapsed)     # compute the remaining time
            CplexSolver(CPX_PARAM_TILIM = time_remain)        #sets time limit parameter for the solver

            status = solve(m2)
            obj_val = getobjectivevalue(m2)
            x_val = getvalue(x)
            m_output = m2
            x_output = x_val
            obj_output = obj_val
            status_output = :Not_Optimal
        end

        iter_ecp += 1

    end

    return m_output, x_output, obj_output, status_output

end
