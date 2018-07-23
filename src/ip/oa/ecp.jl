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

function ECP(n::Int64, m::JuMP.Model, dd_list::Array{DecisionDiagram, 1}, con_list::Array{SeparableFunctionEvaluation, 1}; time_limit::Float64 = 10000.0)::Tuple{JuMP.Model, Array{Float64, 1}, Float64, Symbol}

    @assert (length(dd_list) == length(con_list))

    # stopping criteria for the ECP algorithm
    iter_max_ecp::Int = 100      # the maximum number of iterations to execute the algorithm
    tolerance_ecp::Float64 = 0.0001    # the tolerance for relative error
    tol_num_ecp::Int = 5         # the number of past objective values (including the current one) wrt which the tolerance is computed.
    obj_history_ecp::Array{Float64, 1} = [-Inf for i=1:tol_num_ecp]    # the vector that stores the previous (improved) objective values for tolerance check


    # function to define what stopping criterion must be used:
    # 0: time limit
    # 1: iteration number
    # 2: objective improvement: the relative difference between k consequtive improvements
    # 3: error on the optimality gap (TO DO: need a primal heuristic method)
    function stop_rule_ecp(stp::Array{Int64})::Bool   # multiple stopping criteria can be passed as an array
        p = true
        for v::Int in stp
            if v ==  0      # time limit criterion
                clock_end::Float64 = time()
                time_elapsed::Float64 = float(clock_end - clock_start)
                p = p && (time_elapsed <= time_limit)
                #println("v==0: ", p)
            end
            if v == 1       # the number of iteration criterion
                p = p && (iter_ecp <= iter_max_ecp)
                #println("v==1: ", p)
            end
            if v == 2       # uses the objective improvement criterion
                if cut_added   # if a cut has been added to the model, the `obj_history_ecp` vector must be updated and the criterion is checked
                    push!(obj_history_ecp, obj_val)
                    deleteat!(obj_history_ecp, 1)
                    p = p && ((obj_val - obj_history_ecp[1])/(abs(obj_val) + 1e-5) > tolerance_ecp)
                end
            end
        end
        return p    # is false when the stopping criterion is met
    end

    #*******************************************************************************
    m2::JuMP.Model = copy(m)
    x = getindex(m2,:x)     # allows to use the same variables of the orignal model for the new model


    # initializing cut gnerating paratmeters
    subgrad_flag::Bool = true       # if subgradient method is used to generate cuts
    CGLP_flag::Bool = false             # if CGLP method is used to generate cuts
    use_CGLP::Bool = true          # determines whether CGLP will be used at all in the ECP
    max_cut_num::Int = 10     # maximum number of cuts (one cut for each violated inequality) to be added at each iteration
    st_rule::Array{Int, 1} = [0, 1]        # sets the stopping rule for the algorithm
    cut_added::Bool = false       # indicates whether a cut (of any type) has been added to the model
    con_num::Int = length(con_list)          # the number of constraint
    subgrad_start_point::Array{Float64, 2} = zeros(con_num,n)  # stores the starting point for each constraint to be used in the subgradient method
    violation_val::Array{Float64, 1} = Array{Float64}(con_num)     # stores the violation value for each constraint at the current optimal point
    normalized_violation_val::Array{Float64, 1} = Array{Float64}(con_num)     # stores the violation value after being normalized wrt the rhs
    iter_ecp::Int = 1                                # iteration counter

    # set parameters to switch between solving LP and IP at iterations
    lp_flag::Bool = true          # if the LP method is used to solve the current model
    ip_flag::Bool = false         # if the IP method is used to solve the current model
    lp_to_ip_iter::Int = 100       # number of iterations to solve as LP before switching to IP
    lp_to_ip_tol::Float64 = 0.0001     # relative constraint violation threshold in the LP mode before swithcing to IP

    # initializing time parameters
    clock_start::Float64 = time()
    CplexSolver(CPX_PARAM_TILIM = time_limit)        #sets time limit parameter for the solver

    # get the initial solution of the model
    status::Symbol = solve(m2)
    obj_val::Float64 = getobjectivevalue(m2)
    x_val::Array{Float64, 1} = getvalue(x)
    m_output::JuMP.Model = m2
    x_output::Array{Float64, 1} = x_val
    obj_output::Float64 = obj_val
    status_output::Symbol = :Not_Optimal

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
        for i::Int=1:con_num
            violation_val[i] = con_list[i].single_function(x_val) + con_list[i].constant
            normalized_violation_val[i] = violation_val[i]/max(1e-9, abs(con_list[i].constant))
        end

        sort_ind::Array{Int, 1} = sortperm(normalized_violation_val)     #stores the indices of the ascending sorted elements of violation_val

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
                println("Switching to solve IP...")
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

        for i::Int=con_num:-1:1

            clock_end::Float64 = time()
            time_elapsed::Float64 = float(clock_end - clock_start)
            c_index::Int = sort_ind[i]

            # stop the loop if the inequality is satisfied or the number of allowable cuts is fulfilled or time exceeds the limit
            if (max_cut_num < cut_num) || (violation_val[c_index] <= 0) || (time_elapsed > time_limit)
                break
            end

            time_remain::Float64 = max(0.0, time_limit - time_elapsed)     # compute the remaining time

            if subgrad_flag == true       # uses subgradient method to generate cut
                obj1::Float64, cut_coef1::Array{Float64, 1}, cut_rhs1::Float64 = @time subgradient(dd_list[c_index], x_val, step_rule = 0, starting_point = subgrad_start_point[c_index,:], tilim = time_remain)

                println("\nsubgrad alg used\n")

                if obj1 > 0             # if the inequality cuts off the point
                    cut_num += 1
                    subgrad_cut = true
                    aff::JuMP.AffExpr = AffExpr(x, cut_coef1, 0)      # stores the affine expression of the inequality
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
            time_remain = max(0.0, time_limit - time_elapsed)     # compute the remaining time
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
