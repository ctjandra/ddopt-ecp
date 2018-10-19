include("../../dd/construction.jl")
include("../cut/cutgenerator.jl")

"""
Perform ECP using DD cuts.

# Input
- `m::JuMP.Model`: The current model to which the cuts will be added.
- `dd_list::Array{DecisionDiagram}`: An array of decision diagrams corresponding to inequalities of the current model.
- `constraint_list::Array{Function}`: An array of functions representing variable part of constraints with one-to-one correspondence to the DD array.
- `constant_list`: Array of constant part of constraints where all moved to the lhs of <= 0 form.
- `time_limit::Float64`: Time limit to run the algorithm (optional).
- `iter_max_ecp`: The maximum number of iterations to execute the algorithm (optional).
- `max_cut_num`: The maximum number of cuts (one cut for each violated inequality) to be added at EACH iteration (optional).
- `lp_to_ip_iter`: The number of iterations to solve as LP before switching to IP (optional).
- `lp_to_ip_tol`: The relative constraint violation threshold in the LP mode before swithcing to IP (optional).
- `obj_improvement_tol`: The improvement tolerance for the objective function (optional).
- `tol_num_ecp`: The number of past objective values (including the current one) wrt which the objective improvement tolerance is computed (optional).



# Output
- `x_output`: Optimal solution.
- `obj_output`: Optimal value.
- `status_output`: Status of the output model.
- `initial_lp_obj`: Objective value of the starting LP.
- `iter_ecp`: Total number of ECP iterations.
- `total_cuts`: Total number of added cuts.
- `ecp_time`: Total time elapsed.
- `relative_obj_improvement`: The relative improvement of the objective function value over the last few (tol_num_ecp) iterations.
"""

function ECP!(m::JuMP.Model, dd_list::Array{DecisionDiagram, 1}, constraint_list::Array{Function, 1}, constant_list::Array{Float64, 1};
    time_limit::Float64 = 10000.0, iter_max_ecp::Int = 1000, max_cut_num::Int = 10,
    lp_to_ip_iter::Int = 500, lp_to_ip_tol::Float64 = 0.00001, obj_improvement_tol::Float64 = 1e-8, tol_num_ecp::Int = 10)::
    Tuple{Array{Float64, 1}, Float64, Symbol, Float64, Int, Int, Float64, Float64}

    time_start::Float64 = time()

    @assert (length(dd_list) == length(constraint_list))
    @assert (length(constant_list) == length(constraint_list))

    n::Int = MathProgBase.numvar(m)
    # stopping criteria for the ECP algorithm
    obj_history_ecp::Array{Float64, 1} = [Inf for i=1:tol_num_ecp]    # the vector that stores the previous (improved) objective values for tolerance check
    relative_obj_improvement::Float64 = NaN    # The relative improvement of the objective function value over the last few (tol_num_ecp) iterations

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
                    #popfirst!(obj_history_ecp)         #not supported for Julia .6
                    deleteat!(obj_history_ecp, 1)
                    relative_obj_improvement = (obj_history_ecp[1] - obj_val)/(abs(obj_val) + 1e-8)
                    p = p && (relative_obj_improvement > obj_improvement_tol)
                end
            end
        end
        return p    # is false when the stopping criterion is met
    end

    #*******************************************************************************

    # initializing cut gnerating paratmeters
    subgrad_flag::Bool = true       # if subgradient method is used to generate cuts
    CGLP_flag::Bool = false             # if CGLP method is used to generate cuts
    use_CGLP::Bool = false          # determines whether CGLP will be used at all in the ECP
    st_rule::Array{Int, 1} = [0, 1, 2]        # sets the stopping rule for the algorithm
    cut_added::Bool = false       # indicates whether a cut (of any type) has been added to the model
    con_num::Int = length(constraint_list)          # the number of constraint
    subgrad_start_point::Array{Float64, 2} = zeros(con_num,n)  # stores the starting point for each constraint to be used in the subgradient method
    violation_val::Array{Float64, 1} = Array{Float64}(  con_num)     # stores the violation value for each constraint at the current optimal point
    normalized_violation_val::Array{Float64, 1} = Array{Float64}(con_num)     # stores the violation value after being normalized wrt the rhs
    iter_ecp::Int = 1                                # iteration counter
    total_cuts::Int = 0                                 # total number of cuts added

    # set parameters to switch between solving LP and IP at iterations
    lp_flag::Bool = true          # if the LP method is used to solve the current model
    ip_flag::Bool = false         # if the IP method is used to solve the current model
    use_ip::Bool = false          # if the model is supposed to swotch to IP at some point (must be false for continuous programs)

    # initializing time parameters
    clock_start::Float64 = time()
    CplexSolver(CPX_PARAM_TILIM = time_limit)        #sets time limit parameter for the solver
    #ClpSolver(MaximumSeconds = time_limit)

    # get the initial solution of the model
    x = getindex(m, :x)                 # redefines variables in model m
    status::Symbol = solve(m)
    obj_val::Float64 = getobjectivevalue(m)
    initial_lp_obj = obj_val
    x_val = Array{Float64, 1}(n)    # DD variables and constraint functions are based on the flattened out form in MathProgBase, seo we need to store variable values in that format
    for i=1:n
        x_val[i] = getvalue(Variable(m,i))
        #x_val[i] = getvalue(x[i])
    end
    m_output::JuMP.Model = m
    x_output::Array{Float64, 1} = x_val         #array copy, change of x_val changes x_output
    obj_output::Float64 = obj_val
    status_output::Symbol = :Not_Optimal

    #************************************************************
    #performing the ecp method until the stopping criteria is met

    #println("@@@@@@@@@@@@@@@@@@@@", stop_rule_ecp(st_rule))

    while stop_rule_ecp(st_rule)

        # if the previous iteration has not added any cuts to the model:
        if cut_added == false && iter_ecp > 1
            # 1. if the algorithm is allowed to switch to CGLP
            if use_CGLP == true
                # 1.1. if CGLP method has not been activated: activate it
                if CGLP_flag == false
                    CGLP_flag = true
                    subgrad_flag = false
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
            violation_val[i] = constraint_list[i](x_val) + constant_list[i]
            normalized_violation_val[i] = violation_val[i]/max(1e-9, abs(constant_list[i]))
        end

        sort_ind::Array{Int, 1} = sortperm(normalized_violation_val)     #stores the indices of the ascending sorted elements of violation_val

        # if the current point satisfies all constraints, it is global optimal
        if violation_val[sort_ind[end]] <= 0 + 1e-8
            println("Global optimal solution found!")
            status_output = :Global_Optimal
            break
        end

        # check whether to solve IP or LP based on the violation tolerance and the iteration number
        if (lp_flag == true) && (use_ip == true)
            if (normalized_violation_val[sort_ind[end]] <= lp_to_ip_tol) || (iter_ecp > lp_to_ip_iter)
                lp_flag = false
                ip_flag = true
                println("Switching to solve IP...")
                for i=1:n
                    setcategory(Variable(m,i), :Int)        # makes variable type with "linear index" i in its flattened out form, integer
                end
            end
        end

        cut_num = 0       # number of cuts added so far
        CGLP_cut = false    # CGLP cut has been added
        subgrad_cut = false     # subgradient cut has been added
        cut_added = false      # any cut has been added
        status_CGLP = :Optimal # it is assumed that CGLP will solve to optimality by default

        # analyze each constraint
        for i::Int=con_num:-1:1

            clock_end::Float64 = time()
            time_elapsed::Float64 = float(clock_end - clock_start)
            c_index::Int = sort_ind[i]

            # stop the loop if the inequality is satisfied or the number of allowable cuts is fulfilled or time exceeds the limit
            if (max_cut_num <= cut_num) || (violation_val[c_index] <= 0 + 1e-8) || (time_elapsed > time_limit)
                break
            end

            time_remain::Float64 = max(0.0, time_limit - time_elapsed)     # compute the remaining time

            if subgrad_flag == true       # uses subgradient method to generate cut
                obj1::Float64, cut_coef1::Array{Float64, 1}, cut_rhs1::Float64 = subgradient(dd_list[c_index], x_val, step_rule = 0, starting_point = subgrad_start_point[c_index,:], tilim = time_remain)

                if obj1 > 0 + 1e-7            # if the inequality cuts off the point
                    #println("\nsubgrad alg used\n")
                    cut_num += 1
                    subgrad_cut = true
                    #aff::JuMP.AffExpr = AffExpr(x, cut_coef1, 0)      # stores the affine expression of the inequality
                    aff = sum(Variable(m, i)*cut_coef1[i] for i=1:n)      # operates on variables in their original form (converted from flattened form)
                    @constraint(m, aff <= cut_rhs1)
                    subgrad_start_point[c_index,:] = cut_coef1      # uses the current cut's normal vector as the initial point for the next subgrad iteration
                end
            end

            if (CGLP_flag == true) && (use_CGLP == true)       # uses CGLP to generate cut

                println("SOLVING CGLP NOW")

                obj1, cut_coef1, cut_rhs1, status1 = CGLP(dd_list[c_index], x_val, tilim = time_remain)

                if (status1 == :Optimal) && (obj1 > 0 + 1e-8)
                    cut_num += 1
                    CGLP_cut = true
                    #aff = AffExpr(x, cut_coef1, 0)      #stores the affine expression of the inequality
                    aff = sum(Variable(m, i)*cut_coef1[i] for i=1:n)      # operates on variables in their original form (converted from flattened form)
                    @constraint(m, aff <= cut_rhs1)
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
            #ClpSolver(MaximumSeconds = time_remain)

            status = solve(m)
            println("\n\nLinear CUT:\t", MathProgBase.numlinconstr(m), "\n\n")
            obj_val = getobjectivevalue(m)
            println("Objective values: ", obj_val)
            if status == :Optimal               # the model has been solved to optimality
                for i=1:n
                    x_val[i] = getvalue(Variable(m,i))
                    #x_val[i] = getvalue(x[i])
                end
                x_output = x_val       # the change already takes effect because of array copy
                obj_output = obj_val
                status_output = :Not_Optimal
            end
        end

        total_cuts += cut_num
        iter_ecp += 1

    end

    time_end::Float64 = time()
    total_time_elapsed::Float64 = float(time_end - time_start)

    return x_output, obj_output, status_output, initial_lp_obj, iter_ecp, total_cuts, total_time_elapsed, relative_obj_improvement

end
