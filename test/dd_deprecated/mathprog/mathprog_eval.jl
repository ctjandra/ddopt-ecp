# Evaluation methods for additively separable functions

import Base.Cartesian.lreplace


"""Given an array of expressions, returns the array of negation of expressions"""
function negate_exprs(exprs::Array)::Array
    negated_exprs = []
    for expr in exprs
        if isa(expr, Real)
            push!(negated_exprs, -expr)
        else
            push!(negated_exprs, Expr(:call, :-, expr))
        end
    end
    return negated_exprs
end


"""Return the variable indices involved in an expression"""
function get_vars_expr(expr::Expr)
    vars = Set()
    exprs_to_process = [expr]               #stores remaining (sub)expressions in expr that need to be preocessed
    while !isempty(exprs_to_process)
        cur_expr = pop!(exprs_to_process)
        if cur_expr.head == :call           #if the current expression is a "call" type, e.g., x + y or 2*x, ...
            for i in 2:length(cur_expr.args)
                if isa(cur_expr.args[i], Expr)  #if an argument is exprresion, it is added to the list to be processed
                    push!(exprs_to_process, cur_expr.args[i])
                end
            end
        elseif cur_expr.head == :ref        #if the current expression is a "reference" type, e.g., x[i]
            if cur_expr.args[1] != :x
                error("Variables in expression not supported (", expr, ")")
            end
            push!(vars, cur_expr.args[2]) # found variable, store index
        else
            error("Expression in constraint not supported (", expr, ")")
        end
    end
    return vars
end



"""Decomposes a given expression into separable terms.

Input:
expr: the expression.

Output:
An array that stores decomposed parts (expressions) of the constraint.
These parts are separated by + sing, e.g., 3x[1] - x[1]^2 + sin(x[2]) - 4 <= 0, will be stored in an array with 4 elements for the lhs.
If the constraint is not additively separable, throws error.
"""

separate_additively(expr::Real)::Array = [expr]

function separate_additively(expr::Expr)::Array
    # Examine the expression tree: first split (in)equalities and then check if
    # each term separated by a + or - uses only one type of variable
    println(expr.head, " -- ", expr.args, " -- ", expr.typ)
    if expr.head == :call
        # println("Call")
        # (In)equality or non-unary -
        if expr.args[1] in [:(<=), :(==), :(>=)] || (expr.args[1] == :- && length(expr.args) >= 3)
            # If (in)equality, move right-hand side to left-hand side by negating it
            # If subtraction, negate second expression
            exprs = []
            append!(exprs, separate_additively(expr.args[2]))
            append!(exprs, negate_exprs(separate_additively(expr.args[3])))
            return exprs
        # Non-unary +
        elseif expr.args[1] in [:+] && length(expr.args) >= 3
            exprs = []
            for i in 2:length(expr.args)
                append!(exprs, separate_additively(expr.args[i]))
            end
            return exprs
        # Unary +
        elseif expr.args[1] in [:+] # unary; binary operator already checked
            return separate_additively(expr.args[2])
        # Unary -
        elseif expr.args[1] in [:-] # unary; binary operator already checked
            return negate_exprs(separate_additively(expr.args[2]))
        # Process expression as a single component; must only contain one type of variable
        else
            vars = get_vars_expr(expr)
            if length(vars) >= 2
                error("Non-additively separable constraint in model")
            end
            return [expr]
        end
    elseif expr.head == :comparison
        error("Only one-sided constraints are supported")  # to make DD construction easier
    elseif expr.head == :ref
        return [expr]  # isolated variable
    else
        error("Expression in constraint not supported (", expr, ")")
    end
end


"""
Function: Returns the sense of a given expression in form of an (in)equality.
Returns error if the expression is not (in)equality
"""

function get_constraint_sense(expr::Expr)
    if expr.head == :call && expr.args[1] in [:(<=), :(==), :(>=)]
        return a.args[1]            #Sybmbol for the (in)equality sense
    else
        error("Expression not in constraint form (", expr, ")")
    end
end


#************************************************************************
#*************************************************************************
#*************************************************************************
"""Evaluate a univariate expression replacing the variable by the given value"""
function eval_var(expr::Expr, val::Real)
    # Replace references by the value
    if expr.head == :ref        #if the expression is composed of a single variable assigns the value
        return val
    else                        #otherwise decomposes the expression down to variables
        copy_expr = copy(expr)
        exprs_to_process = [copy_expr]
        while !isempty(exprs_to_process)
            cur_expr = pop!(exprs_to_process)
            if cur_expr.head == :call
                for i in 2:length(cur_expr.args)
                    if isa(cur_expr.args[i], Expr)
                        if cur_expr.args[i].head == :ref
                            cur_expr.args[i] = val          #replacing the variable with its value
                        else
                            push!(exprs_to_process, cur_expr.args[i])
                        end
                    end
                end
            else
                error("Expression in constraint not supported (", expr, ")")
            end
        end
        return eval(copy_expr)
    end
end

"""
TODO: Another method to evaluate a given expression, where we allow several variables (for nonseparable case!)
"""


"""
Function: calculates min and max value of a univariate expression over a given domain

Input:
expr: Expression
lb: lower bound
ub: upper bound

Output:
min_val
max_val
"""
function min_max_calculator(expr::Expr, doms::Array{Range})
    min_val = Inf
    max_val = -Inf
    for j in doms
        val = eval_var(expr, j)
        min_val = min(val, min_val)
        max_val = max(val, max_val)
    end
    return min_val, max_val
end

#An alternative method, where it takes a univarate function as input and computes the min and max over given domain
function min_max_calculator(func::Function, doms::Array{Range})
    min_val = Inf
    max_val = -Inf
    for j in doms
        val = func(j)
        min_val = min(val, min_val)
        max_val = max(val, max_val)
    end
    return min_val, max_val
end




"""Oracle that yields decomposition terms of a given function"""
mutable struct SeparableFunctionEvaluation
    separate_exprs::Array{Expr}     # each element stores the univariate expression for a variable
    separate_functions::Array{Function} # each element stores the univariate function for a variable
    constant::Real         # constant term of function
#    values::Array{Dict}    # values of each term per variable-value pair, each element corresponding to a variable (depends on the variables domain)
#    max_vals::Array{Real}  # maximum values of each term per variable (depends on the variables domain)
#    min_vals::Array{Real}  # minimum values of each term per variable (depends on the variables domain)
end


"""
Function: Returns decomposition terms of constraints of a model

Input:
m: the model

Output:
evals: An array of function specifications for constraints
"""

function evaluate_separable_constraints(m::Model)::Array{SeparableFunctionEvaluation}
    d = JuMP.NLPEvaluator(m)        #obtains an NLP evaluator object from a JuMP model, which enables query of function characteristics in the model
    MathProgBase.initialize(d, [:ExprGraph])    #enables query of the functions of the model in their expression graph form

    # Gather parts for separable constraints
    nconstrs = MathProgBase.numconstr(m)
    nvars = MathProgBase.numvar(m)
    evals = Array{SeparableFunctionEvaluation}(undef, nconstrs)

    for i in 1:nconstrs
        expr = MathProgBase.constr_expr(d, i)       #gets the expression of the ith constraint
        parts = separate_additively(expr)         #gets the array of decomposed part of the constraint

        # Put all constraints in f(x) <= 0 form by negating the lhs when >=.
        con_sense = get_constraint_sense(expr)
        if con_sense == :>=
            parts = negate_exprs(parts)
        # elseif con_sense == :==           #need to implement for == case: create two constraints with <=...
        end

        # Separate constant part from expressions
        constant_indices = (j for (j, p) in enumerate(parts) if isa(p, Real))
        if isempty(constant_indices)
            constant = 0
        else
            constant = sum(parts[j] for j in constant_indices)
            deleteat!(parts, constant_indices)
        end

        # Organize parts per variables
        var_parts = Array{Expr}(undef, nvars)
        for p in parts
            v = first(get_vars_expr(p))
            if !isassigned(var_parts, v)
                var_parts[v] = p
            else
                var_parts[v] = Expr(:call, :+, var_parts[v], p)
            end
        end

        # Create function corresponding to each univariate expression
        function_parts = Array{Function}(undef, nvars)
        for j in 1:nvars
            if isassigned(var_parts, j)
                function_parts[j] = (xx::Int64) -> eval_var(var_parts[j], xx)
            end
        end

        #Computes the attributes corresponding to function value over variable domain
        #Assumes the domain consists of all integer values within the bounds of the variable
        """
        constr_evals = Array{Dict}(undef, nvars)
        min_evals = Array{Real}(undef, nvars)
        max_evals = Array{Real}(undef, nvars)
        for v in 1:nvars
            var = Variable(m, v)
            lb = ceil(Int, getlowerbound(var))
            ub = floor(Int, getupperbound(var))
            if !isassigned(var_parts, v)
                constr_evals[v] = Dict(j => 0 for j in lb:ub)
                min_evals[v] = 0
                max_evals[v] = 0
            else
                constr_evals[v] = Dict()
                min_eval = Inf
                max_eval = -Inf
                for j in lb:ub
                    val = eval_var(var_parts[v], j)
                    min_eval = min(val, min_eval)
                    max_eval = max(val, max_eval)
                    constr_evals[v][j] = val
                end
                min_evals[v] = min_eval
                max_evals[v] = max_eval
            end
        end
        """
        evals[i] = SeparableFunctionEvaluation(var_parts, function_parts, constant, constr_evals, min_evals, max_evals)

    end

    return evals
end
