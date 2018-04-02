# Evaluation methods for additively separable functions

import Base.Cartesian.lreplace

"""Oracle that yields values of the parts of a separable function"""
mutable struct SeparableFunctionEvaluation
    values::Array{Dict}    # values of each term per variable-value pair
    max_vals::Array{Real}  # maximum values of each term per variable
    min_vals::Array{Real}  # minimum values of each term per variable
    constant::Real         # constant term of function
end

"""Negate expressions"""
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

separate_additively(expr::Real)::Array = [expr]

"""Return the variable indices involved in an expression"""
function get_vars_expr(expr::Expr)
    vars = Set()
    exprs_to_process = [expr]
    while !isempty(exprs_to_process)
        cur_expr = pop!(exprs_to_process)
        if cur_expr.head == :call
            for i in 2:length(cur_expr.args)
                if isa(cur_expr.args[i], Expr)
                    push!(exprs_to_process, cur_expr.args[i])
                end
            end
        elseif cur_expr.head == :ref
            if cur_expr.args[1] != :x
                error("Expression in constraint not supported (", expr, ")")
            end
            push!(vars, cur_expr.args[2]) # found variable, store index
        else
            error("Expression in constraint not supported (", expr, ")")
        end
    end
    return vars
end

"""Evaluate an expression replacing any variable by the given value"""
function eval_var(expr::Expr, val::Real)
    # Replace references by the value
    if expr.head == :ref
        return val
    end
    copy_expr = copy(expr)
    exprs_to_process = [copy_expr]
    while !isempty(exprs_to_process)
        cur_expr = pop!(exprs_to_process)
        if cur_expr.head == :call
            for i in 2:length(cur_expr.args)
                if isa(cur_expr.args[i], Expr)
                    if cur_expr.args[i].head == :ref
                        cur_expr.args[i] = val
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

"""Return the parts of an additively separable constraint. If the constraint is
not additively separable, throws error."""
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

function evaluate_separable_constraints(m::Model)::Array{SeparableFunctionEvaluation}
    d = JuMP.NLPEvaluator(m)
    MathProgBase.initialize(d, [:ExprGraph])

    # Gather parts for separable constraints
    constr_parts = []
    constr_constants = []
    nconstrs = MathProgBase.numconstr(m)
    nvars = MathProgBase.numvar(m)
    evals = Array{SeparableFunctionEvaluation}(nconstrs)
    for i in 1:nconstrs
        expr = MathProgBase.constr_expr(d, i)
        parts = separate_additively(expr)

        # Separate constant part from expressions
        constant_indices = (i for (i, p) in enumerate(parts) if isa(p, Real))
        constant = sum(parts[i] for i in constant_indices)
        deleteat!(parts, constant_indices)

        # Organize parts per variables
        var_parts = Array{Expr}(nvars)
        for p in parts
            v = first(get_vars_expr(p))
            if !isassigned(var_parts, v)
                var_parts[v] = p
            else
                var_parts[v] = Expr(:call, :+, var_parts[v], p)
            end
        end

        # Assume the domain consists of all integer values within the bounds of the variable
        constr_evals = Array{Dict}(nvars)
        min_evals = Array{Real}(nvars)
        max_evals = Array{Real}(nvars)
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

        evals[i] = SeparableFunctionEvaluation(constr_evals, min_evals, max_evals, constant)
    end

    return evals
end
