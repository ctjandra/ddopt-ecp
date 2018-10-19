# Evaluation methods for additively separable functions

import Base.Cartesian.lreplace

using JuMP
using Memoize


"""
Return the array of negation of expressions.

# Input
- `exprs::Array`: An array of expressions.
"""
function negate_exprs(exprs::Array{Expr, 1})::Array{Expr, 1}
    negated_exprs = Array{Expr, 1}()
    for expr in exprs
        push!(negated_exprs, Expr(:call, :-, expr))
    end
    return negated_exprs
end


"""
Return the variable indices involved in an expression.
"""
# NOTE: The extracted expressions are already flattened out by MathProgBase, i.e., multidimentional variables are represented as 1 dimensional
# NOTE: Also, all variable names are renamed as x automatically
function get_vars_expr(expr::Expr)::IntSet
    vars = IntSet()
    exprs_to_process::Array{Expr, 1} = [expr]               #stores remaining (sub)expressions in expr that need to be preocessed
    while !isempty(exprs_to_process)
        cur_expr::Expr = pop!(exprs_to_process)
        if cur_expr.head == :call           #if the current expression is a "call" type, e.g., x + y or 2*x, ...
            for i::Int in 2:length(cur_expr.args)
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



"""
Decompose a given expression into separable terms.

# Output:
- An array that stores decomposed parts (expressions) of a given expression.
- The constant value of the expression
These parts are separated by ``+`` sign, e.g., ``3x[1] - x[1]^2 + sin(x[2]) - 4 <= 0``, will be stored in an array with 3 elements for the lhs and a constant -4
If the constraint is not additively separable, throws error.
"""

function separate_additively(expr::Expr)::Tuple{Array{Expr, 1}, Float64}
    # Examine the expression tree: first split (in)equalities and then check if
    # each term separated by a + or - uses only one type of variable
    # println(expr.head, " -- ", expr.args, " -- ", expr.typ)
    constant_val::Float64 = 0
    if expr.head == :call
        # println("Call")
        # (In)equality or non-unary -
        if expr.args[1] in [:(<=), :(==), :(>=)] || (expr.args[1] == :- && length(expr.args) >= 3)
            # If (in)equality, then move right-hand side to left-hand side by negating it
            # If subtraction, then negate second expression
            exprs1 = Array{Expr, 1}()
            if isa(expr.args[2], Float64)
                constant_val += expr.args[2]
            else
                a, b = separate_additively(expr.args[2])
                append!(exprs1, a)
                constant_val += b
            end
            if isa(expr.args[3], Float64)
                constant_val -= expr.args[3]
            else
                a, b = separate_additively(expr.args[3])
                append!(exprs1, negate_exprs(a))
                constant_val -= b
            end
            return exprs1, constant_val
        # Non-unary +
        elseif expr.args[1] in [:+] && length(expr.args) >= 3
            exprs = Array{Expr, 1}()
            for i in 2:length(expr.args)
                if isa(expr.args[i], Float64)
                    constant_val += expr.args[3]
                else
                    a, b = separate_additively(expr.args[i])
                    append!(exprs, a)
                    constant_val += b
                end
            end
            return exprs, constant_val
        # Unary +
        elseif expr.args[1] in [:+] # unary; binary operator already checked
            a, b = separate_additively(expr.args[2])
            return a, b
        # Unary -
        elseif expr.args[1] in [:-] # unary; binary operator already checked
            a, b = separate_additively(expr.args[2])
            return negate_exprs(a), -b
        # Process expression as a single component; must only contain one type of variable
        else
            vars::IntSet = get_vars_expr(expr)
            if length(vars) >= 2
                error("Non-additively separable constraint in model")
            end
            return [expr], 0
        end
    elseif expr.head == :comparison
        error("Only one-sided constraints are supported")  # to make DD construction easier
    elseif expr.head == :ref
        return [expr], 0  # isolated variable
    else
        error("Expression in constraint not supported (", expr, ")")
    end
end


"""
Return the sense of a given expression in form of an (in)equality.
Returns error if the expression is not (in)equality
"""

function get_constraint_sense(expr::Expr)::Symbol
    if expr.head == :call && expr.args[1] in [:(<=), :(==), :(>=)]
        return expr.args[1]            #Symbol for the (in)equality sense
    else
        error("Expression not in constraint form (", expr, ")")
    end
end


#************************************************************************
#*************************************************************************
#*************************************************************************

# Evaluate a univariate expression replacing the variable by the given value.
@memoize function eval_var(expr::Expr, val::Float64)::Float64
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

# TODO: Another method to evaluate a given expression, where we allow several variables (for nonseparable case!)



"""
Oracle that yields functional specifications of a constraint.
"""
struct SeparableFunctionEvaluation
#    separate_exprs::Array{Expr}        # each element stores the univariate expression for a variable
    single_function::Function           # evaluate the whole multivariate function in the constraint at any point (without the constant part)
    separate_functions::Array{Function, 1} # each element stores the univariate function for a variable
    constant::Float64                   # constant term of function
end


"""
Return decomposition terms of constraints of a model.

# Output
- `evals`: An array of instances that stores function specification for constraints
"""
function evaluate_separable_constraints(m::Model)::Array{SeparableFunctionEvaluation}
    d = JuMP.NLPEvaluator(m)        #obtains an NLP evaluator object from a JuMP model, which enables query of function characteristics in the model
    MathProgBase.initialize(d, [:ExprGraph])    #enables query of the functions of the model in their expression graph form

    # Gather parts for separable constraints
    nlinearcon = MathProgBase.numlinconstr(m)           # number of linear constraints (first ordered in NLPEvaluator)
    nquadcon = MathProgBase.numquadconstr(m)            # number of quadratic constraints (second ordered in NLPEvaluator)
    nlpcon = JuMP.numnlconstr(m)                        # number of nonlinear constraints, besides quadratic (third ordered in NLPEvaluator)
    nconstrs = nquadcon + nlpcon                    # total number of nonlinear constraints

    nvars = MathProgBase.numvar(m)
    evals = Array{SeparableFunctionEvaluation}(  nconstrs)

    for i in 1:nconstrs
        # NOTE: The extracted expressions are already flattened out by MathProgBase, i.e., multidimentional variables are represented as 1 dimensional
        # NOTE: Also, all variable names are renamed as x automatically
        expr = MathProgBase.constr_expr(d, i + nlinearcon)       #gets the expression of the ith constraint among the nonlinear ones
        parts, constant = separate_additively(expr)         #gets the array of decomposed part of the constraint

        # Put all constraints in f(x) <= 0 form by negating the lhs when >=.
        con_sense = get_constraint_sense(expr)
        if con_sense == :>=
            parts = negate_exprs(parts)
            constant = -constant
        # elseif con_sense == :==           #need to implement for == case: create two constraints with <=...
        end

        # Organize parts per variables
        var_parts = Array{Expr}(  nvars)
        for p in parts
            v = first(get_vars_expr(p))
            if !isassigned(var_parts, v)
                var_parts[v] = p
            else
                var_parts[v] = Expr(:call, :+, var_parts[v], p)
            end
        end

        # Create function corresponding to each univariate expression
        function_parts = Array{Function, 1}(  nvars)
        for j in 1:nvars
            if isassigned(var_parts, j)
                function_parts[j] = (x::Real) -> eval_var(var_parts[j], x)
            end
        end

        # Create a single multivariate function representing lhs of a constraint
        single_function = function (x::Array{T} where T<:Real)
            @assert(length(x) == nvars)
            sum_v = sum(function_parts[j](x[j]) for j in 1:nvars if isassigned(var_parts, j))
            return sum_v
        end

        evals[i] = SeparableFunctionEvaluation(single_function, function_parts, constant)

    end

    return evals
end
