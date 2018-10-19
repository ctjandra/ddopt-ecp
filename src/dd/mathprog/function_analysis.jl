
using Memoize

"""
Store functional properties of univariate functions over a domain interval
"""
struct UVFuncProp
    convex::Bool
    concave::Bool
    monotone::Bool
    lb::Float64             #lower bound on the variable range (domain) over which the properties hold
    ub::Float64             #upper bound on the variable range (domain) over which the properties hold
end


"""
Type to store a univariate function together with its properties
"""
mutable struct UVFunc
    f::Function
    prop::Array{UVFuncProp, 1}         # properties are stored in an array that represent intervals of domain over which a property holds
end


#=
"""
Compute convexity of a univariate function over a given domain

# Input:
- `f`: univariate function
- `dom`: a tuple of (lowerbound, upperbound, variable type)
"""
function compute_convexity(f::Function, dom::Tuple{Float64, Float64, Symbol})::Bool
end


"""
Compute monotonicity of a univariate function over a given domain

# Input:
- `f`: univariate function
- `dom`: a tuple of (lowerbound, upperbound, variable type)
"""
function compute_convexity(f::Function, dom::Tuple{Float64, Float64, Symbol})::Bool
end
=#

#=
"""
Identify the functional properties of the univariate functions for each variable in the constraint.

# Input:
- `func`
"""
function identify_property!(func::UVFunc)
    # each identified property over an interval is added at the field `prop` of func
end
=#


"""
Update the functional properties of the univariate function using the information given as input.

# Input:
- `func`
- `p`: Array of interval properties
"""
function update_property!(func::UVFunc, p::Array{UVFuncProp, 1})
    func.prop = copy(p)
end



"""
Compute a lower bound (minimum) of a univariate function over a domain

# Input:
- `f`: univariate function
- `prop`: array of intervals in the function domain, where each interval has functional properties
- `dom`: a tuple of (lowerbound, upperbound, variable type)
"""
function uvf_lb(f::Function, prop::Array{UVFuncProp, 1}, dom::Tuple{Float64, Float64, Symbol})::Float64
    # prop stores consequitive intervals on variable domain together with function properties over that interval
    lb = -Inf
    # find the intervals that cover the input domain
    for i=1:length(prop)
        if (dom[1] >= prop[i].lb) && (dom[1] <= prop[i].ub)     # the position of lb of domain is found,
            lb = Inf
            j = i
            while (dom[2] >= prop[j].lb)
                lb1 = max(prop[j].lb, dom[1])           # finds the startpoint of interval
                ub1 = min(prop[j].ub, dom[2])           # finds the endpoint of interval
                # if the function is monotone or concave, evaluate the lower and upper bound only
                if prop[j].monotone || prop[j].concave
                    if (dom[3] == :Int)
                        lb0 = min(f(ceil(lb1)), f(floor(ub1)))
                    else
                        lb0 = min(f(lb1), f(ub1))
                    end
                # if the variable is integer and the range is not too large, an enumeration yeilds an exact minimum over domain
                elseif (dom[3] == :Int) && (ub1 - lb1 <= 1000)
                    lb0 = Inf
                    for j in ceil(Int,lb1):floor(Int,ub1)
                        lb0 = min(f(float(j)), lb0)
                    end
                #TODO: if the function is convex a lowerbound is computed using univariate minimization packages and then using convexity inequality
                #TODO: otherwise, a univariate minimization should be run using optim.jl, howerver there is no guarantee for lower bound!
                #TODO: one way is to decompose the univariate function into separable underestimator (if it is a product), or serveral simple functions with above properties, and take lb of each
                else
                    lb0 = -Inf
                end
                lb = min(lb, lb0)
                j = j+1
                if j > length(prop)
                    break
                end
            end
            break
        end
    end
    return lb
end



"""
Compute an upper bound (maximum) of a univariate function over a domain

# Input:
- `f`: univariate function
- `prop`: array of intervals in the function domain, where each interval has functional properties
- `dom`: a tuple of (lowerbound, upperbound, variable type)
"""
function uvf_ub(f::Function, prop::Array{UVFuncProp, 1}, dom::Tuple{Float64, Float64, Symbol})::Float64
    # prop stores consequitive intervals on variable domain together with function properties over that interval
    ub = Inf
    # find the intervals that cover the input domain
    for i=1:length(prop)
        if (dom[1] >= prop[i].lb) && (dom[1] <= prop[i].ub)     # the position of lb of domain is found,
            ub = -Inf
            j = i
            while (dom[2] >= prop[j].lb)
                lb1 = max(prop[j].lb, dom[1])           # finds the startpoint of interval
                ub1 = min(prop[j].ub, dom[2])           # finds the endpoint of interval
                # if the function is monotone or concave, evaluate the lower and upper bound only
                if prop[j].monotone || prop[j].convex
                    if (dom[3] == :Int)
                        ub0 = max(f(ceil(lb1)), f(floor(ub1)))
                    else
                        ub0 = max(f(lb1), f(ub1))
                    end
                # if the variable is integer and the range is not too large, an enumeration yeilds an exact minimum over domain
            elseif (dom[3] == :Int) && (ub1 - lb1 <= 1000)
                    ub0 = -Inf
                    for j in ceil(Int,lb1):floor(Int,ub1)
                        lb0 = max(f(float(j)), ub0)
                    end
                #TODO: if the function is convex a lowerbound is computed using univariate minimization packages and then using convexity inequality
                #TODO: otherwise, a univariate minimization should be run using optim.jl, howerver there is no guarantee for lower bound!
                #TODO: one way is to decompose the univariate function into separable underestimator (if it is a product), or serveral simple functions with above properties, and take lb of each
                else
                    ub0 = Inf
                end
                ub = max(ub, ub0)
                j = j+1
                if j > length(prop)
                    break
                end
            end
            break
        end
    end
    return ub
end
