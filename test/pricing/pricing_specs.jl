
"""
Pricing problem structure.

# Description
The separable functions in the constraints is of the form ``f(x_i) = a_i exp(-x_i^(k_i))``, and the objective is linear.
"""
mutable struct PricingSpecs
	var_num::Int64				#number of variables
    pricing_constr_num::Int64	#number of pricing constraints
	lb::Array{Int, 1}				#vector of lower bounds
	ub::Array{Int, 1}				#vector of upper bounds
	obj_c::Array{Float64, 1}		#vector of objective coeffcients
	obj_d::Array{Float64, 1}		#vector of degrees of variables in the objective
	constr_c::Matrix{Float64}	#coefficient matrix of variables in constraints
	constr_d::Matrix{Float64}	#degree matrix of variables in constraints
	constr_rhs::Array{Float64, 1}	#rhs value of constraints
end
