
using JuMP
using Clp

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




"""
Create a pricing model given an instance of PricingSpecs.

# Input
- `pricing_spec`: An instance of PricingSpecs.

# Output
- `pricing_model`: Pricing optimization model.
"""
function create_pricing_model(pricing_spec::PricingSpecs)::JuMP.Model

	#m = Model(solver = ClpSolver())
	m = Model()

	#@variable(m, x[i=1:pricing_spec.var_num], lowerbound = pricing_spec.lb[i], upperbound = pricing_spec.ub[i])
	@variable(m, x[i=1:pricing_spec.var_num], lowerbound = pricing_spec.lb[i], upperbound = pricing_spec.ub[i], Int)

	scale::Array{Int, 1} = pricing_spec.ub											# variable x is divided by scale everywhere in the model
	obj::JuMP.AffExpr = AffExpr(x, pricing_spec.obj_c, 0.0)				#the objective function is linear
	@objective(m, Max, obj)

	@constraintref c[1:5]
	@NLconstraint(m, c[i=1:pricing_spec.pricing_constr_num], sum(pricing_spec.constr_c[i,j]*(x[j]/scale[j])*e^(-(x[j]/scale[j])^pricing_spec.constr_d[i,j]) for j=1:pricing_spec.var_num) <= pricing_spec.constr_rhs[i])

	return m

end



"""
Create a the initial LP of the pricing model.

# Input
- `pricing_spec`: An instance of PricingSpecs.
"""
function create_pricing_LP(pricing_spec::PricingSpecs)::JuMP.Model

	#m = Model(solver = ClpSolver())
	m = Model(solver = CplexSolver())

	@variable(m, x[i=1:pricing_spec.var_num], lowerbound = pricing_spec.lb[i], upperbound = pricing_spec.ub[i])

	scale::Array{Int, 1} = pricing_spec.ub											# variable x is divided by scale everywhere in the model
	obj::JuMP.AffExpr = AffExpr(x, pricing_spec.obj_c, 0.0)				#the objective function is linear
	@objective(m, Max, obj)

	return m

end
