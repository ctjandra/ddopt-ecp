include("pricing_specs.jl")

using JuMP
using Clp

"""
Read a text file containing pricing problem data. Assumes file is correctly formatted.

# Input
- `filename`: Address of the file.

# Output
- `price_spec`: An instance of the PricingSpecs type.
"""
function read_pricing(filename::AbstractString)::PricingSpecs
	f = open(filename)
	for i=1:15			#this is the description part of the text file, so we skip all these lines
		readline(f)
	end

	line = readline(f)
	var_num = parse(Int, line)		#number of variables
	readline(f)						#an empty line

	readline(f)					#total number of constraints (no field to assign to, so we skip)
	readline(f)

	line = readline(f)
	pricing_constr_num = parse(Int, line)		#number of pricing constraints
	readline(f)

	line = readline(f)
	pieces = split(line)					#split strings into pieces separated by space (default delimiter)
	@assert(length(pieces) == var_num)
	lb = Array{Int}(var_num)
	for i=1:var_num
		lb[i] = parse(Int, pieces[i])		#lower bound vector
	end
	readline(f)

	line = readline(f)
	pieces = split(line)
	@assert(length(pieces) == var_num)
	ub = Array{Int}(var_num)
	for i=1:var_num
		ub[i] = parse(Int, pieces[i])		#upper bound vector
		#ub[i] = 2
	end
	readline(f)

	line = readline(f)
	pieces = split(line)
	@assert(length(pieces) == var_num)
	obj_c = Array{Float64}(var_num)
	for i=1:var_num
		obj_c[i] = parse(Float64, pieces[i])		#objective function coefficients
	end
	readline(f)

	line = readline(f)
	pieces = split(line)
	@assert(length(pieces) == var_num)
	obj_d = Array{Float64}(var_num)
	for i=1:var_num
		obj_d[i] = parse(Float64, pieces[i])		#objective function degree of variables
	end
	readline(f)

	constr_c = Array{Float64}(pricing_constr_num, var_num)
	for i=1:pricing_constr_num
		line = readline(f)
		pieces = split(line)
		@assert(length(pieces) == var_num)
		for j=1:var_num
			constr_c[i,j] = parse(Float64, pieces[j])		#coefficient of variable j in constraint i
		end
	end
	readline(f)

	constr_d = Array{Float64}(pricing_constr_num, var_num)
	for i=1:pricing_constr_num
		line = readline(f)
		pieces = split(line)
		@assert(length(pieces) == var_num)
		for j=1:var_num
			constr_d[i,j] = parse(Float64, pieces[j])		#degree of the exp power for variable j in constraint i
		end
	end
	readline(f)

	line = readline(f)
	pieces = split(line)
	@assert(length(pieces) == pricing_constr_num)
	constr_rhs = Array{Float64}(pricing_constr_num)
	for i=1:pricing_constr_num
		constr_rhs[i] = parse(Float64, pieces[i])		#rhs value of constraints
	end
	readline(f)

	price_spec = PricingSpecs(var_num, pricing_constr_num, lb, ub, obj_c, obj_d, constr_c, constr_d, constr_rhs)	#create an instance of PricingSpecs from the read data
	return price_spec

end



"""
Create a pricing model given an instance of PricingSpecs.

# Input
- `pricing_spec`: An instance of PricingSpecs.

# Output
- `pricing_model`: Pricing optimization model.
"""
function create_pricing_model(pricing_spec::PricingSpecs)::JuMP.Model

	m = Model(solver = ClpSolver())

	@variable(m, x[i=1:pricing_spec.var_num], lowerbound = pricing_spec.lb[i], upperbound = pricing_spec.ub[i])

	scale = pricing_spec.ub											# variable x is divided by scale everywhere in the model
	obj = AffExpr(x, pricing_spec.obj_c./scale, 0.0)				#the objective function is linear
	@objective(m, Min, obj)

	@NLconstraint(m, pricing_constr[i=1:pricing_spec.pricing_constr_num], sum(pricing_spec.constr_c[i,j]*(x[j]/scale[j])*e^(-(x[j]/scale[j])^pricing_spec.constr_d[i,j]) for j=1:pricing_spec.var_num) <= pricing_spec.constr_rhs[i])

	return m

end
