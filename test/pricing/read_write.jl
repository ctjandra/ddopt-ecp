include("pricing_specs.jl")


"""
Read a text file containing pricing problem data. Assumes file is correctly formatted.

# Input
- `filename`: Address of the file.

# Output
- `price_spec`: An instance of the PricingSpecs type.
"""
function read_pricing(filename::AbstractString)::PricingSpecs
	f = open(filename)
	for i::Int=1:15			#this is the description part of the text file, so we skip all these lines
		readline(f)
	end

	line = readline(f)
	var_num::Int = parse(Int, line)		#number of variables
	readline(f)						#an empty line

	readline(f)					#total number of constraints (no field to assign to, so we skip)
	readline(f)

	line = readline(f)
	pricing_constr_num::Int = parse(Int, line)		#number of pricing constraints
	readline(f)

	line = readline(f)
	pieces = split(line)					#split strings into pieces separated by space (default delimiter)
	@assert(length(pieces) == var_num)
	lb::Array{Int, 1} = Array{Int, 1}(  var_num)
	for i::Int=1:var_num
		lb[i] = parse(Int, pieces[i])		#lower bound vector
	end
	readline(f)

	line = readline(f)
	pieces = split(line)
	@assert(length(pieces) == var_num)
	ub::Array{Int, 1} = Array{Int}(  var_num)
	for i::Int=1:var_num
		ub[i] = parse(Int, pieces[i])		#upper bound vector
		#ub[i] = 2
	end
	readline(f)

	line = readline(f)
	pieces = split(line)
	@assert(length(pieces) == var_num)
	obj_c::Array{Float64, 1} = Array{Float64, 1}(  var_num)
	for i::Int=1:var_num
		obj_c[i] = parse(Float64, pieces[i])		#objective function coefficients
	end
	readline(f)

	line = readline(f)
	pieces = split(line)
	@assert(length(pieces) == var_num)
	obj_d::Array{Float64, 1} = Array{Float64, 1}(  var_num)
	for i::Int=1:var_num
		obj_d[i] = parse(Float64, pieces[i])		#objective function degree of variables
	end
	readline(f)

	constr_c::Array{Float64, 2} = Array{Float64, 2}(  pricing_constr_num, var_num)
	for i::Int=1:pricing_constr_num
		line = readline(f)
		pieces = split(line)
		@assert(length(pieces) == var_num)
		for j::Int=1:var_num
			constr_c[i,j] = parse(Float64, pieces[j])		#coefficient of variable j in constraint i
		end
	end
	readline(f)

	constr_d::Array{Float64, 2} = Array{Float64, 2}(  pricing_constr_num, var_num)
	for i::Int=1:pricing_constr_num
		line = readline(f)
		pieces = split(line)
		@assert(length(pieces) == var_num)
		for j::Int=1:var_num
			constr_d[i,j] = parse(Float64, pieces[j])		#degree of the exp power for variable j in constraint i
		end
	end
	readline(f)

	line = readline(f)
	pieces = split(line)
	@assert(length(pieces) == pricing_constr_num)
	constr_rhs::Array{Float64, 1} = Array{Float64, 1}(  pricing_constr_num)
	for i::Int=1:pricing_constr_num
		constr_rhs[i] = parse(Float64, pieces[i])		#rhs value of constraints
	end
	readline(f)

	price_spec = PricingSpecs(var_num, pricing_constr_num, lb, ub, obj_c, obj_d, constr_c, constr_d, constr_rhs)	#create an instance of PricingSpecs from the read data
	return price_spec

end





"""
Write an output text file

# Input
- `filename`: Address of the file.
- `time_limit`: Time limit.
- `iter_max_ecp`: Max limit on ecp iterations
- `max_cut_num`: Max cut number added at EACH iteration of ecp
- `lp_to_ip_iter`: Number of ecp iterations to solve LP before switching to IP (if applicable)
- `lp_to_ip_tol`: Constraint violation tolerance below which ecp switches to solving IP (if applicable)
- `obj_improvement_tol`: Objective improvement tolerance over the last few iterations
- `tol_num_ecp`: The number of last iterations over which the relative obj improvement is computed
- `obj_val`: Objective value of the relaxation
- `ecp_status`: Status of the ECP at termination
- `initial_lp_obj`: Starting LP objective
- `iter_ecp`: Number of ECP iteration executed
- `total_cuts`: Total number of added DD cuts
- `ecp_time`: Total ECP elapsed time
- `width`: Array of widths of DDs in the model
- `dd_time`: Time to construct DD (and form root node)
- `relative_obj_improvement`: relative obj improvement over the last few iterations at termination
"""
function write_pricing(filename::AbstractString, time_limit::Float64, iter_max_ecp::Int, max_cut_num::Int, lp_to_ip_iter::Int, lp_to_ip_tol::Float64, obj_val::Float64, obj_improvement_tol::Float64, tol_num_ecp::Int,
	ecp_status::Symbol, initial_lp_obj::Float64, iter_ecp::Int, total_cuts::Int, ecp_time::Float64, width::Array{Int, 1}, dd_time::Float64, relative_obj_improvement::Float64)

	open(filename, "w") do f
		write(f, "This file contains the result of ECP for pricing problem\n")
		write(f, "********************************************************\n\n")
		write(f, "Starting objective value:\t$initial_lp_obj\n\n")
		write(f, "Final objective value:\t$obj_val\n\n")
		write(f, "ECP status:\t$ecp_status\n\n")
		write(f, "DD construction time:\t$dd_time\n\n")
		write(f, "ECP time:\t$ecp_time\n\n")
		write(f, "Time limit:\t$time_limit\n\n")
		write(f, "Relative obj improvement tolerance:\t$obj_improvement_tol\n\n")
		write(f, "Number of last few iterations to compute relative obj improvement:\t$tol_num_ecp\n\n")
		write(f, "Relative obj improvement at termination:\t$relative_obj_improvement\n\n")
		write(f, "Number of ECP iterations:\t$iter_ecp\n\n")
		write(f, "Number of added DD cuts:\t$total_cuts\n\n")
		write(f, "ECP iteration limit:\t$iter_max_ecp\n\n")
		write(f, "Max cut number at each iteration:\t$max_cut_num\n\n")
		write(f, "DD widths:\n")
		for i=1:length(width)
			write(f, "$(width[i])\t")
		end
	end
end
