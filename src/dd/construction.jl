include_dependency("ordering.jl")
include_dependency("core.jl")
include("mathprog/function_analysis.jl")
include("mathprog/decomposition.jl")


"""
Compute transition function for additively separable constraints.

# Input:
- `evals`: An array of decomposed terms of a constraint.
- `func`: An array of univariate functions that stores the separate functions and their properties for each variable (array elements must be undefined as input)
- `prop`: Keyword that contains an initial assessment for properties of the univariate functions (optional)
"""
# ultimately, we would want this function to decompose the separable functions in evals into an algebraic functional form
function create_constraint_transition_function!(func::Array{UVFunc, 1}, evals::SeparableFunctionEvaluation; prp::Array{Array{UVFuncProp, 1}, 1} = [])
	@assert(length(evals.separate_functions) == length(func))
	for i=1:length(evals.separate_functions)
		if isassigned(evals.separate_functions, i)
			# if no initial information is given for the function property
			if !isassigned(prp, i)
				#TODO: run a functional analysis, to obtain function properties on the domain
				a::UVFuncProp = UVFuncProp(false, false, false, -Inf, Inf)		# sets a general functional property
				func[i] = UVFunc(evals.separate_functions[i], [a])
			# otherwise, the given info is passed as the properties
			else
				func[i] = UVFunc(evals.separate_functions[i], prp[i])		#NOTE: prop[i] is an array, so changing its elements will change func[i].prop and vice-versa
			end
		end
	end
	return
end



"""
Compute initial state for additively separable constraints.

# Input:
- `eval`: Decomposed terms of a constraint.

"""
function create_constraint_initial_state(eval::SeparableFunctionEvaluation)::Float64
    return eval.constant  	# constant value of the constraint as moved to lhs of <= constraints
end


"""
Compute single function representing the variable part of constraints.

# Input:
- `eval`: Decomposed terms of a constraint.
"""
function create_constraint_single_function(eval::SeparableFunctionEvaluation)::Function
    return eval.single_function 	# variable part of the constraint as moved to lhs of <= constraints
end



"""
Compute domain range for additively separable constraints based on a given math prog model.

# Input:
- `domain`: store an array whose elements return a tuple (lowerbound, upperbound, variable type).
- `m`: math prog model.
"""
function create_constraint_domain_function!(domain::Array{Tuple{Float64, Float64, Symbol}, 1}, m::JuMP.Model)
	nvars::Int64 = MathProgBase.numvar(m)
	@assert(length(domain) == nvars)
	for v in 1:nvars
		v_obj = Variable(m, v)
		v3 = getcategory(v_obj)		# available categories: :Cont, :Int, :Bin
		if v3 == :Bin				# for now, handle Binary as Integer
			v3 == :Int
		end
		v1 = getlowerbound(v_obj)
		v2 = getupperbound(v_obj)
		domain[v] = (v1, v2, v3)
	end
end




"""
Construct a decision diagram.

# Input
- `func`: Array of transition functions for each variable.
- `initial`: Initial state value.
- `domain`: Array of variable domains
- `lb_ub`: Array of tuples that contains (lower bound, upper bound) of univariate functions over the domain of variables defined for this DD
- `ordering`: Ordering of the variables in DD.
- `reduced_arc`: Determines whether parallel arcs are merged to their lower and upper bound label values
"""
function construct_DD(func::Array{UVFunc, 1},  initial::Float64, domain::Array{Tuple{Float64, Float64, Symbol}, 1}, lb_ub::Array{Tuple{Float64, Float64}, 1}; ordering::Ordering = NoOrdering(), reduced_arc::Bool = false)::DecisionDiagram

	n::Int = length(func)	# number of variables
	# TODO
	# - This is not a top-down search, but a simple depth-first one
	# - Implement long arcs
	# - Implement mapping between problem variables and DD variables
	# - Implement reduction at the end

	# preprocessing the min max values function at each level
	lb_ahead = Array{Float64}(  n)		# For each variable, stores the lb value that can be obtained starting from the next layer to the terminal
	ub_ahead = Array{Float64}(  n)		# For each variable, stores the ub value that can be obtained starting from the next layer to the terminal
	lb_ahead[get_var(ordering,n)] = 0		#initialize values for the last layer
	ub_ahead[get_var(ordering,n)] = 0
	for i::Int in n-1:-1:1
		var::Int = get_var(ordering,i)		# gets the variable of current layer i
		var_next::Int = get_var(ordering,i+1)		# gets the variable in layer i+1
		if isassigned(func, var_next)		#checks if the function is defined for the variable
			lb_ahead[var] = lb_ub[var_next][1] + lb_ahead[var_next]
			ub_ahead[var] = lb_ub[var_next][2] + ub_ahead[var_next]
		else
			lb_ahead[var] = lb_ahead[var_next]	#if the function is not defined for the next layer, it does not affect the lb ub values
			ub_ahead[var] = ub_ahead[var_next]
		end
	end

	S = Float64
	# constructing dd
	dd::DecisionDiagram{S} = DecisionDiagram(n)

	# TODO For top-down, this should probably only be a single Dict
	layer_states::Array{Dict{S, Int}, 1} = [Dict{S, Int}() for i=1:n+1]		# stores pairs of (state value, node id) at each layer

	# Add source node
	source::Int = add_node!(dd, 1, initial)	#the index of the created node
	layer_states[1][initial] = source

	unexplored_nodes::Array{Tuple{Int, Int}, 1} = Array{Tuple{Int, Int}, 1}(  0)	# a node tuple contains (layer, node id)
	push!(unexplored_nodes, (1, source))
	niters::Int = 0

	# stores an array of possible arc labels with their computed lb value for each variable
	# this saves multiple evaluation of the variable function at similar values
	#NOTE: this can only be used, if the same grid/values of arc labels will be evaluated as outgoing arcs of all nodes at each layer.
	memo::Array{Array{Float64, 1}, 1} = [Array{Float64, 1}() for i=1:n]

	while !isempty(unexplored_nodes)

		# Read next unexplored node
		node::Tuple{Int, Int} = pop!(unexplored_nodes)
		niters += 1
		layer::Int = node[1]
		id::Int = node[2]

		if layer > n
			continue # terminal node
		end

		state::S = get_node_state(dd, layer, id)
		var = get_var(ordering, layer)

		# println("Layer $(layer): exploring $(node) with state $(state)")

		# Create children
		# if the variable is integer
		if domain[var][3] == :Int
			v1 = ceil(Int, domain[var][1])
			v2 = floor(Int, domain[var][2])
			for val0::Int in v1:v2
				new_state::S = state
				val = float(val0)			# converts Int values to float
				if isassigned(func, var)
					#if !isassigned(memo[var], val0-v1+1)		# the label has not been evaluated at the function yet
					#	@assert(length(memo[var]) == val0-v1)
					#	push!(memo[var], func[var].f(val))
					#end
					#new_state += memo[var][val0-v1+1]		# computes the new state value
					new_state += func[var].f(val)
				end

				#To determine feasibility, it is assumed that the constraint form is f(x) <= 0
				if new_state <= -ub_ahead[var]			#new state is feasible no matter what variable assignments are chosen from now on
					new_state = -ub_ahead[var]
				elseif new_state > -lb_ahead[var]			#new state is infeasible no matter what variable assignments are chosen
					continue							# do nothing and go to the next val in the loop
				end

				# Relaxation step
				#new_state = floor(new_state/10)
				new_state = floor(new_state)
				#new_state = floor(10*new_state)/10		# rounds down the fractional state values: ONLY use to construct RELAXED DD when state values are fractional, e.g., for the pricing problem

				next_layer::Int = layer + 1
				states_to_nodes::Dict{S, Int} = layer_states[next_layer]

				# Equivalence check
				if haskey(states_to_nodes, new_state)
					# Equivalent node found: update it
					t::Int = states_to_nodes[new_state]		#the index of the node with that state value
				else
					# No equivalent node: create new child
					t = add_node!(dd, next_layer, new_state)
					states_to_nodes[new_state] = t
					# Insert new node into list of unexplored nodes
					push!(unexplored_nodes, (next_layer,t))
				end

				if reduced_arc
					add_reduced_arc!(dd, layer, id, next_layer, t, val)
				else
					add_arc!(dd, layer, id, next_layer, t, val)
				end
				# println("Added $val-arc from $node to $t")
			end
		# if the variable is continuous
		elseif domain[var][3] == :Cont
			# TODO: implement an elaborate grid desgin (uniform, or non-uniform)
			dom_eps::Float64 = 1.00e-8
			if domain[var][2] - domain[var][1] <= dom_eps	# if the range is very small, then pick a single grid
				grid = 1
			else
				grid::Int = 30					# the number of intervals to divide the continuous range
			end
			frac::Float64 = (domain[var][2] - domain[var][1])/grid
			interval_lb = domain[var][1]
			for i=1:grid
				interval_ub = interval_lb + frac
				#println("\n", interval_lb, "\t", frac, "\n")
				new_state::S = state
				if isassigned(func, var)
					if !isassigned(memo[var], i)		# the label has not been evaluated at the function yet
						@assert(length(memo[var]) == i-1)
						# TODO: find the properties of the function over the new interval.
						# compute the lb of the function over the interval
						func_lb = uvf_lb(func[var].f, func[var].prop, (interval_lb, interval_ub, :Cont))
						push!(memo[var], func_lb)
					end
					new_state += memo[var][i]		# computes the new state value (for merged nodes on the inverval)
				end

				#To determine feasibility, it is assumed that the constraint form is f(x) <= 0
				if new_state <= -ub_ahead[var]			#new state is feasible no matter what variable assignments are chosen from now on
					new_state = -ub_ahead[var]
				elseif new_state > -lb_ahead[var]			#new state is infeasible no matter what variable assignments are chosen
					interval_lb = interval_ub
					continue							# do nothing and go to the next val in the loop
				end

				# Relaxation step
				#new_state = floor(new_state/10)
				new_state = floor(new_state)
				#new_state = floor(10*new_state)/10		# rounds down the fractional state values: ONLY use to construct RELAXED DD when state values are fractional, e.g., for the pricing problem

				next_layer::Int = layer + 1
				states_to_nodes::Dict{S, Int} = layer_states[next_layer]

				# Equivalence check
				if haskey(states_to_nodes, new_state)
					# Equivalent node found: update it
					t::Int = states_to_nodes[new_state]		#the index of the node with that state value
				else
					# No equivalent node: create new child
					t = add_node!(dd, next_layer, new_state)
					states_to_nodes[new_state] = t
					# Insert new node into list of unexplored nodes
					push!(unexplored_nodes, (next_layer,t))
				end

				if reduced_arc
					add_reduced_arc!(dd, layer, id, next_layer, t, interval_lb)
					add_reduced_arc!(dd, layer, id, next_layer, t, interval_ub)
				else
					# add both lb and ub label arcs between the nodes
					add_arc!(dd, layer, id, next_layer, t, interval_lb)
					add_arc!(dd, layer, id, next_layer, t, interval_ub)
				end
				# println("Added $val-arc from $node to $t")

				interval_lb = interval_ub
			end

		end
	end

	return dd
end
