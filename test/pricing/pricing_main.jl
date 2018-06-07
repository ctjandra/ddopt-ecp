include("read_write.jl")
include("pricing_specs.jl")
include("../../src/Dexter.jl")
#include("../../src/dd/construction.jl")

using Dexter

using LightGraphs
using Gadfly
using GraphPlot
using JuMP
using Clp
using CPLEX

# read problem instance and construct the model
filename = joinpath(@__DIR__, "instances", "n100_c10_u10_s1.prc")
specs = read_pricing(filename)
m = create_pricing_model(specs)

print(m)

n = specs.var_num
constr_num = specs.pricing_constr_num
dd_list = Array{Dexter.DecisionDiagram}(constr_num)

# Determine constraint decomposition, ordering and variable domains
evals = Dexter.evaluate_separable_constraints(m)
domains = Dexter.create_constraint_domain_function(m)
ordering = Dexter.NoOrdering()

# Create specs for decision diagram
for i=1:constr_num

    eval = evals[i]
    initial_state = Dexter.create_constraint_initial_state(eval)
    transition = Dexter.create_constraint_transition_function(eval, domains, ordering)
    mathprog_specs = Dexter.ProblemSpecs(initial_state, transition, domains)

    println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

    # Construct decision diagram
    dd_list[i] = Dexter.construct_DD(n, mathprog_specs, ordering=ordering)

    println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")

    # Plot decision diagram
    #locs_x = [Dexter.get_node_idlayer(dd, i) * 20.0 for i in 1:nv(dd.graph)]
    #locs_y = [Dexter.get_node_layer(dd, i) / 1.0 for i in 1:nv(dd.graph)]
    #gplot(dd.graph, locs_x, locs_y) |> PNG("test_mp1.png", 20cm, 50cm)

    # Print the layer widths in file
    #fname = "Constraint$i"
    #f = open(fname, "w")
    #   for j=1:n+1
    #       write(f, "Layer $j: $(length(dd_list[i].layers[j]))\n")
    #   end
    #close(f)

end

# Initializing parameters
time_lim = 100


# Build the initial model
m2 = Model(solver = CplexSolver(CPX_PARAM_TILIM = time_lim))

@variable(m2, x[i=1:n], lowerbound = specs.lb[i], upperbound = specs.ub[i])

obj = AffExpr(x, specs.obj_c, 0.0)				#the objective function is linear
@objective(m2, Max, obj)


# Running ECP algorithm
m3, x_val, obj_val, ecp_status = Dexter.ECP(n, m2, dd_list, evals, time_limit = time_lim)

print("\n\n\n")
println("Objective value: ", obj_val)
println("ECP status: ", ecp_status)
