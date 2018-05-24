include("read_write.jl")
include("pricing_specs.jl")
include("../../src/Dexter.jl")
using Dexter

using LightGraphs
using Gadfly
using GraphPlot
using JuMP
using Clp

# read problem instance and construct the model
filename = joinpath(@__DIR__, "instances", "n100_c10_u10_s1.prc")
specs = read_pricing(filename)
m = create_pricing_model(specs)

print(m)

# TODO: Right now we are assuming <=. Consider >= case.

# Compute evaluations of parts of additively separable constraints
evals = Dexter.evaluate_separable_constraints(m)

# TODO: Transition function should support dynamic ordering

# Set ordering and objective
n = MathProgBase.numvar(m)
ordering = Dexter.NoOrdering()

# Create specs for decision diagram
eval = evals[7]
initial_state = Dexter.create_constraint_initial_state(eval)
domains = Dexter.create_constraint_domain_function(m)
transition = Dexter.create_constraint_transition_function(eval, domains, ordering)
mathprog_specs = Dexter.ProblemSpecs(initial_state, transition, domains)

println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

# Construct decision diagram
dd = Dexter.construct_DD(n, mathprog_specs, ordering=ordering)

println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")

# Plot decision diagram
locs_x = [Dexter.get_node_idlayer(dd, i) * 20.0 for i in 1:nv(dd.graph)]
locs_y = [Dexter.get_node_layer(dd, i) / 1.0 for i in 1:nv(dd.graph)]
gplot(dd.graph, locs_x, locs_y) |> PNG("test_mp1.png", 20cm, 50cm)
