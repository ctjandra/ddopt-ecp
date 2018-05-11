include("../../../src/Dexter.jl")
using Dexter

using LightGraphs
using Gadfly
using GraphPlot
using JuMP
using Clp

m = Model(solver = ClpSolver())
@variable(m, 0 <= x <= 5)
@variable(m, 0 <= y <= 5)
@variable(m, 0 <= z <= 5)

@objective(m, Max, 5x + 3y + 2z)
# @NLconstraint(m, cos(x)/y * x * e^(-x^2) <= 3.0 )
# @NLconstraint(m, cos(x) * e^(-x^2) + y - x*x + 1 <= 5.0 - x^2)
# @NLconstraint(m, cos(x) * e^(-x^2) + y - x*x + -(x/2 + 2*y) + 1 <= x^2 + 5)
# @NLconstraint(m, 2x <= 3.0 )
# @NLconstraint(m, x + y + z <= 8)
# @NLconstraint(m, - 2 * x * e^(-x) - 5 * y * e^(-y^2) - 3 * z * e^(-z) <= 30)
@NLconstraint(m, 2 * x * e^(-x) + 5 * y * e^(-y^2) + 3 * z * e^(-z) <= 3)

print(m)

# TODO: Right now we are assuming <=. Consider >= case.

# Compute evaluations of parts of additively separable constraints
evals = Dexter.evaluate_separable_constraints(m)

# TODO: Transition function should support dynamic ordering

# Set ordering and objective
n = MathProgBase.numvar(m)
ordering = Dexter.NoOrdering()

# Create specs for decision diagram
eval = evals[1]
initial_state = Dexter.create_constraint_initial_state(eval)
domains = Dexter.create_constraint_domain_function(m)
transition = Dexter.create_constraint_transition_function(eval, domains, ordering)
mathprog_specs = Dexter.ProblemSpecs(initial_state, transition, domains)

# Construct decision diagram
dd = Dexter.construct_DD(n, mathprog_specs, ordering=ordering)

# Plot decision diagram
locs_x = [Dexter.get_node_idlayer(dd, i) * 1.0 for i in 1:nv(dd.graph)]
locs_y = [Dexter.get_node_layer(dd, i) / n for i in 1:nv(dd.graph)]
gplot(dd.graph, locs_x, locs_y) |> PNG("test_mp.png", 16cm, 16cm)
