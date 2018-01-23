# TODO list

## Implement decision diagrams [CT]

We should allow problem flexibility (i.e. not write the code specifically for additively separable constraints); see previous implementations of decision diagrams.

* Class structures
  * DP formulation structure
  * Decision diagram, nodes, etc.
  * Decision diagram construction
* Top-down construction (doing it first because it is simpler)
* Ordering
* Converter that takes an additively separable JuMP constraint and returns a DP formulation
* Depth-first construction from the DP formulation
  * Issue: can the DP formulation provide equivalence class information for the depth-first construction?


## Implement cut generation procedure

Given a decision diagram, return a JuMP constraint consisting of the DD cut.

> CT: For now, assume the decision diagram is a graph following the [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl) module; we'll rearrange the structure later. This will also make it easier to test your code independently from the decision diagram part, as you can create a small example and add edges yourself.

* Cut generating LP
  * Use the LP solver given by the user (which is of type MathProgBase.AbstractMathProgSolver)
* Projected subgradient algorithm


## Implement branching

* Variable selection
* Node selection
* Perhaps borrow ideas of class design for the branching tree from SCIP?


## Create test set

A test set would be a set of instances with additively separable constraints in JuMP model format.
