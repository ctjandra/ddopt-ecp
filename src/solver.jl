export DDECPSolver

# Solver object to indicate a solver is not set
type UnsetSolver <: MathProgBase.AbstractMathProgSolver
end

# DD-ECP solver
type DDECPSolver <: MathProgBase.AbstractMathProgSolver
    lp_solver::MathProgBase.AbstractMathProgSolver   # LP solver
end

function DDECPSolver(;
    lp_solver = UnsetSolver(),
    )

    if lp_solver == UnsetSolver()
        error("LP solver not specified; set lp_solver")
    end
end
