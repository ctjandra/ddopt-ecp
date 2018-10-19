include("read_write.jl")
include("pricing_specs.jl")
include("../../src/Dexter.jl")
#include("../../src/dd/construction.jl")

using Dexter

using LightGraphs
#using Gadfly
using GraphPlot
using JuMP
using Clp
using CPLEX

# read problem instance and construct the model
filename = joinpath(@__DIR__, "instances", "C_Input_Mon_s1.txt")
specs = read_pricing(filename)
m1 = create_pricing_model(specs)
m2 = create_pricing_LP(specs)


# This part computes the functional property of separable functions
aux_p = Array{Array{Array{Dexter.UVFuncProp, 1}, 1}, 1}(specs.pricing_constr_num)
for i=1:specs.pricing_constr_num
    aux_p[i] = Array{Array{Dexter.UVFuncProp, 1}, 1}(specs.var_num)
    for j=1:specs.var_num
        a = specs.ub[j]/(specs.constr_d[i,j]^(1/specs.constr_d[i,j]))
        rng1 = Dexter.UVFuncProp(0, 0, 1, 0, a)
        rng2 = Dexter.UVFuncProp(0, 0, 1, a, Inf)
        aux_p[i][j] = [rng1, rng2]
    end
end


bb_tree = Dexter.create_BB_tree(m1, m2)
root_time = @elapsed Dexter.create_BB_root!(bb_tree, prop = aux_p)


# Initializing parameters
time_limit0 = 200.0
iter_max_ecp0 = 1000        # max limit on ecp iterations
max_cut_num0 = 3            # max cut number added at EACH iteration of ecp
lp_to_ip_iter0 = 500         # number of ecp iterations to solve LP before switching to IP (if applicable)
lp_to_ip_tol0 = 0.00001      # constraint violation tolerance below which ecp switches to solving IP (if applicable)
tol_num_ecp0 = 10                # number of last iterations over which the relative objective improvement is computed
obj_improvement_tol0 = 1e-4         # relative objective value improvememtn over the last few (tol_num_ecp) iterations

# Running ECP algorithm
x_val, obj_val, ecp_status, initial_lp_obj, iter_ecp, total_cuts, ecp_time, relative_obj_improvement =
Dexter.ECP!(bb_tree.tree[1].model_info.m, bb_tree.tree[1].dd_info.dd, bb_tree.tree[1].dd_info.con_function, bb_tree.tree[1].dd_info.initial,
time_limit = time_limit0, iter_max_ecp = iter_max_ecp0, max_cut_num = max_cut_num0,
lp_to_ip_iter = lp_to_ip_iter0, lp_to_ip_tol = lp_to_ip_tol0, obj_improvement_tol = obj_improvement_tol0, tol_num_ecp = tol_num_ecp0)

# computing DD widths for all DDs
width = Array{Int, 1}(specs.pricing_constr_num)
for i=1:specs.pricing_constr_num
    width[i] = Dexter.get_width(bb_tree.tree[1].dd_info.dd[i])
end

# write results in txt file
fname = joinpath(@__DIR__, "instances", "Dexter_Result_Int_s1.txt")
write_pricing(fname, time_limit0, iter_max_ecp0, max_cut_num0, lp_to_ip_iter0, lp_to_ip_tol0, obj_val, obj_improvement_tol0, tol_num_ecp0,
 ecp_status, initial_lp_obj, iter_ecp, total_cuts, ecp_time, width, root_time, relative_obj_improvement)

#print("\n\n\n")
#println("Objective value: ", obj_val)
#println("ECP status: ", ecp_status)
