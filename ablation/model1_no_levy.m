function [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = model1_no_levy(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities)
    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
        cSA_GOA_main_ablation(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities, 'no_levy');
end
