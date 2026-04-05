function run_all_visualizations(output_dir)
% run_all_visualizations - 运行所有可视化脚本
% 消融实验 -> figures/ablation
% 对比实验 -> figures/comparison
fprintf('========== Running All Visualizations ==========\n\n');

fprintf('1. Generating Ablation Charts...\n');
try
    plot_ablation_all();
    fprintf('   -> Done\n');
catch ME
    fprintf('   -> Error: %s\n', ME.message);
end

fprintf('\n2. Generating Comparison Charts...\n');
try
    plot_comparison_all();
    fprintf('   -> Done\n');
catch ME
    fprintf('   -> Error: %s\n', ME.message);
end

fprintf('\n3. Generating Deployment Topology...\n');
try
    plot_deployment_topology();
    fprintf('   -> Done\n');
catch ME
    fprintf('   -> Error: %s\n', ME.message);
end

fprintf('\n========== All Visualizations Complete ==========\n');
fprintf('Output directories:\n');
fprintf('  - figures/ablation (消融实验)\n');
fprintf('  - figures/comparison (对比实验)\n');
fprintf('  - figures (部署拓扑图)\n');
end