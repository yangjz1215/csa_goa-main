function run_all_visualizations(output_dir)
if nargin < 1
    output_dir = fullfile('..', 'figures');
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('========== Running All Visualizations ==========\n\n');

fprintf('1. Generating Pareto Comparison (Figures A & B)...\n');
try
    plot_pareto_comparison([], output_dir);
    fprintf('   -> Done\n');
catch ME
    fprintf('   -> Error: %s\n', ME.message);
end

fprintf('\n2. Generating Convergence & HV Curves (Figure C)...\n');
try
    plot_convergence_hv([], output_dir);
    fprintf('   -> Done\n');
catch ME
    fprintf('   -> Error: %s\n', ME.message);
end

fprintf('\n3. Generating Ablation Results (Figure D)...\n');
try
    plot_ablation_results([], output_dir);
    fprintf('   -> Done\n');
catch ME
    fprintf('   -> Error: %s\n', ME.message);
end

fprintf('\n4. Generating Boxplots (Figure E)...\n');
try
    plot_boxplot([], output_dir);
    fprintf('   -> Done\n');
catch ME
    fprintf('   -> Error: %s\n', ME.message);
end

fprintf('\n5. Generating Statistical Tests (Table F)...\n');
try
    plot_statistical_tests([], output_dir);
    fprintf('   -> Done\n');
catch ME
    fprintf('   -> Error: %s\n', ME.message);
end

fprintf('\n6. Generating Deployment Topology (Figure G)...\n');
try
    plot_deployment_topology([], [], output_dir);
    fprintf('   -> Done\n');
catch ME
    fprintf('   -> Error: %s\n', ME.message);
end

fprintf('\n========== All Visualizations Complete ==========\n');
fprintf('Output directory: %s\n', fullfile(pwd, output_dir));
end