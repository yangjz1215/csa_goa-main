function table_str = plot_statistical_tests(data_file, output_dir)
if nargin < 1 || isempty(data_file)
    data_file = fullfile('..', 'experiments', 'ablation_results_para_Map1_Medium_20260405_184023.mat');
end
if nargin < 2 || isempty(output_dir)
    output_dir = fullfile('..', 'figures');
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

load(data_file);

variants = {'proposed', 'no_subpop', 'no_adaptive', 'no_levy', 'no_stop'};
variant_labels = {'Proposed', 'w/o Multi-Subpop', 'w/o Adaptive', 'w/o E-Levy', 'w/o MO-SmartStop'};

if ~isfield(results, 'proposed')
    fprintf('Error: proposed variant not found in data\n');
    table_str = '';
    return;
end

proposed_hv = results.proposed.hv_values;

n_variants = length(variants);
p_values = zeros(n_variants, 1);
h_stats = zeros(n_variants, 1);
mean_improvement = zeros(n_variants, 1);
effect_sizes = zeros(n_variants, 1);

for v = 1:n_variants
    variant = variants{v};
    if isfield(results, variant) && strcmp(variant, 'proposed')
        p_values(v) = 1.0;
        h_stats(v) = 0;
        mean_improvement(v) = 0;
        effect_sizes(v) = 0;
    elseif isfield(results, variant)
        variant_hv = results.(variant).hv_values;
        [p_values(v), h_stats(v)] = signrank(proposed_hv, variant_hv);
        mean_improvement(v) = (mean(proposed_hv) - mean(variant_hv)) / mean(variant_hv) * 100;
        diffs = proposed_hv - variant_hv;
        effect_sizes(v) = mean(diffs) / std(diffs);
    else
        p_values(v) = NaN;
        h_stats(v) = NaN;
        mean_improvement(v) = NaN;
        effect_sizes(v) = NaN;
    end
end

fig = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.75, 0.4]);
set(gcf, 'Color', 'w');

col_labels = {'Variant', 'Mean HV', 'Std HV', 'Improvement (%)', 'Effect Size', 'p-value', 'Significant?'};
table_data = {};
row_names = {};

for v = 2:n_variants
    variant = variants{v};
    row_names = [row_names, variant_labels{v}];
    if isfield(results, variant)
        mean_val = mean(results.(variant).hv_values);
        std_val = std(results.(variant).hv_values);
        imp_val = mean_improvement(v);
        es_val = effect_sizes(v);
        p_val = p_values(v);
        h_val = h_stats(v);
    else
        mean_val = NaN; std_val = NaN; imp_val = NaN; es_val = NaN; p_val = NaN; h_val = NaN;
    end

    if h_val == 1
        sig_str = 'Yes (*, p<0.05)';
    elseif ~isnan(p_val) && p_val < 0.10
        sig_str = 'Marginal (+, p<0.10)';
    else
        sig_str = 'No';
    end

    if isnan(mean_val)
        row = {variant_labels{v}, 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A'};
    else
        row = {variant_labels{v}, sprintf('%.4f', mean_val), sprintf('%.4f', std_val), ...
            sprintf('%.2f', imp_val), sprintf('%.3f', es_val), sprintf('%.4f', p_val), sig_str};
    end
    table_data = [table_data; row];
end

uitable('Data', table_data, ...
    'ColumnName', col_labels, ...
    'RowName', row_names, ...
    'FontSize', 9, ...
    'FontName', 'Times New Roman', ...
    'BackgroundColor', [1, 1, 1], ...
    'ForegroundColor', [0, 0, 0], ...
    'Units', 'normalized', ...
    'Position', [0.02, 0.08, 0.96, 0.88]);

title({'Wilcoxon Signed-Rank Test (Paired): Proposed vs Ablation Variants'; ...
    sprintf('Proposed HV: %.4f +/- %.4f | n=30 runs', mean(proposed_hv), std(proposed_hv))}, ...
    'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

table_str = sprintf('=== Wilcoxon Signed-Rank Test Results ===\n');
table_str = [table_str, sprintf('Baseline (Proposed): HV=%.4f +/- %.4f\n\n', mean(proposed_hv), std(proposed_hv))];
table_str = [table_str, sprintf('%-20s | %10s | %10s | %12s | %10s | %10s | %s\n', ...
    'Variant', 'Mean HV', 'Std HV', 'Improvement', 'Effect Sz', 'p-value', 'Sig?')];
table_str = [table_str, repmat('-', 1, 100), '\n'];
for v = 2:n_variants
    if isfield(results, variants{v})
        if h_stats(v) == 1
            ss = '*';
        elseif p_values(v) < 0.10
            ss = '+';
        else
            ss = '-';
        end
        table_str = [table_str, sprintf('%-20s | %10.4f | %10.4f | %11.2f%% | %10.3f | %10.4f | %s\n', ...
            variant_labels{v}, mean(results.(variants{v}).hv_values), ...
            std(results.(variants{v}).hv_values), mean_improvement(v), ...
            effect_sizes(v), p_values(v), ss)];
    end
end
table_str = [table_str, '\n* p<0.05 significant, + p<0.10 marginal\n'];

fprintf('\n%s\n', table_str);

saveas(fig, fullfile(output_dir, 'statistical_tests.fig'));
saveas(fig, fullfile(output_dir, 'statistical_tests.png'));
fprintf('Statistical tests saved to %s\n', output_dir);
close(fig);
end