function plot_deployment_topology(result_file, map_file, output_dir)
if nargin < 1
    result_file = fullfile('..', 'experiments', 'ablation_results_para_Map1_Medium_20260405_184023.mat');
end
if nargin < 2
    map_file = fullfile('..', 'maps', 'Map1_Medium.mat');
end
if nargin < 3
    output_dir = fullfile('..', 'figures');
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

load(map_file);
load(result_file);

N_User_actual = size(User, 1);
if N_User_actual == 200
    config_suffix = 'Small';
elseif N_User_actual == 500
    config_suffix = 'Medium';
else
    config_suffix = 'Large';
end

config_file = fullfile('..', 'maps', ['Map_', config_suffix, '_Config.mat']);
if exist(config_file, 'file')
    config_data = load(config_file);
    RRH = config_data.RRH;
    RRH_type = config_data.RRH_type;
else
    error('Config file not found');
end

fig = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.85, 0.45]);
set(gcf, 'Color', 'w');

high_priority_users = User(priorities >= 3, :);
normal_users = User(priorities < 3, :);

subplot(1, 2, 1);
hold on;

scatter(normal_users(:,1), normal_users(:,2), 30, [0.6, 0.6, 0.6], 'filled', 'MarkerFaceAlpha', 0.4);

scatter(high_priority_users(:,1), high_priority_users(:,2), 60, [1, 0, 0], 'filled', 'Marker', '^', 'LineWidth', 1);

rrh_normal = RRH(RRH_type == 0, :);
rrh_enhanced = RRH(RRH_type == 1, :);
scatter(rrh_normal(:,1), rrh_normal(:,2), 100, [0.3, 0.3, 0.3], 'square', 'filled', 'LineWidth', 2);
scatter(rrh_enhanced(:,1), rrh_enhanced(:,2), 120, [0, 0.5, 0], 'square', 'filled', 'LineWidth', 2);

if isfield(results, 'proposed') && ~isempty(results.proposed.pareto_fronts)
    all_pf = results.proposed.pareto_fronts{1};
    if ~isempty(all_pf)
        [~, max_idx] = max(all_pf(:,1));
        [~, min_idx] = min(all_pf(:,2));

        if max_idx ~= min_idx
            uav_max_cov = squeeze(results.proposed.final_uavs{max_idx});
            uav_min_eng = squeeze(results.proposed.final_uavs{min_idx});
        else
            uav_max_cov = squeeze(results.proposed.final_uavs{max_idx});
            uav_min_eng = uav_max_cov;
        end

        scatter(uav_max_cov(:,1), uav_max_cov(:,2), 300, [1, 0, 0], 'o', 'filled', 'LineWidth', 2.5);
        for i = 1:size(uav_max_cov, 1)
            circle(uav_max_cov(i,1), uav_max_cov(i,2), 150, 'Color', [1, 0, 0], 'LineWidth', 1, 'LineStyle', '--');
        end
    end
end

xlabel('X (m)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Y (m)', 'FontSize', 11, 'FontWeight', 'bold');
title('(a) Maximum Coverage Solution', 'FontSize', 13, 'FontWeight', 'bold');
legend({'Normal Users', 'High-Priority Users', 'RRH', 'Enhanced RRH', 'UAVs', 'Coverage Radius'}, ...
    'Location', 'best', 'FontSize', 9);
axis([0 1000 0 1000]);
axis square;
grid on;

subplot(1, 2, 2);
hold on;

scatter(normal_users(:,1), normal_users(:,2), 30, [0.6, 0.6, 0.6], 'filled', 'MarkerFaceAlpha', 0.4);
scatter(high_priority_users(:,1), high_priority_users(:,2), 60, [1, 0, 0], 'filled', 'Marker', '^', 'LineWidth', 1);

scatter(rrh_normal(:,1), rrh_normal(:,2), 100, [0.3, 0.3, 0.3], 'square', 'filled', 'LineWidth', 2);
scatter(rrh_enhanced(:,1), rrh_enhanced(:,2), 120, [0, 0.5, 0], 'square', 'filled', 'LineWidth', 2);

if exist('uav_min_eng', 'var')
    scatter(uav_min_eng(:,1), uav_min_eng(:,2), 300, [0, 0.4, 0.8], 'o', 'filled', 'LineWidth', 2.5);
    for i = 1:size(uav_min_eng, 1)
        circle(uav_min_eng(i,1), uav_min_eng(i,2), 150, 'Color', [0, 0.4, 0.8], 'LineWidth', 1, 'LineStyle', '--');
    end
end

xlabel('X (m)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Y (m)', 'FontSize', 11, 'FontWeight', 'bold');
title('(b) Minimum Energy Solution', 'FontSize', 13, 'FontWeight', 'bold');
legend({'Normal Users', 'High-Priority Users', 'RRH', 'Enhanced RRH', 'UAVs', 'Coverage Radius'}, ...
    'Location', 'best', 'FontSize', 9);
axis([0 1000 0 1000]);
axis square;
grid on;

annotation('textbox', [0.02, 0.02, 0.3, 0.08], 'String', ...
    {'UAV Deployment Topology (Map1, Medium Scale)', ...
     'Circle: UAV coverage area | Triangle: High-priority users'}, ...
    'FontSize', 9, 'EdgeColor', 'none', 'BackgroundColor', [1, 1, 1]);

saveas(fig, fullfile(output_dir, 'deployment_topology.fig'));
saveas(fig, fullfile(output_dir, 'deployment_topology.png'));
fprintf('Deployment topology saved to %s\n', output_dir);
close(fig);

    function h_circle = circle(x, y, r, varargin)
        theta = linspace(0, 2*pi, 100);
        h_circle = plot(x + r*cos(theta), y + r*sin(theta), varargin{:});
    end
end