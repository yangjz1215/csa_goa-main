function metrics = calculate_all_metrics(obtained_front, true_front, reference_point)
    if nargin < 3
        reference_point = [];
    end
    if nargin < 2
        true_front = [];
    end

    metrics = struct();

    if ~isempty(obtained_front)
        if isempty(reference_point)
            if isempty(true_front)
                reference_point = max(obtained_front, [], 1) * 1.1;
            else
                reference_point = max([obtained_front; true_front], [], 1) * 1.1;
            end
        end
        metrics.hv = hypervolume(obtained_front, reference_point);
    else
        metrics.hv = 0;
    end

    if ~isempty(obtained_front) && ~isempty(true_front)
        metrics.igd = igd(obtained_front, true_front);
        metrics.gd = gd(obtained_front, true_front);
    else
        metrics.igd = nan;
        metrics.gd = nan;
    end

    if ~isempty(obtained_front)
        metrics.spread = spread(obtained_front, true_front);
    else
        metrics.spread = nan;
    end

    if ~isempty(obtained_front)
        metrics.n_solutions = size(obtained_front, 1);
    else
        metrics.n_solutions = 0;
    end
end
