function [bx, by, bstd, slope, x_fit, y_fit] = compute_binned_stats(x_data, y_data, x_min, x_max, num_bins, std_thresh)
    if nargin < 3 || isempty(x_min), x_min = min(x_data(:), [], 'omitnan'); end
    if nargin < 4 || isempty(x_max), x_max = max(x_data(:), [], 'omitnan'); end
    if nargin < 5 || isempty(num_bins), num_bins = 20; end
    if nargin < 6 || isempty(std_thresh), std_thresh = 4; end

    % Flatten arrays
    x = x_data(:);
    y = y_data(:);

    % Remove NaN/Inf
    validIdx = ~isnan(x) & ~isnan(y) & ~isinf(x) & ~isinf(y);
    x_valid = x(validIdx);
    y_valid = y(validIdx);

    % Outlier exclusion (4 std)
    mu_x = mean(x_valid, 'omitnan'); std_x = std(x_valid, 'omitnan');
    mu_y = mean(y_valid, 'omitnan'); std_y = std(y_valid, 'omitnan');

    inlierIdx = (abs(x_valid - mu_x) <= std_thresh * std_x) & ...
                (abs(y_valid - mu_y) <= std_thresh * std_y);

    x_valid = x_valid(inlierIdx);
    y_valid = y_valid(inlierIdx);

    % Discretize into bins
    x_edges = linspace(x_min, x_max, num_bins + 1);
    bin_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;

    [bin_idx, ~] = discretize(x_valid, x_edges);

    y_mean = nan(1, num_bins);
    y_std  = nan(1, num_bins);

    for k = 1:num_bins
        y_in_bin = y_valid(bin_idx == k);
        if ~isempty(y_in_bin)
            y_mean(k) = mean(y_in_bin, 'omitnan');
            y_std(k)  = std(y_in_bin, 'omitnan');
        end
    end

    % Remove empty bins for fitting
    valid_bins = ~isnan(y_mean);
    bx = bin_centers(valid_bins);
    by = y_mean(valid_bins);
    bstd = y_std(valid_bins);

    % Linear regression
    p = polyfit(bx, by, 1);
    slope = p(1);
    x_fit = linspace(x_min*5, x_max*5, 100);
    y_fit = polyval(p, x_fit);
end