% --- Filter out simulations with downrange error > 100% ---
valid_idx = abs(s_errors_pct) <= 100;
filtered_lats = lats(valid_idx);
filtered_lons = lons(valid_idx);
filtered_errors_pct = s_errors_pct(valid_idx);

% === 1. Filtered Landing Dispersion Plot ===
figure;
hold on; grid on;

filtered_lons = lon_c * ones(size(filtered_lons));  % Fix longitudes for visualization

plot(filtered_lons, filtered_lats, 'g.', 'MarkerSize', 10, 'DisplayName', 'Successful Landings');

lat_mean = mean(filtered_lats);
plot(lon_c, lat_mean, 'ko', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Mean Landing');

plot(lon_c, lat_c, 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Nominal Target');

% Draw 50 km radius circle
r_nominal = 50000;
deg_per_m_lat = 1 / 111320;
deg_per_m_lon = 1 / (111320 * cosd(lat_c));
theta = linspace(0, 2*pi, 300);
lat_circle_nom = lat_c + r_nominal * deg_per_m_lat * sin(theta);
lon_circle_nom = lon_c + r_nominal * deg_per_m_lon * cos(theta);
plot(lon_circle_nom, lat_circle_nom, 'r--', 'LineWidth', 2, 'DisplayName', '50 km Radius');

title('Filtered Landing Dispersion (Downrange Error ≤ 100%)');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
legend('Location', 'bestoutside');
axis equal;

% === 2. Filtered Histogram ===
figure;
histogram(abs(filtered_errors_pct), 'BinWidth', 1, 'FaceColor', [0.1 0.6 0.2]);
xlabel('Downrange Error [%]');
ylabel('Frequency');
title('Filtered Histogram of Downrange Error');
grid on;

% === 3. Print Filtered Statistics ===
mean_filtered_err = mean(abs(filtered_errors_pct));
std_filtered_err = std(abs(filtered_errors_pct));
fprintf('Filtered Downrange Error: Mean = %.2f%% | Std = %.2f%%\n', mean_filtered_err, std_filtered_err);

% === 4. Count and Report Failure Cases ===
num_total = length(s_errors_pct);
num_failures = sum(abs(s_errors_pct) > 100);
num_success = num_total - num_failures;

fprintf('Total Simulations: %d\n', num_total);
fprintf('Successful Landings (≤ 100%% error): %d\n', num_success);
fprintf('Failures (> 100%% error): %d\n', num_failures);
fprintf('Failure Rate: %.2f%%\n', 100 * num_failures / num_total);

% === 5. Determine Most Frequent Percentile Bin ===
edges = 0:5:100;  % 5% bin width up to 100%
[counts, bin_edges] = histcounts(abs(filtered_errors_pct), edges);

[max_count, max_idx] = max(counts);
peak_bin_start = bin_edges(max_idx);
peak_bin_end = bin_edges(max_idx+1);

% Convert 5% of nominal downrange into meters
nominal_downrange = 1.75364e6;  % meters
error_5pct_meters = 0.05 * nominal_downrange;

fprintf('Most frequent downrange error bin: %.1f%% – %.1f%% (%d cases)\n', ...
    peak_bin_start, peak_bin_end, max_count);
fprintf('A 5%% downrange error corresponds to approximately %.0f meters.\n', ...
    error_5pct_meters);
