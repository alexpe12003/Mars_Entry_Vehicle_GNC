% === 5. Plot Landing Dispersion (Longitude Fixed) ===
figure;
hold on; grid on;

% Force all longitudes to the nominal value for plotting
lons = lon_c * ones(size(lons));  % Override longitudes

% Plot Monte Carlo landing points
plot(lons, lats, 'b.', 'MarkerSize', 10, 'DisplayName', 'Landing Points');

% Plot centroid (mean landing)
lat_mean = mean(lats);
lon_mean = lon_c;
plot(lon_mean, lat_mean, 'ko', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Mean Landing');

% Plot nominal target point
plot(lon_c, lat_c, 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Nominal Target');

% === Draw 100 km radius circle around nominal point ===
r_nominal = 100000;  % 10 km in meters
deg_per_m_lat = 1 / 111320;
deg_per_m_lon = 1 / (111320 * cosd(lat_c));

theta = linspace(0, 2*pi, 300);
lat_circle_nom = lat_c + r_nominal * deg_per_m_lat * sin(theta);
lon_circle_nom = lon_c + r_nominal * deg_per_m_lon * cos(theta);
plot(lon_circle_nom, lat_circle_nom, 'r--', 'LineWidth', 2, 'DisplayName', '100 km Nominal Circle');

% === Final plot styling ===
title('Landing Dispersion with Bank Angle Modulation');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
legend('Location', 'bestoutside');
axis equal;


% === 1. Carregar resultados guardados ===
% Espera-se que contenha: results(i).s_final para i = 1:N

% === 2. Definir s_star (nominal downrange) ===
s_star_nominal = 2.209934992639319e+06;  % metros — substitui com o teu valor nominal real

% === 3. Calcular erros percentuais ===
N = length(results);
s_errors_pct = zeros(1, N);

for i = 1:N
    s_final = results(i).s_final;
    s_errors_pct(i) = 100 * (s_final - s_star_nominal) / s_star_nominal;
end

% === 4. Plotar histograma ===
figure;
histogram(s_errors_pct, 'BinWidth', 1, 'FaceColor', [0.2 0.4 0.8]);
xlabel('Downrange Error [%]');
ylabel('Frequency');
title('Histogram of Downrange Error (% of Nominal)');
grid on;

% === 5. Estatísticas ===
mean_err = mean(s_errors_pct);
std_err = std(s_errors_pct);
fprintf('Downrange Error: Mean = %.2f%% | Std = %.2f%%\n', mean_err, std_err);


