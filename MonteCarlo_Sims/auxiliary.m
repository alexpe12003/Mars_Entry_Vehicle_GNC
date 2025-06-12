% Organize data from results structure
s_errors_pct = abs([results.s_error_pct]);
s_finals = [results.s_final];
lats = [results.final_lat];
lons = [results.final_lon];

% === 1. Best, Worst, and Average Cases ===
% Best Case: minimum error
best_case = min(s_errors_pct);
best_case_index = find(s_errors_pct == best_case);
best_case_s_final = s_finals(best_case_index);

% Worst Success Case: maximum error (below 100%)
worst_success_case = max(s_errors_pct(s_errors_pct < 100));
worst_success_case_index = find(s_errors_pct == worst_success_case);
worst_success_case_s_final = s_finals(worst_success_case_index);

% Worst Failure Case: maximum error (above 100%)
worst_failure_case = max(s_errors_pct(s_errors_pct > 100));
worst_failure_case_index = find(s_errors_pct == worst_failure_case);
worst_failure_case_s_final = s_finals(worst_failure_case_index);

% === 2. Successful and Failure Cases Statistics ===
% Successful cases (error < 100%)
successful_errors = s_errors_pct(s_errors_pct < 100);
mean_success_err = mean(successful_errors);
std_success_err = std(successful_errors);

% Failure cases (error > 100%)
failure_errors = s_errors_pct(s_errors_pct > 100);
mean_failure_err = mean(failure_errors);
std_failure_err = std(failure_errors);

% === 3. Count the number of successes and failures ===
num_successes = length(successful_errors);
num_failures = length(failure_errors);

% === 4. Print Key Statistics ===
fprintf('Best Case: %.2f%% error | Final Landing Distance = %.2f m\n', best_case, best_case_s_final);
fprintf('Worst Success Case: %.2f%% error | Final Landing Distance = %.2f m\n', worst_success_case, worst_success_case_s_final);
fprintf('Worst Failure Case: %.2f%% error | Final Landing Distance = %.2f m\n', worst_failure_case, worst_failure_case_s_final);

% Print the number of successes and failures
fprintf('\nNumber of Successful Landings: %d\n', num_successes);
fprintf('Number of Failed Landings: %d\n', num_failures);

fprintf('\nSuccessful Landing Stats:\n');
fprintf('  Mean Error = %.2f%% | Std Dev = %.2f%%\n', mean_success_err, std_success_err);

fprintf('\nFailure Landing Stats:\n');
fprintf('  Mean Error = %.2f%% | Std Dev = %.2f%%\n', mean_failure_err, std_failure_err);

% === 5. Additional Summary Information ===
% Mean and Std Dev of All Errors
mean_error = mean(abs(s_errors_pct));
std_error = std(abs(s_errors_pct));
fprintf('\nOverall Error Stats:\n');
fprintf('  Mean Error = %.2f%% | Std Dev = %.2f%%\n', mean_error, std_error);

