function J = simulate_qc_cost(sigma_vector)

    t_total = 500;                         % Total sim time
    t_ref = linspace(0, t_total, length(sigma_vector));
    sigma_ref_ga = timeseries(sigma_vector, t_ref);

    ts_sigma.Name = 'sigma_ref_ga';

    % Wrap it into a Dataset
    signalDataset = Simulink.SimulationData.Dataset;
    signalDataset = signalDataset.addElement(ts_sigma);

    save('sigma_ref_dataset_ga.mat', 'signalDataset');

    simOut = sim('entryVehicle', 'StopTime', '700', ...
        'SrcWorkspace','base');

    % Assume you log q_c as a signal
    qc = simOut.q_c.signals.values;
    impulse= simOut.impulse.signals.values;
    t = simOut.q_c.time;

    J = trapz(t, qc) + trapz(t,impulse);    % Total heat flux
end

