file_name = '2a';
exp_data = readmatrix(file_name + ".csv");
sim_data = readmatrix(file_name + ".txt");

exp_time = exp_data(:, 1);
sim_time = sim_data(:, 1);

for idx = 1:3
    subplot(3, 1, idx)
    hold on
    plot(exp_time, exp_data(:, idx+1), lineWidth=2);
    plot(sim_time, sim_data(:, idx+1), lineWidth=2, lineStyle='--');
    hold off
    xlabel("Time [Sec]")
    ylabel("Voltage [V]")
    title("Voltage at Channel # " + idx)
    legend("Measured Data", "Experimental Data")
end
