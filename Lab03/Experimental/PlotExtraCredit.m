Z0 = 50;
data = struct();
data.("smith") = readmatrix(fullfile("mystery_smith.csv"));

% Filter valid columns
sm = data.smith(all(isfinite(data.smith(:,1:3)),2), :);

% Keep frequency in Hz for calculations, create MHz version for plotting
freq_Hz = sm(:,1);          
freq_MHz = sm(:,1) / 1e6;   

% Impedance and Admittance calculations
Z = sm(:,2) + 1i*sm(:,3);
Y = 1 ./ Z;                 % Admittance (Y = G + jB)
Gamma = (Z - Z0) ./ (Z + Z0);

conductance = real(Y);      % G (Real part of admittance)
susceptance = imag(Y);      % B (Imaginary part of admittance)

% Extracting parallel component values based on susceptance sign:
% We have positive susceptance so the reactive part is probably a capacitor
capacitance = susceptance ./ (2 * pi * freq_Hz); 
inductance  = 1 ./ (2 * pi * freq_Hz .* abs(susceptance));

%% --- Plotting ---
figure(1)
% Plot capacitance where susceptance is capacitive (positive), converted to pF
plot(freq_MHz, capacitance * 1e12, 'LineWidth', 2); 
xlabel("Frequency (MHz)")
ylabel("Parallel Capacitance (pF)")
grid on;

figure(2)
tiledlayout(2, 2)

Tile1 = nexttile(1);
plot(freq_MHz, conductance, 'LineWidth', 2);
xlabel("Frequency (MHz)")
ylabel("Conductance (S)")
grid on;

Tile2 = nexttile(3);
plot(freq_MHz, susceptance, 'LineWidth', 2);
xlabel("Frequency (MHz)")
ylabel("Susceptance (S)")
grid on;

Tile3 = nexttile(2, [2, 1]);
% Use 'ZY' GridType to display admittance lines for parallel components on
% the Smith chart
smithplot(Gamma, 'GridType', 'ZY');

% Based on the plots, susceptance seems reasonable (steady relative to the
% rest of the data and not subject to much noise) between ~10 and ~40MHz

keep = find((10e+6 < freq_Hz & freq_Hz < 40e+6));

fprintf("Median capacitance in acceptable range is %d\n", median(capacitance(keep)))
fprintf("Mean capacitance in acceptable range is %d\n", mean(capacitance(keep)))

