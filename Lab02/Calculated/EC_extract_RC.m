% Mystery load from three sweeps at the tee: load, open, short.
% V = (Vpp/2)*exp(j*phase). Divide by the open run so generator, cable and
% sync skew cancel. Subtract the short run to remove the grabber leads.
clear; clc; close all;

Z0 = 50;
M  = readmatrix('Data Collection - Mystery Load', 'NumHeaderLines', 1);

load_data  = M(:, 1:3);   load_data  = load_data (all(~isnan(load_data ), 2), :);
open_data  = M(:, 5:7);   open_data  = open_data (all(~isnan(open_data ), 2), :);
short_data = M(:, 9:11);  short_data = short_data(all(~isnan(short_data), 2), :);

freq_MHz = load_data(:,1);
w        = 2*pi*freq_MHz*1e6;

V_load  = (load_data(:,2)/2)  .* exp(1j*deg2rad(load_data(:,3)));
V_open  = (open_data(:,2)/2)  .* exp(1j*deg2rad(open_data(:,3)));
V_short = (short_data(:,2)/2) .* exp(1j*deg2rad(short_data(:,3)));

% Open and short were swept on a coarser list; put them on the load's list
V_open  = interp1(open_data(:,1),  V_open,  freq_MHz, 'linear', 'extrap');
V_short = interp1(short_data(:,1), V_short, freq_MHz, 'linear', 'extrap');

% 50-ohm source divider: V/V_open = Z/(Z0+Z)  ->  Z = Z0*r/(1-r)
r_load  = V_load  ./ V_open;   
Z_tee   = Z0 * r_load  ./ (1 - r_load);

r_short = V_short ./ V_open;   
Z_leads = Z0 * r_short ./ (1 - r_short);

Z_load    = Z_tee - Z_leads;
Gamma_load = (Z_load - Z0) ./ (Z_load + Z0);
L_leads   = imag(Z_leads) ./ w;

% Fit Z = (R || C) + L_res over 1-40 MHz (short run is not a clean L above that)
band  = freq_MHz <= 40;

% p = [R, C, L] from the load
model = @(p, w) 1 ./ (1/p(1) + 1j*w*p(2)*1e-12) + 1j*w*p(3)*1e-9;   % [R ohm, C pF, L nH]
err   = @(p) norm((model(p, w(band)) - Z_load(band)) ./ abs(Z_load(band)));
p     = fminsearch(err, [110 40 150]);
R = p(1);  C_pF = p(2);  L_res_nH = p(3);

fprintf('Grabber leads  L = %.0f nH\n', mean(L_leads(freq_MHz <= 35))*1e9);
fprintf('Load fit       R = %.1f ohm,  C = %.1f pF,  L_res = %.0f nH\n', R, C_pF, L_res_nH);
fprintf('Low-freq check Z_load(1 MHz) = %.1f %+.1fj ohm\n', real(Z_load(1)), imag(Z_load(1)));

Results = table(freq_MHz, real(Z_load), imag(Z_load), abs(Gamma_load), rad2deg(angle(Gamma_load)), ...
    'VariableNames', {'Freq_MHz', 'Re_Z', 'Im_Z', 'Gamma_mag', 'Gamma_deg'});
disp(Results)

fm = linspace(1, 80, 400)';  Zm = model(p, 2*pi*fm*1e6);
figure
plot(freq_MHz, real(Z_load), 'o', freq_MHz, imag(Z_load), 's', fm, real(Zm), '-', fm, imag(Zm), '--')
grid on; xlim([0 50]); xlabel('Freq (MHz)'); ylabel('\Omega')
legend('Re Z meas', 'Im Z meas', 'Re Z fit', 'Im Z fit')
title(sprintf('Z_{load}: %.0f \\Omega || %.0f pF, + %.0f nH leads', R, C_pF, L_res_nH))

figure
subplot(2,1,1); plot(freq_MHz, abs(Gamma_load), '-o'); grid on; ylabel('|\Gamma|'); xlim([0 50])
subplot(2,1,2); plot(freq_MHz, rad2deg(angle(Gamma_load)), '-o'); grid on; ylabel('\angle\Gamma (deg)'); xlabel('Freq (MHz)'); xlim([0 50])
