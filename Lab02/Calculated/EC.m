%% EC.m  --  Extract an unknown parallel R||C load from a tee-voltage sweep
%  Model: 50-ohm source -> cables -> tee (C_O || R_O) -> grabber leads (L_G) -> (R || C)
%  Matches EC_model_matched.asc, refit to the "Mystery Load" sweep.
clear; clc; close all;

%% ---------------- Settings ----------------
CsvFile = 'EC Sweep - Sheet1.csv';

Z0     = 50;        % generator Rser and cable impedance [ohm]
V_inc  = 1;         % incident-wave AMPLITUDE at the tee [V]
                    % (4 Vpp open-circuit -> 2 Vpp incident -> 1 V amplitude)
C_O    = 13e-12;    % Scope_Tee_C
R_O    = 1e6;       % Scope_R
L_G    = 290e-9;    % Grabber_Leads (refit; was 345n)

GenSyncSkew  = 2.24e-9;   % refit (was 2.495n)
BNC_6in      = 0.762e-9;
BNC_2ft_Sync = 3.048e-9;
DeltaDelay   = GenSyncSkew + BNC_6in - BNC_2ft_Sync;  % tee-path minus ref-path delay [s]

bCorrectCableSkew = true;

% Recorded-phase convention:
%   +1 -> Phase = ph(V_tee) - ph(V_ref)   <-- this sweep
%   -1 -> Phase = ph(V_ref) - ph(V_tee)   <-- LTSpice comment / original script
PhaseSign = +1;

bFitFixture = true;      % global fit also refits L_G and DeltaDelay
FitBand_MHz = [1 70];    % points used by the global fit (80 MHz is outside the model's reach)
AvgBand_MHz = [14 35];   % points averaged in the point-by-point extraction

%% ---------------- Load data ----------------
M = readmatrix(CsvFile, 'NumHeaderLines', 2);
M = M(all(~isnan(M(:,1:3)), 2), 1:3);
f_MHz    = M(:,1);
Vpp      = M(:,2);
PhaseDeg = M(:,3);

f = f_MHz*1e6;
w = 2*pi*f;

% Measured (1+Gamma) as recorded, i.e. still including the path skew:
%   a_raw = (1+Gamma) * exp(-j*w*DeltaDelay)
a_raw = (Vpp/2)/V_inc .* exp(1j*PhaseSign*deg2rad(PhaseDeg));

%% ---------------- Global least-squares fit of the full model ----------------
inFit = f_MHz >= FitBand_MHz(1) & f_MHz <= FitBand_MHz(2);
if bFitFixture
    p0 = [log(130), log(40e-12), log(L_G), DeltaDelay*1e9];
else
    p0 = [log(130), log(40e-12)];
end
cost = @(p) fitCost(p, w(inFit), a_raw(inFit), bFitFixture, L_G, DeltaDelay, C_O, R_O, Z0);
opts = optimset('TolX',1e-10, 'TolFun',1e-12, 'MaxFunEvals',2e4, 'MaxIter',2e4);
p = fminsearch(cost, p0, opts);
p = fminsearch(cost, p,  opts);   % restart once to polish
[R_fit, C_fit, L_fit, dT_fit] = unpackParams(p, bFitFixture, L_G, DeltaDelay);
fprintf('Fit cost: start %.4g -> final %.4g\n', cost(p0), cost(p));

% Sanity check on the phase convention: refit with the opposite sign
costFlip = @(q) fitCost(q, w(inFit), conj(a_raw(inFit)), bFitFixture, L_G, DeltaDelay, C_O, R_O, Z0);
pFlip = fminsearch(costFlip, fminsearch(costFlip, p0, opts), opts);
if costFlip(pFlip) < 0.5*cost(p)
    warning(['Opposite phase sign fits much better (cost %.4g vs %.4g). ' ...
             'Check PhaseSign vs. whether the CSV phase column was already negated.'], ...
             costFlip(pFlip), cost(p));
end

fprintf('Global fit (%g-%g MHz):\n', FitBand_MHz);
fprintf('  Unknown_R     = %.1f ohm\n', R_fit);
fprintf('  Unknown_C     = %.1f pF\n',  C_fit*1e12);
fprintf('  Grabber_Leads = %.0f nH\n',  L_fit*1e9);
fprintf('  Gen_sync_skew = %.3f ns  (net DeltaDelay = %.3f ns)\n\n', ...
        (dT_fit - BNC_6in + BNC_2ft_Sync)*1e9, dT_fit*1e9);

%% ---------------- Point-by-point extraction (original method) ----------------
L_use  = L_fit;      % same fixture as the global fit, so both methods share one model
dT_use = dT_fit;

if bCorrectCableSkew
    a = a_raw .* exp(1j*w*dT_use);
else
    a = a_raw;
end
Gamma = a - 1;
Z_T   = Z0*(1+Gamma)./(1-Gamma);        % impedance seen at the tee
Y_PT  = 1./Z_T - (1j*w*C_O + 1/R_O);    % remove scope + tee
Z_L   = 1./Y_PT - 1j*w*L_use;           % remove grabber leads
Y_L   = 1./Z_L;

R = 1./real(Y_L);
C = imag(Y_L)./w;

% Impedance and reflection coefficient of the unknown load alone (referenced to Z0)
%   measured: Z_L de-embedded above (tee, scope and leads removed)
%   fitted:   smooth R||C from the global fit
Gamma_L     = (Z_L - Z0)./(Z_L + Z0);
Z_L_fit     = R_fit ./ (1 + 1j*w*R_fit*C_fit);
Gamma_L_fit = (Z_L_fit - Z0)./(Z_L_fit + Z0);

Results = table(f_MHz, R, C*1e12, ...
    real(Z_L), imag(Z_L), abs(Gamma_L), rad2deg(angle(Gamma_L)), ...
    real(Z_L_fit), imag(Z_L_fit), abs(Gamma_L_fit), rad2deg(angle(Gamma_L_fit)), ...
    'VariableNames', {'Freq_MHz','R_ohm','C_pF', ...
                      'ReZL_ohm','ImZL_ohm','GammaL_mag','GammaL_deg', ...
                      'ReZL_fit_ohm','ImZL_fit_ohm','GammaL_fit_mag','GammaL_fit_deg'});
disp(Results)
writetable(Results, 'EC_Results.csv');   % optional: save alongside the input

inAvg = f_MHz >= AvgBand_MHz(1) & f_MHz <= AvgBand_MHz(2);
R_Avg = mean(R(inAvg));
C_Avg = mean(C(inAvg));
fprintf('Point-by-point average (%g-%g MHz): R = %.1f ohm, C = %.1f pF\n', ...
        AvgBand_MHz, R_Avg, C_Avg*1e12);

%% ---------------- Plots ----------------
fm = linspace(1e6, 80e6, 800).';
wm = 2*pi*fm;
a_mod = teeModel(wm, R_fit, C_fit, L_fit, C_O, R_O, Z0) .* exp(-1j*wm*dT_fit);

figure('Name','Model vs data');
subplot(2,1,1);
plot(f_MHz, Vpp, 'o', fm/1e6, 2*V_inc*abs(a_mod), '-'); grid on;
ylabel('V_{pp} at tee (V)'); legend('Measured','Fitted model','Location','best');
subplot(2,1,2);
plot(f_MHz, PhaseDeg, 'o', fm/1e6, PhaseSign*rad2deg(angle(a_mod)), '-'); grid on;
xlabel('Freq (MHz)'); ylabel('Recorded phase (deg)');

figure('Name','Point-by-point extraction');
subplot(2,1,1); plot(f_MHz, R, '-o'); yline(R_fit, '--', 'global fit');
ylabel('R (\Omega)'); grid on;
subplot(2,1,2); plot(f_MHz, C*1e12, '-o'); yline(C_fit*1e12, '--', 'global fit');
ylabel('C (pF)'); xlabel('Freq (MHz)'); grid on;

ZLm = R_fit ./ (1 + 1j*wm*R_fit*C_fit);
GLm = (ZLm - Z0)./(ZLm + Z0);
figure('Name','Unknown load: Z_L and Gamma_L');
subplot(3,1,1);
plot(f_MHz, real(Z_L), 'o', f_MHz, imag(Z_L), 's', fm/1e6, real(ZLm), '-', fm/1e6, imag(ZLm), '--');
grid on; ylabel('Z_L (\Omega)'); legend('Re meas','Im meas','Re fit','Im fit','Location','best');
subplot(3,1,2);
plot(f_MHz, abs(Gamma_L), 'o', fm/1e6, abs(GLm), '-'); grid on; ylabel('|\Gamma_L|');
subplot(3,1,3);
plot(f_MHz, rad2deg(angle(Gamma_L)), 'o', fm/1e6, rad2deg(angle(GLm)), '-'); grid on;
ylabel('\angle\Gamma_L (deg)'); xlabel('Freq (MHz)');

%% ---------------- Local functions ----------------
function a = teeModel(w, R, C, L, C_O, R_O, Z0)
% (1 + Gamma) at the tee for the fixture + R||C load
    Z_L = R ./ (1 + 1j*w*R*C);
    Y_T = 1j*w*C_O + 1/R_O + 1 ./ (1j*w*L + Z_L);
    Z_T = 1 ./ Y_T;
    a   = 2*Z_T ./ (Z_T + Z0);
end

function [R, C, L, dT] = unpackParams(p, bFitFixture, L_G, DeltaDelay)
    R = exp(p(1));
    C = exp(p(2));
    if bFitFixture
        L  = exp(p(3));
        dT = p(4)*1e-9;
    else
        L  = L_G;
        dT = DeltaDelay;
    end
end

function J = fitCost(p, w, a_raw, bFitFixture, L_G, DeltaDelay, C_O, R_O, Z0)
    [R, C, L, dT] = unpackParams(p, bFitFixture, L_G, DeltaDelay);
    a_mod = teeModel(w, R, C, L, C_O, R_O, Z0) .* exp(-1j*w*dT);
    J = sum(abs(a_mod - a_raw).^2);
end