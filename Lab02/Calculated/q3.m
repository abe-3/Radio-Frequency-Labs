% Constants
t = linspace(0, 100, 5000); % Time in nanoseconds
td_BNC12 = 18.288;          % Time Delays of going through once
td_BNC18 = 27.432;
Z_shunt = 0;
ZO = 50;
ZL = inf;

if (Z_shunt == 0), Z_j = 0; else, Z_j = (Z_shunt * ZO) / (Z_shunt + ZO); end
gamma = (Z_j - ZO) / (Z_j + ZO);
T = 1 + gamma;

% Voltage Step
v = @(t, amp1, td1, amp2, td2) amp1*(t >= 2 + td1) + amp2*(t >= 2 + td2);

% Plotting
subplot(3,1,1)
plot(t, v(t, 0.5, 0, gamma*0.5, 2*td_BNC12), lineWidth=2)
title("Channel 1")
xlabel("Time [ns]")
ylabel("Voltage [V]")

subplot(3,1,2)
plot(t, v(t, T*0.5, td_BNC12, 0, 0), lineWidth=2)
title("Channel 2")
xlabel("Time [ns]")
ylabel("Voltage [V]")

subplot(3,1,3)
plot(t, v(t, T*0.5, td_BNC12+td_BNC18, 0, 0), lineWidth=2)
title("Channel 3")
xlabel("Time [ns]")
ylabel("Voltage [V]")
sgtitle(Z_shunt + "\Omega Shunt, " + ZL + "\Omega Termination", 'FontWeight', "Bold")