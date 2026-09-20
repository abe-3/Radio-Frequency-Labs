% Constants
t = linspace(0, 100, 5000); % Time in Nano Seconds
td_BNC12 = 18.288;          % Time Delays of going through once
td_BNC18 = 27.432;
ZL = 200;
ZO = 50;
if (ZL == inf), gamma=1; else, gamma = (ZL - ZO) / (ZL + ZO); end

% Voltage Step
v = @(t, td_start, td_end) 0.5*(t >= 2 + td_start) + gamma*0.5*(t >= 2 + td_end);

% Plotting
subplot(3,1,1)
plot(t, v(t, 2, 2*(td_BNC12 + td_BNC18)), lineWidth=2)
title("Channel 1")
xlabel("Time [ns]")
ylabel("Voltage [V]")

subplot(3,1,2)
plot(t, v(t, td_BNC12, td_BNC12+2*td_BNC18), lineWidth=2)
title("Channel 2")
xlabel("Time [ns]")
ylabel("Voltage [V]")

subplot(3,1,3)
plot(t, v(t, td_BNC12+td_BNC18, td_BNC12+td_BNC18), lineWidth=2)
title("Channel 3")
xlabel("Time [ns]")
ylabel("Voltage [V]")
sgtitle("Channel Voltages with Load of " + ZL + "\Omega", 'FontWeight', "Bold")