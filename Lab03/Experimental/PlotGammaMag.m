Z0 = 50;

data = struct();
data.("smith") = readmatrix(fullfile("lmatch_add_180n.csv"));

sm    = data.smith(all(isfinite(data.smith(:,1:3)),2), :);
freq = sm(:,1);
Z     = sm(:,2) + 1i*sm(:,3);
Gamma = (Z - Z0) ./ (Z + Z0);

GammaMag = abs(Gamma);

markerIndex = find((GammaMag == min(GammaMag)));

figure()
plot(freq, GammaMag, '-o', 'MarkerIndices', markerIndex, 'MarkerEdgeColor', 'black', 'MarkerSize', 8, LineWidth=2)
xlabel("Frequency (MHz)")
ylabel("|\Gamma|")
title("Reflection Coefficient vs. Frequency")