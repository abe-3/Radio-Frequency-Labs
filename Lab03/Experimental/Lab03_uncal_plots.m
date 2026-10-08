%% Lab03_uncal_plots.m
% One figure per measurement type (open, load50, short), each with
% magnitude, phase and Smith chart tiles. Expects files such as
% uncal_open_mag.csv, uncal_open_phase.csv, uncal_open_smith.csv
% somewhere under the current folder.

dataType = "cal";

Z0       = 50;     % reference impedance (ohms)
gammaMax = 1.5;    % drop Smith points with |Gamma| above this (uncal artifacts)

types = ["open","load50","short"];
kinds = ["mag","phase","smith"];
titles = ["Open Termination", "50 \Omega Termination", "Short Termination"];

% Filename pattern for each type (allows load50 / load_50 / load-50)
typePat = struct("open","open", "load50","load[_\- ]?50", "short","short");

% Collect all CSVs (current folder and subfolders) once
files = dir(fullfile(pwd, "**", "*.csv"));
files = files(~[files.isdir]);
names = lower(string({files.name}));

for k = 1:numel(types)
    type = types(k);
    data = struct();
    plotTitle = titles(k);

    for j = 1:numel(kinds)
        kind = kinds(j);

        % Each keyword is checked independently, so order doesn't matter
        hasType  = ~cellfun("isempty", regexp(cellstr(names), typePat.(type), "once"));
        hasUncal = contains(names, dataType);
        hasKind  = contains(names, kind);
        idx = find(hasType & hasUncal & hasKind);

        if isempty(idx)
            fprintf("CSV files found:\n");
            fprintf("  %s\n", fullfile(string({files.folder}), string({files.name})));
            error("No file found for type '%s', kind '%s'.", type, kind);
        elseif numel(idx) > 1
            warning("Multiple matches for %s/%s; using %s", ...
                type, kind, files(idx(1)).name);
        end

        f = files(idx(1));
        data.(kind) = readmatrix(fullfile(f.folder, f.name));
    end

    % Magnitude and phase: [freq, value]
    mag   = data.mag(all(isfinite(data.mag(:,1:2)),2), :);
    phase = data.phase(all(isfinite(data.phase(:,1:2)),2), :);

    % Smith file: [freq, R, X, ...] with R+jX in ohms -> convert to Gamma
    sm    = data.smith(all(isfinite(data.smith(:,1:3)),2), :);
    Z     = sm(:,2) + 1i*sm(:,3);
    Gamma = (Z - Z0) ./ (Z + Z0);
    % keep  = isfinite(Gamma) & abs(Gamma) <= gammaMax;
    % if any(~keep)
        % fprintf("%s: dropped %d of %d Smith points with |Gamma| > %.1f\n", ...
            % type, nnz(~keep), numel(keep), gammaMax);
    % end
    % Gamma = Gamma(keep);

    figure("Name", type);
    tiledlayout(1,3);

    nexttile;
    plot(mag(:,1)/1e6, mag(:,2), "LineWidth", 1.2, LineWidth=2);
    grid on; 
    xlabel("Frequency (MHz)"); 
    ylabel("Magnitude");
    title(plotTitle + " - Magnitude");

    nexttile;
    plot(phase(:,1)/1e6, phase(:,2), "LineWidth", 1.2, LineWidth=2);
    grid on; 
    xlabel("Frequency (MHz)"); 
    ylabel("Phase");
    title(plotTitle + " - Phase");

    nexttile;
    if exist("smithplot", "file")
        smithplot(Gamma);
    else
        % Fallback without RF Toolbox: Gamma plane with unit circle
        th = linspace(0, 2*pi, 400);
        plot(cos(th), sin(th), "k-"); hold on;
        plot(real(Gamma), imag(Gamma), "LineWidth", 1.2);
        axis equal; grid on;
        xlabel("Re(\Gamma)"); ylabel("Im(\Gamma)");
    end
    title(plotTitle + " - Smith Chart");
end
