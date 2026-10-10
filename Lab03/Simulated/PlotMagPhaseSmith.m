function Output = PlotMagPhaseSmith(bUseCorrectionFactor)

    arguments

        bUseCorrectionFactor logical = false;

    end

    Files = dir(fullfile('SimulationData', '*.txt'));
    Files = Files(~[Files.isdir]);

    figure(1);
    tiledlayout(2,4);

    Tile1 = nexttile(1, [1 2]); 
    hold(Tile1,'on'); 
    grid(Tile1,'on');
    title(Tile1,'Reflection Coefficient Magnitude'); 
    ylabel(Tile1,'|\Gamma|');
    xlabel(Tile1,'Frequency (MHz)');

    Tile2 = nexttile(5, [1 2]); 
    hold(Tile2,'on'); 
    grid(Tile2,'on');
    title(Tile2,'Reflection Coefficient Phase'); 
    ylabel(Tile2,'Phase (degrees)');
    xlabel(Tile2,'Frequency (MHz)');

    Colors = colororder('gem');

    for FileIndex = 1:numel(Files)
        Data = readmatrix(fullfile(Files(FileIndex).folder, Files(FileIndex).name), ...
        'FileType','text', 'Delimiter',{'\t', ','});
        
        FoundLengths = regexp(Files(FileIndex).name, '(?<=_)[0-9]+(?=m_)', 'match'); % Look for underscore followed by any number of digits, followed by another underscore
        CableLength = str2double(FoundLengths{1});

        Frequency = Data(:,1);

        Frequency = Frequency / 1e6;

        Gamma = complex(Data(:,2), Data(:,3));
    
        bInFreqRange = Frequency >= 70 & Frequency <= 130;

        Frequency = Frequency(bInFreqRange);
        CorrectionFactor = exp(2j*2*pi*Frequency*CableLength/(0.8*300));

        Gamma = Gamma(bInFreqRange);

        if (bUseCorrectionFactor)
            Gamma = Gamma .* CorrectionFactor;
        end

        DataName = erase(Files(FileIndex).name, {'.txt', '2a_'});
        DataName = replace(DataName, "_", " ");
        DataName = replace(DataName, "open", "Open Termination");
        DataName = replace(DataName, "short", "Short Termination");
        DataName = replace(DataName, "50", "50 Ohm Termination");

        plot(Tile1, Frequency, abs(Gamma), 'DisplayName', DataName, LineWidth=2);
        plot(Tile2, Frequency, wrapToPi(angle(Gamma))*180/pi, 'DisplayName', DataName, LineWidth=2);

        CurrentTile = nexttile;
        title(DataName);
        SmithPlot = smithplot(CurrentTile, Gamma, 'Color', Colors(FileIndex, :), 'MarkerSize', 15, 'Marker', '.');
        SmithPlot.TitleTop = DataName;
        SmithPlot.TitleTopFontSizeMultiplier = 2.5;
    end

    legend(Tile1, 'Location', 'best');
    legend(Tile2, 'Location', 'best');
end