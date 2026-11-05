%% Laminar parallel-plate channel flow: analytical solution and CFD post-processing
% Run this file from the folder containing the exported Fluent profile files.
% The script evaluates the parallel-plate solution, compares the six outlet
% profiles with that solution, and saves the figures and summary tables.

clear;
close all;
clc;

%% Figure appearance
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 14);
set(groot, 'defaultTextFontSize', 14);
set(groot, 'defaultAxesLineWidth', 1.0);
set(groot, 'defaultLineLineWidth', 1.8);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'on');
set(groot, 'defaultAxesXMinorTick', 'on');
set(groot, 'defaultAxesYMinorTick', 'on');
set(groot, 'defaultLegendFontName', 'Times New Roman');
set(groot, 'defaultLegendFontSize', 11);

figureFolder = fullfile(pwd, 'figures_matlab');
if ~exist(figureFolder, 'dir')
    mkdir(figureFolder);
end

meshColours = [0.0000 0.4470 0.7410; ...
               0.8500 0.3250 0.0980; ...
               0.4660 0.6740 0.1880];
black = [0.10 0.10 0.10];
markers = {'o', 's', 'd'};

%% Fluid properties and channel dimensions
H = 0.02;                 % Full channel height [m]
L = 0.10;                 % Channel length [m]
Re = 12;
rho = 1000;               % Density [kg/m^3]
mu = 0.00102;              % Dynamic viscosity [Pa s]
alpha = 0.143e-6;          % Thermal diffusivity [m^2/s]
cp = 4185;                 % Specific heat capacity [J/(kg K)]
T_inlet = 300;             % Thermal-case inlet temperature [K]
T_lowerWall = 50;          % Lower-wall temperature [K]
T_upperWall = 350;         % Upper-wall temperature [K]

nu = mu / rho;
Dh = 2 * H;
Umean = Re * nu / Dh;

% The analytical formulation uses y = 0 at the channel centreline and walls at y = +/-H/2.
% From the analytical mean-velocity relation:
% Umean = DeltaP*H^2/(12*mu*L).
% This gives the pressure drop for the ideal fully developed solution.
DeltaP_ideal = 12 * mu * L * Umean / H^2;

% Analytical fully developed velocity profile:
% u(y) = DeltaP/(2*mu*L) * (H^2/4 - y^2).
yPlot = linspace(-H/2, H/2, 500)';
Yplot = yPlot + H/2;      % Convert to the lower-wall coordinate used by Fluent.
Uplot = (DeltaP_ideal / (2 * mu * L)) .* (H^2/4 - yPlot.^2);

% Fluent outlet files use Y = 0 at the lower wall, so convert back to y.
Utheory = @(Y) (DeltaP_ideal / (2 * mu * L)) .* ...
                 (H^2/4 - (Y - H/2).^2);

uMax = DeltaP_ideal * H^2 / (8 * mu * L);
kThermal = alpha * rho * cp;

fprintf('\nANALYTICAL REFERENCE VALUES\n');
fprintf('Reynolds number                 : %.4f\n', Re);
fprintf('Hydraulic diameter              : %.6f m\n', Dh);
fprintf('Kinematic viscosity             : %.6e m^2/s\n', nu);
fprintf('Mean velocity                   : %.6e m/s\n', Umean);
fprintf('Maximum velocity                : %.6e m/s\n', uMax);
fprintf('Ideal fully developed pressure drop: %.6e Pa\n', DeltaP_ideal);
fprintf('Thermal conductivity             : %.6f W/(m K)\n\n', kThermal);

% Check the wall conditions and the centreline maximum.
if abs(Utheory(0)) > 1e-12 || abs(Utheory(H)) > 1e-12
    error('The analytical profile does not satisfy the no-slip wall conditions.');
end
if abs(max(Uplot) - 1.5 * Umean) > 1e-10
    error('The analytical maximum does not equal 1.5 times the mean velocity.');
end

% Keep the analytical reference values with the processed results.
Parameter = {'Reynolds number'; 'Channel height'; 'Channel length'; ...
    'Density'; 'Dynamic viscosity'; 'Kinematic viscosity'; ...
    'Hydraulic diameter'; 'Mean velocity'; 'Maximum velocity'; ...
    'Ideal pressure drop'; 'Thermal diffusivity'; ...
    'Specific heat capacity'; 'Thermal conductivity'; ...
    'Thermal inlet temperature'; 'Lower-wall temperature'; ...
    'Upper-wall temperature'};
Value = [Re; H; L; rho; mu; nu; Dh; Umean; uMax; ...
    DeltaP_ideal; alpha; cp; kThermal; T_inlet; T_lowerWall; T_upperWall];
Unit = {'-'; 'm'; 'm'; 'kg/m^3'; 'Pa s'; 'm^2/s'; 'm'; ...
    'm/s'; 'm/s'; 'Pa'; 'm^2/s'; 'J/(kg K)'; 'W/(m K)'; ...
    'K'; 'K'; 'K'};
AnalyticalValues = table(Parameter, Value, Unit);
writetable(AnalyticalValues, fullfile(pwd, 'Analytical_Reference_Values.csv'));

%% Files exported from Fluent
profileFiles = { ...
    'coarse_PO.txt'; ...
    'medium_PO.txt'; ...
    'fine_PO.txt'; ...
    'coarse_OF.txt'; ...
    'medium_OF.txt'; ...
    'fine_OF.txt'};

meshNames = {'Coarse'; 'Medium'; 'Fine'; 'Coarse'; 'Medium'; 'Fine'};
outletNames = {'Pressure Outlet'; 'Pressure Outlet'; 'Pressure Outlet'; ...
               'Outflow'; 'Outflow'; 'Outflow'};
fileTags = {'coarse_PO'; 'medium_PO'; 'fine_PO'; ...
            'coarse_OF'; 'medium_OF'; 'fine_OF'};
cellCounts = [8000; 32000; 128000; 8000; 32000; 128000];

firstOrderFile = 'coarse_PO_1st_order_scheme.txt';
thermalFile = 'thermal_outlet_temperature.txt';
pressureFile = 'centerline_pressure_medium_PO.txt';

%% Outlet velocity profiles and error calculations
nCases = numel(profileFiles);
Results = cell(nCases, 1);

nPoints = zeros(nCases, 1);
maxCFD = zeros(nCases, 1);
maxTheoryAtSamples = zeros(nCases, 1);
meanAbsError = zeros(nCases, 1);
maxAbsError = zeros(nCases, 1);
normalisedL2 = zeros(nCases, 1);

for iCase = 1:nCases
    [Y, Ucfd] = readTwoColumnFile(profileFiles{iCase});

    if any(Y < -1e-8 | Y > H + 1e-8)
        error('Coordinates in %s must lie between 0 and %.4f m.', ...
            profileFiles{iCase}, H);
    end

    Uana = Utheory(Y);
    absErr = abs(Ucfd - Uana);

    % Relative error is not evaluated where the analytical velocity is zero.
    relErr = nan(size(absErr));
    nonzero = abs(Uana) > 1e-12;
    relErr(nonzero) = 100 * absErr(nonzero) ./ abs(Uana(nonzero));

    Results{iCase} = struct('Y', Y, 'Ucfd', Ucfd, 'Uana', Uana, ...
        'absError', absErr, 'relError', relErr);

    nPoints(iCase) = numel(Y);
    maxCFD(iCase) = max(Ucfd);
    maxTheoryAtSamples(iCase) = max(Uana);
    meanAbsError(iCase) = mean(absErr);
    maxAbsError(iCase) = max(absErr);

    % Discrete relative L2 error, evaluated at the exported cell-centre points.
    normalisedL2(iCase) = 100 * norm(Ucfd - Uana) / norm(Uana);

    Pointwise = table(Y, Ucfd, Uana, absErr, relErr, ...
        'VariableNames', {'Y_m', 'CFD_XVelocity_m_s', ...
        'AnalyticalVelocity_m_s', 'AbsoluteError_m_s', ...
        'RelativeError_percent'});
    writetable(Pointwise, fullfile(pwd, ...
        ['pointwise_error_' fileTags{iCase} '.csv']));
end

Summary = table(meshNames, outletNames, cellCounts, nPoints, maxCFD, ...
    maxTheoryAtSamples, meanAbsError, maxAbsError, normalisedL2, ...
    'VariableNames', {'Mesh', 'OutletCondition', 'CellCount', ...
    'ProfilePoints', 'MaxCFD_m_s', 'MaxAnalyticalAtSamples_m_s', ...
    'MeanAbsoluteError_m_s', 'MaximumAbsoluteError_m_s', ...
    'NormalisedL2Error_percent'});

disp('HYDRODYNAMIC ERROR SUMMARY');
disp(Summary);
writetable(Summary, fullfile(pwd, 'CFD_Comparison_Summary.csv'));

%% Figure 1: analytical fully developed velocity profile
fig = figure('Name', 'Analytical velocity profile', ...
    'Position', [100 100 900 600]);
plot(Yplot, Uplot, '-', 'Color', black, 'LineWidth', 2.2);
formatAxes('Transverse coordinate, Y [m]', 'Axial velocity, u [m/s]');
xlim([0 H]);
exportFigure(fig, figureFolder, 'analytical_velocity_profile');

%% Figure 2: Pressure Outlet profiles
fig = figure('Name', 'Pressure Outlet velocity comparison', ...
    'Position', [100 100 900 600]);
hold on;
plot(Yplot, Uplot, '-', 'Color', black, 'LineWidth', 2.4, ...
    'DisplayName', 'Analytical solution');

for iCase = 1:3
    markerIndex = chooseMarkerIndices(numel(Results{iCase}.Y), 3, iCase);
    plot(Results{iCase}.Y, Results{iCase}.Ucfd, ...
        'Color', meshColours(iCase,:), 'LineStyle', '-', ...
        'Marker', markers{iCase}, 'MarkerIndices', markerIndex, ...
        'MarkerSize', 8.5, 'MarkerEdgeColor', black, ...
        'MarkerFaceColor', meshColours(iCase,:), ...
        'DisplayName', [meshNames{iCase} ' mesh']);
end

formatAxes('Transverse coordinate, Y [m]', 'Axial velocity, u [m/s]');
xlim([0 H]);
legend('Location', 'best', 'Box', 'off');
exportFigure(fig, figureFolder, 'velocity_profiles_pressure_outlet');

%% Figure 3: Outflow profiles
fig = figure('Name', 'Outflow velocity comparison', ...
    'Position', [100 100 900 600]);
hold on;
plot(Yplot, Uplot, '-', 'Color', black, 'LineWidth', 2.4, ...
    'DisplayName', 'Analytical solution');

for iCase = 4:6
    meshIndex = iCase - 3;
    markerIndex = chooseMarkerIndices(numel(Results{iCase}.Y), 3, meshIndex);
    plot(Results{iCase}.Y, Results{iCase}.Ucfd, ...
        'Color', meshColours(meshIndex,:), 'LineStyle', '-', ...
        'Marker', markers{meshIndex}, 'MarkerIndices', markerIndex, ...
        'MarkerSize', 8.5, 'MarkerEdgeColor', black, ...
        'MarkerFaceColor', meshColours(meshIndex,:), ...
        'DisplayName', [meshNames{iCase} ' mesh']);
end

formatAxes('Transverse coordinate, Y [m]', 'Axial velocity, u [m/s]');
xlim([0 H]);
legend('Location', 'best', 'Box', 'off');
exportFigure(fig, figureFolder, 'velocity_profiles_outflow');

%% Figure 4: first- and second-order schemes on the coarse mesh
if isfile(firstOrderFile)
    [Yfirst, Ufirst] = readTwoColumnFile(firstOrderFile);

    if any(Yfirst < -1e-8 | Yfirst > H + 1e-8)
        error('Coordinates in %s must lie between 0 and %.4f m.', ...
            firstOrderFile, H);
    end

    UanaFirst = Utheory(Yfirst);
    absFirst = abs(Ufirst - UanaFirst);
    l2First = 100 * norm(Ufirst - UanaFirst) / norm(UanaFirst);

    % The coarse Pressure Outlet profile is the second-order result.
    Ysecond = Results{1}.Y;
    Usecond = Results{1}.Ucfd;
    UanaSecond = Utheory(Ysecond);
    absSecond = abs(Usecond - UanaSecond);
    l2Second = 100 * norm(Usecond - UanaSecond) / norm(UanaSecond);

    fig = figure('Name', 'Coarse mesh scheme comparison', ...
        'Position', [100 100 900 600]);
    hold on;
    plot(Yplot, Uplot, '-', 'Color', black, 'LineWidth', 2.4, ...
        'DisplayName', 'Analytical solution');

    markerSecond = chooseMarkerIndices(numel(Ysecond), 3, 1);
    markerFirst = chooseMarkerIndices(numel(Yfirst), 3, 2);

    plot(Ysecond, Usecond, '-o', 'Color', meshColours(1,:), ...
        'MarkerIndices', markerSecond, 'MarkerSize', 8.5, ...
        'MarkerEdgeColor', black, 'MarkerFaceColor', meshColours(1,:), ...
        'DisplayName', 'Second-order upwind');

    plot(Yfirst, Ufirst, '--s', 'Color', meshColours(2,:), ...
        'MarkerIndices', markerFirst, 'MarkerSize', 8.5, ...
        'MarkerEdgeColor', black, 'MarkerFaceColor', meshColours(2,:), ...
        'DisplayName', 'First-order upwind');

    formatAxes('Transverse coordinate, Y [m]', 'Axial velocity, u [m/s]');
    xlim([0 H]);
    legend('Location', 'best', 'Box', 'off');
    exportFigure(fig, figureFolder, 'scheme_comparison_coarse');

    Scheme = {'First order'; 'Second order'};
    ProfilePoints = [numel(Yfirst); numel(Ysecond)];
    MeanAbsoluteError = [mean(absFirst); mean(absSecond)];
    MaximumAbsoluteError = [max(absFirst); max(absSecond)];
    NormalisedL2Error = [l2First; l2Second];

    SchemeSummary = table(Scheme, ProfilePoints, MeanAbsoluteError, ...
        MaximumAbsoluteError, NormalisedL2Error);
    disp('COARSE MESH SCHEME COMPARISON');
    disp(SchemeSummary);
    writetable(SchemeSummary, fullfile(pwd, 'Scheme_Comparison_Summary.csv'));
else
    warning('The first-order profile file %s was not found. Figure 4 was skipped.', ...
        firstOrderFile);
end

%% Figure 5: normalised L2 error against cell count
fig = figure('Name', 'Mesh error comparison', ...
    'Position', [100 100 900 600]);
hold on;

plot(cellCounts(1:3), normalisedL2(1:3), '-o', ...
    'Color', meshColours(1,:), 'MarkerSize', 9, ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', meshColours(1,:), ...
    'DisplayName', 'Pressure Outlet');

plot(cellCounts(4:6), normalisedL2(4:6), '--s', ...
    'Color', meshColours(2,:), 'MarkerSize', 9, ...
    'MarkerEdgeColor', black, 'MarkerFaceColor', meshColours(2,:), ...
    'DisplayName', 'Outflow');

set(gca, 'XScale', 'log');
formatAxes('Number of cells [-]', 'Normalised L_2 velocity error [%]');
legend('Location', 'best', 'Box', 'off');
exportFigure(fig, figureFolder, 'mesh_error_L2');

%% Figure 6: outlet temperature profile, when the file is available
if isfile(thermalFile)
    [Ytemp, Tout] = readTwoColumnFile(thermalFile);

    if any(Ytemp < -1e-8 | Ytemp > H + 1e-8)
        warning('Coordinates in %s are outside 0 to %.4f m; plot skipped.', ...
            thermalFile, H);
    else
        fig = figure('Name', 'Thermal outlet temperature profile', ...
            'Position', [100 100 900 600]);
        markerIndex = chooseMarkerIndices(numel(Ytemp), 5, 1);
        plot(Ytemp, Tout, '-o', 'Color', meshColours(1,:), ...
            'MarkerIndices', markerIndex, 'MarkerSize', 8, ...
            'MarkerEdgeColor', black, 'MarkerFaceColor', meshColours(1,:));
        formatAxes('Transverse coordinate, Y [m]', 'Static temperature, T [K]');
        xlim([0 H]);
        exportFigure(fig, figureFolder, 'thermal_outlet_temperature');
    end
else
    warning('Optional file %s was not found; thermal profile plot was skipped.', ...
        thermalFile);
end

%% Figure 7: centreline pressure distribution, when the file is available
if isfile(pressureFile)
    [X, P] = readTwoColumnFile(pressureFile);

    if any(X < -1e-8 | X > L + 1e-8)
        warning('Coordinates in %s are outside 0 to %.4f m; plot skipped.', ...
            pressureFile, L);
    else
        fig = figure('Name', 'Centreline static pressure distribution', ...
            'Position', [100 100 900 600]);
        plot(X, P, '-', 'Color', meshColours(1,:), 'LineWidth', 1.8);
        formatAxes('Streamwise coordinate, x [m]', 'Static pressure [Pa]');
        xlim([0 L]);
        exportFigure(fig, figureFolder, 'centerline_pressure_medium_PO');
    end
else
    warning('Optional file %s was not found; pressure plot was skipped.', ...
        pressureFile);
end

%% Mesh statistics from the Fluent mesh reports
Mesh = {'Coarse'; 'Medium'; 'Fine'};
NominalElementSize = [0.0005; 0.00025; 0.000125];
Cells = [8000; 32000; 128000];
Faces = [16240; 64480; 256960];
Nodes = [8241; 32481; 128961];
MinimumOrthogonalQuality = [0.999691; 0.998463; 0.998801];
MaximumAspectRatio = [1.44151; 1.46573; 1.45933];

MeshStatistics = table(Mesh, NominalElementSize, Cells, Faces, Nodes, ...
    MinimumOrthogonalQuality, MaximumAspectRatio);
disp('MESH STATISTICS');
disp(MeshStatistics);
writetable(MeshStatistics, fullfile(pwd, 'Mesh_Statistics.csv'));

fprintf('\nPost-processing finished.\n');
fprintf('Figures (PNG and FIG) are in: %s\n', figureFolder);
fprintf('Reference values and result tables are saved in the simulation folder.\n');

%% Functions used by the post-processing sections
function [x, values] = readTwoColumnFile(fileName)
    if ~isfile(fileName)
        error('Required file not found: %s. Check the MATLAB Current Folder.', ...
            fileName);
    end

    data = readmatrix(fileName);

    if size(data, 2) < 2
        error('File %s must contain at least two numeric columns.', fileName);
    end

    data = data(:, 1:2);
    data = data(all(isfinite(data), 2), :);

    if isempty(data)
        error('No numeric two-column data could be read from %s.', fileName);
    end

    data = sortrows(data, 1);
    x = data(:, 1);
    values = data(:, 2);
end

function formatAxes(xLabelText, yLabelText)
    xlabel(xLabelText, 'FontName', 'Times New Roman', 'FontSize', 16);
    ylabel(yLabelText, 'FontName', 'Times New Roman', 'FontSize', 16);

    ax = gca;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 14;
    ax.LineWidth = 1.0;
    ax.TickDir = 'out';
    ax.Box = 'on';
    ax.Layer = 'top';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.GridAlpha = 0.15;
    ax.MinorGridAlpha = 0.07;

    grid on;
    grid minor;
end

function index = chooseMarkerIndices(nPoints, maxMarkers, phase)
    if nPoints <= maxMarkers
        index = 1:nPoints;
        return;
    end

    % Offset marker positions slightly between curves to keep them distinct.
    lowerFraction = 0.08 + (phase - 1) * 0.08;
    upperFraction = min(0.98, 0.82 + (phase - 1) * 0.08);
    fractions = linspace(lowerFraction, upperFraction, maxMarkers);
    index = unique(round(1 + fractions * (nPoints - 1)));
end

function exportFigure(figHandle, folder, baseName)
    print(figHandle, fullfile(folder, [baseName '.png']), '-dpng', '-r300');
    savefig(figHandle, fullfile(folder, [baseName '.fig']));
end
