function growthCurves
% This function processes tables of Algal Growth Curves data and presents
% them as combined plots.
%
% It requires the UnivarScatter function from the Mathworks file exchange:
% https://uk.mathworks.com/matlabcentral/fileexchange/54243-univarscatter
% Written by Manuel Lera Ramírez
% 
% It requires findchangepts function from the Signal Processing Toolbox.
% 
% It requires tinv function from the Statistics Toolbox.
%
% The function produces graphs that are saved as PNG files.


% Jakub Nedbal, King's College London, 2019
% File Created: 18. Nov, 2019
%
% The file is distributed under the BSD License
% Copyright 2019 Jakub Nedbal
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are
% met:
%
% 1. Redistributions of source code must retain the above copyright notice, 
%    this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright 
%    notice, this list of conditions and the following disclaimer in the 
%    documentation and/or other materials provided with the distribution.
%
% 3. Neither the name of the copyright holder nor the names of its 
%    contributors may be used to endorse or promote products derived from 
%    this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
% THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% Global variable used in plotting the data
global graphPar

% Populate the cell with the parameters of the plot.
% Each dataset gets one line in the cell.
% The columns have the following meaning:
%   1. scatter plot marker shape
%   2. linear regression line style
%   3. Color of the dataset
%   4. Legend label for each dataset
graphPar = {'o', ':', [0, 0.4470, 0.7410], 'medium'; ...
            'x', '-.', [0.8500, 0.3250, 0.0980], 'low'};% ...
            %'s', [0.9290, 0.6940, 0.1250], 'high'};

% Cell with the names of the XLSX (Excel files) containing the cell counts
% of the various dates.
% Each file contains one experiment. Do not include the '.xlsx' extension
powers = {'medium', 'low'};

% Get rid of any old figures
close all

% Create an array of figure and axes handles
hf = gobjects(3, 2);
% Form the figures with axes used throughout
for i = 1 : size(hf, 1)
    hf(i, 1) = figure('Units', 'pixels', 'Position', [0 0 750, 750]);
    hf(i, 2) = axes('Units', 'pixels');
end

% Run through all datasets one-by-one
for pwr = 1 : numel(powers)

    % Load the growth data
    % The sheet name in the table is critical to parsing the XLSX file.
    table.Cvulgaris = readtable([powers{pwr} '.xlsx'], ...
                                'Sheet', 'C. vulgaris');
    table.Dquadricauda = readtable([powers{pwr} '.xlsx'], ...
                                   'Sheet', 'D. quadricauda');

    % Create structs for the data in the XLSX tables
    Chl = struct('date', [], ...
                 'density', [], ...
                 'avgDensity', [], ...
                 'temperature', []);
    DDq = Chl;
    for i = {'Cvulgaris', 'Dquadricauda'}
        % extract the string with the organism
        org = i{1};
        % Get the dates
        fns = fieldnames(table.(org));
        % read those fields which contain dates
        dateIn = find(~cellfun(@isempty, ...
                               regexp(fns, 'x\d{1,2}_\d{1,2}_\d{4}')));
        % Create a struct with all the results for the dataset
        result = struct('date', [], 'doub', [], 'quad', [], 'oct', [], ...
                        'count', [], 'ss', [], 'density', [], ...
                        'avgDensity', []);
        % Go through the data day-by-day
        for j = 1 : numel(dateIn)
            % Read the date and convert it to a Matlab number format
            result(j).date = datenum(fns{dateIn(j)}, 'xdd_mm_yyyyHH_MM');
            % Calculate the columns for the given date
            cols = dateIn(j) + (0 : (min(diff(dateIn)) - 1));
            % Read the table for a given date.
            block = table.(org)(1 : end, cols).Variables;
            % Go through the results column-by-column
            for k = 1 : size(block, 2)
                % Check what is in the title of each column
                switch block{1, k}
                    case 'Square'   % Hemocytometer square coordinates
                        % Check if there is temperature in the left row
                        if ~isnan(str2double(block{end, k}))
                            % If there is temperature, store is as a number
                            temperature = str2double(block{end, k});
                        else
                            % Without any temperature data, store NaN
                            temperature = NaN;
                        end
                    case 'Doublets' % Cell doublet coenobia count
                        % Convert to numbers
                        result(j).doub = cellfun(@(x) str2double(x), ...
                                                 block(2 : end, k));
                    case 'Quadruples'   % Cell quadruples coenobia count
                        % Convert to numbers
                        result(j).quad = cellfun(@(x) str2double(x), ...
                                                 block(2 : end, k));
                    case 'Octuples'     % Cell octuples coenobia count
                        % Convert to numbers
                        result(j).oct = cellfun(@(x) str2double(x), ...
                                                block(2 : end, k));
                    case 'Count'        % Chlorella cell count
                        % Convert to numbers
                        result(j).count = cellfun(@(x) str2double(x), ...
                                                block(2 : end, k));
                    case 'Sub-squares'  % Chlorella subsquares counted
                        % Convert to numbers
                        result(j).ss = cellfun(@(x) str2double(x), ...
                                               block(2 : end, k));
                    otherwise
                        error('Cannot parse the table')
                end
            end
            % Calculate the totals for both organisms
            switch org
                case 'Cvulgaris'
                    % Work out which results are not NaNs
                    countOK = ~isnan(result(j).count);
                    ssOK = ~isnan(result(j).ss);
                    % Check that NaNs in counts and sub-squares match
                    % First create error message
                    errmsg = ['Empty fields in C. vulgaris ', ...
                              'do not match on %s.'];
                    % Check the empty fields match
                    assert(isequal(countOK, ssOK), ...
                           errmsg, datestr(result(j).date))
                    % Calculate the density
                    % The density is 1e4 times the cells in a square. The
                    % density is scaled by the number of squares counted.
                    result(j).density = ...
                        1e4 * result(j).count(countOK) * ...
                              20 ./ result(j).ss(ssOK);
                    % Calculate the mean density
                    result(j).avgDensity = mean(result(j).density);
                    % Store the data
                    % Date
                    Chl.date(j) = result(j).date;
                    % Day
                    Chl.day{j} = num2str(round(Chl.date(j) - Chl.date(1)));
                    % Density
                    Chl.density(1 : numel(result(j).density), j) = result(j).density;
                    % Average density
                    Chl.avgDensity(j) = result(j).avgDensity;
                    % Culture temperature
                    Chl.temperature(j) = temperature;
                case 'Dquadricauda'
                    % Work out which results are not NaNs
                    doubOK = ~isnan(result(j).doub);
                    quadOK = ~isnan(result(j).quad);
                    octOK = ~isnan(result(j).oct);
                    % Check that NaNs in counts and sub-squares match
                    % First create error message
                    errmsg = ['Empty fields in D. quadricauda ', ...
                              'do not match on %s.'];
                    % Check the empty fields match
                    assert(isequal(countOK, ssOK), ...
                           errmsg, datestr(result(j).date))
                    % Calculate the density.
                    % The density is 1e4 times the cells in a square. The
                    % density is calculated by multiplying the number of
                    % coenobia with 2 cells by two, coenobia with 4 cells
                    % by four, and coenobia with 8 cells by eight.
                    result(j).density = ...
                        1e4 * (2 * result(j).doub(doubOK) + ...
                               4 * result(j).quad(quadOK) + ...
                               8 * result(j).oct(octOK));
                    % Calculate the mean density
                    result(j).avgDensity = mean(result(j).density);
                    % Store the data
                    % Date
                    DDq.date(j) = result(j).date;
                    % Day
                    DDq.day{j} = num2str(round(DDq.date(j) - DDq.date(1)));
                    % Density
                    DDq.density(1 : numel(result(j).density), j) = result(j).density;
                    % Average density
                    DDq.avgDensity(j) = result(j).avgDensity;
                    % Culture temperature
                    DDq.temperature(j) = temperature;
            end
        end
    end
    % Fill in gaps in the data
    % Sometimes the counting is not done on every day. These missing dates
    % need to be inserted into the matrices of data as NaNs
    gaps = diff(cellfun(@str2double, DDq.day)) - 1;
    % Find days when there is more than a day between itself and the next
    % experimental time point.
    gapsIn = find(gaps);
    % Go through all these days in the direction from the last gap to the
    % first gap in the data.
    for i = fliplr(gapsIn)
        % Fill in the names of the dates as an empty string ('') for days
        % when no measurements has taken place.
        DDq.day = horzcat(DDq.day(1 : i), ...
                          repmat({''}, 1, gaps(i)), ...
                          DDq.day(i + 1 : end));
        % Fill in the counts on the dates when experiments were skipped by
        % NaNs
        DDq.density = horzcat(DDq.density(:, 1 : i), ...
                              NaN(size(DDq.density, 1), gaps(i)), ...
                              DDq.density(:, i + 1 : end));
    end
    % Fill in gaps in the data in the Chlorella dataset. Do the same as
    % above
    gaps = diff(cellfun(@str2double, Chl.day)) - 1;
    gapsIn = find(gaps);
    for i = fliplr(gapsIn)
        Chl.day = horzcat(Chl.day(1 : i), ...
                          repmat({''}, 1, gaps(i)), ...
                          Chl.day(i + 1 : end));
        Chl.density = horzcat(Chl.density(:, 1 : i), ...
                              NaN(size(Chl.density, 1), gaps(i)), ...
                              Chl.density(:, i + 1 : end));
    end

    % Plot the results
    % This populates the graphs with the data and labels them accordingly
    % Plot Chlorella vulgaris data
    makePlot(Chl.day, ...       % X Tick labels
             Chl.density, ...   % Cell count data
             'C. vulgaris', ... % Species named in the graph title
             pwr, ...           % Current experiment index (power)
             hf(1, :))          % Handles to appropriate figure and axes

    % Plot Desmodesmus quadricauda data
    makePlot(DDq.day, ...       % X Tick labels
             DDq.density, ...   % Cell count data
             'D. quadricauda', ...  % Species named in the graph title
             pwr, ...           % Current experiment index (power)
             hf([2, 3], :), ... % Handles to appropriate figures and axes
             result)            % Data with coenobia cell counts
end


function makePlot(day, density, species, pwr, hfig, result)
% This function plots the data from the growth curve cell counting
% experiments
%   day:        X Tick labels
%   density:    Cell count data
%   species:    Species in the graph title
%   pwr:        Current experiment index (power)
%   hfig:       Handles to appropriate figures and axes
%   result:     Data with coenobia cell counts (only for D. quadricauda)

% global variable with graph parameters
global graphPar

% Switch to the appropriate figure
figure(hfig(1, 1))
% Switch to the appropriate axes
axes(hfig(1, 2))
% Plot the univariate scatter plots of the cell count results
UnivarScatter(density, ...                              % Data
              'Label', day, ...                         % X Tick labels
              'MarkerFaceColor', 'none', ...            % No Marker fill
              'MarkerEdgeColor', graphPar{pwr, 3}, ...  % Marker line color
              'PointStyle', graphPar{pwr, 1}, ...       % Marker shape
              'PointSize', 12);                         % Marker size

% Change the Y axis to logarithmic scale
set(gca, 'YScale', 'log')
% Provide X axis label
xlabel('Time [day]')
% Provide Y axis label
ylabel('Cell Density [ml^{-1}]')
% Title the graph, naming hte species
title(sprintf('{\\it %s} Culture Growth Timeline', species))
% Ensure the axes have consistant size
set(gca, 'Position', [100, 100, 600, 600])
% Change axes to square - the axes with the univariate scatter plots in
% logarithmic scale seem to misbehave and be buggy. Position parameter does
% not really change the position the way it should. using this command
% helps making a square-ish graph, though not a perfect square.
axis square
% Keep drawing over the existing data
hold on
% Change the day labels, which are strings into numbers
day = repmat(cellfun(@str2double, day), size(density, 1), 1);
% Find abrupt change in the growth curve, i.e. where the cells culture
% starts being limited by light or nutrients. Apply to logarithm of the
% data. The exponential growth phase is assumed to be linear on the
% logarithmic scale.
% Requires the signal processing toolbox
ipt = findchangepts(log10(density(~isnan(density))), 'Statistic', 'std');
% Calculate the linear regression of the data, between the first 
% measurement and the point of the abrupt change
[linCoef, errEst] = polyfit(day(1 : ipt), log10(density(1 : ipt)), 1);
% Calculate the linear regression values 
Yf = polyval(linCoef, day(1 : ipt), errEst);
% Plot the linear part of the curve
plot(1 + day(1 : ipt), 10 .^ Yf, ...    % Linear Regression data points
     graphPar{pwr, 2}, ...              % Line style
     'LineWidth', 2, ...                % Line thickness
     'Color', graphPar{pwr, 3});        % Line color
% Calculate the errors of the linear regression coefficients and the cell
% doubling time.
% The calculations below follow the instructions from the document below:
% http://physics.ujep.cz/~ehejnova/UTM/materialy_studium/linearni_regrese.pdf
%
% Calculate the 95 % confidence interval (CI) of the linear regression
% coefficient
% Select datapoints for the linear least sqaures from the linear section 
% of the cell count data
% X axis data (days of measurement)
x = day(1 : ipt);
% Y axis data (base 2 logarithm of the cell count
y = log2(density(1 : ipt));
% Slope of the linear regression
k = (ipt * sum(x .* y) - sum(x) * sum(y)) /...
    (ipt * sum(x .^ 2) - sum(x) ^ 2);
% Zero offset of the linear regression (not required)
% q = (sum(x .^ 2) * sum(y) - sum(x) * sum(x .* y)) / ...
%     (ipt * sum(x .^ 2) - sum(x) ^ 2);
% temporary variable to calculate the linear regression slope error
S0 = sum(y .^ 2) - sum(y) ^ 2 / ipt - ...
     k * (sum(x .* y) - sum(x) * sum(y) / ipt);
% Standard deviation of the slope
sigmak = sqrt(S0 / ((ipt - 2) * (sum(x .^ 2) - sum(x) ^ 2 / ipt)));
% 95% confidence interval
% This requires the tinv from the Statistics toolbox to calculate the
% Inverse of Student's T cumulative distribution function
CI95k = sigmak * tinv(0.95, ipt - 2);

% Doubling time is the inverse of the linear regression slope coefficient
tD = 1 / k;
%sigmatD = tD * sigmak / k;
% Error of the doubling sime is calculated according to error propagation
% rules
CI95tD = tD * CI95k / k;

% Print the fit results
fprintf('%s doubling time (tD): (%.2f %s %.2f) days.\n', ...
        species, ...    % Organism species
        tD, ...         % Doubling time
        char(177), ...  % Plus minus symbol
        CI95tD)         % 95% confidence interval of the doubling time
fprintf('%s doubling time (tD): (%.0f %s %.0f) hours.\n', ...
        species, ...    % Organism species
        24 * tD, ...    % Doubling time in hours
        char(177), ...  % Plus minus symbol
        24 * CI95tD)    % 95% CI of the doubling time in hours

% Remove whiskers
delete(findobj(get(gca, 'Children'), 'Type', 'Rectangle'))
% Change the short black bars into colored bars that are longer
% Find the handles of all these bars
hl = findobj(get(gca, 'Children'), 'Type', 'Line', ...
                                   'LineStyle', '-', ...
                                   'Color', 'k');
% Change the color line-by-line
for i = 1 : numel(hl)
    % Extend the mean lines and change their color to match the scatter
    % plot color
    set(hl(i), 'XData', get(hl(i), 'XData') + [-0.2, 0.2], ...
               'Color', graphPar{pwr, 3})
end

%% Create a legend
% Create an empty array of handles for unique graphical objects
hleg = gobjects(3, pwr);
% Create a cell array of strings for legend entries
sleg = cell(size(hleg));
for i = 1 : pwr
    % Get the handles of all scatter plots with the marker shape of the
    % current one
    hs = findobj(get(gca, 'Children'), ...
                 'Type', 'Scatter', ...     % Scatter plots
                 'Marker', graphPar{i, 1}); % Objects with the right marker
    % Store only the first identified scatter plot, they are all the same
    hleg(1, i) = hs(1);
    % Store the legend entry, including the illumination power
    sleg{1, i} = sprintf('Cell Count (%s)', graphPar{i, 4});
    
    % Get the handles of all lines with the line style of the current one
    hl = findobj(get(gca, 'Children'), ...
                 'Type', 'Line', ...            % Lines
                 'LineStyle', graphPar{i, 2});  % Right line style objects
    % Store only the first identified line, they are all the same
    hleg(2, i) = hl(1);
    % Store the legend entry, including the illumination power
    sleg{2, i} = sprintf('Linear Regression (%s)', graphPar{i, 4});
    % Get the handles of all lines with the solid line style and color of 
    % the current one
    hl = findobj(get(gca, 'Children'), ...
                 'Type', 'Line', ...        % Lines only
                 'LineStyle', '-', ...      % with solid line style
                 'Color', graphPar{i, 3});  % and the right color
    % Store only the first identified line, they are all the same
    hleg(3, i) = hl(1);
    % Store the legend entry, including the illumination power
    sleg{3, i} = sprintf('Mean (%s)', graphPar{i, 4});
end
% Create the legend in the lower right corner of the graph
legend(hleg(:), sleg(:), 'Location', 'SouthEast')
% Increase the size of the letters in the graph
set(gca, 'FontSize', 16)

% Create a filename from the species to save the figure
fname = sprintf('%s.png', regexprep(species, ' ', ''));
% Store the figure as a PNG file
saveas(gcf, fname, 'png')

% Only for Desmodesmus quadricauda, make a plot of the relative abundance
% of doublets, quadruplets, and octuplets
if nargin > 5
    % Switch to the appropriate figure
    figure(hfig(2, 1))
    % Switch to the appropriate axes
    axes(hfig(2, 2))
    % Calculate the number of cells in coenobia with doublets
    doub = mean([result.doub]) * 2;
    % Calculate the number of cells in coenobia with quadruplets
    quad = mean([result.quad]) * 4;
    % Calculate the number of cells in coenobia with octuplets
    oct = mean([result.oct]) * 8;
    % Calculate the total cells
    tot = doub + quad + oct;
    % Extract the days, which have experimental data, i.e. not NaNs
    day = day(1, ~isnan(day(1, :)));
    % Plot the relative abundance of cells in the different coenobia types
    plot(day, 100 * doub ./ tot, ...
         day, 100 * quad ./ tot, ...
         day, 100 * oct ./ tot)
    % Keep drawing over the existing data
    hold on
end