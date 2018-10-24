function bssePrimedTrackValidateCellsInEmbryo(imarisApplicationId)
% Validate cell segmentation in embryos.
%
% Code for the paper:
%
% Welling et al. "High fidelity lineage tracing in mouse pre-implantation 
% embryos using primed conversion of photoconvertible proteins".
%
% This Imaris XTension required IceImarisConnector to run.
% See: https://github.com/aarpon/IceImarisConnector
%
%    <CustomTools>
%      <Menu>
%       <Submenu name="BSSE">
%        <Submenu name="Primed Track">
%         <Item name="Validate cells in embryo" icon="Matlab">
%          <Command>MatlabXT::bssePrimedTrackValidateCellsInEmbryo(%i)</Command>
%         </Item>
%        </Submenu>
%       </Submenu>
%      </Menu>
%    </CustomTools>
%
% Aaron Ponti (BSSE) 2017, 2018

% Instantiate IceImarisConnector object
conn = IceImarisConnector(imarisApplicationId);

% Is there something loaded?
mDataSet = conn.mImarisApplication.GetDataSet();
if isempty(mDataSet)
    return
end

% Get the spots
spots = conn.getAllSurpassChildren(1, 'Spots');
if isempty(spots)
    uiwait(errordlg('No Spots objects found in the scene!'));
    return;
end

if numel(spots) < 2
    uiwait(errordlg('You need two Spots objects!'));
    return;
end

% Get spot names
spotNames = cell(1, numel(spots));
for i = 1 : numel(spots)
    spotNames{i} = char(spots{i}.GetName());
end

% Ask the user to specify the Spots object with all the cells
[s, v] = listdlg('PromptString', ...
    'Plese pick the Spots object with all the cells',...
    'SelectionMode', 'single', ...
    'ListSize', [600 100], ...
    'ListString', spotNames);
if v == 0
    return
end

% Get the "all cells" Spots object
allSpots = spots{s};
spots(s) = [];
spotNames(s) = [];

if numel(spots) == 1
    % Get the "activated cells" Spots object
    actSpots = spots{1};
else
    % Ask the user to specify the Spots object with the converted cells
    [s, v] = listdlg('PromptString', ...
        'Plese pick the Spots object with the converted cells',...
        'SelectionMode', 'single', ...
        'ListSize', [600 100], ...
        'ListString', spotNames);
    if v == 0
        return
    end
    
    % Get the "activated cells" Spots object
    actSpots = spots{s};
    clear('spots');
    clear('spotNames');
end

% Ask the user to specify the min distance between cell centers to
% resolve double-segmentations
answer = inputdlg( ...
    ['Please specify the min distance between cell centers to ', ...
    'resolve overlapping segmentations'], ...
    'Question', 1, {'2.0'});
if isempty(answer)
    return
end
minCellDist = str2double(answer{1});
if isnan(minCellDist) || minCellDist <= 0
    uiwait(errordlg('The max distance must be a positive scalar.'));
    return
end

% Ask the user to specify the factor that multiplies the estimated diameter
% of the embryo
answer = inputdlg( ...
    ['The algorithm is likely to underestimate the diameter of the ', ...
    'embryo, if the cells are not homogeneously distributed. ', ...
    'Please specify a multiplicative factor'], ...
    'Question', 1, {'2.0'});
if isempty(answer)
    return
end
embryoMultFactor = str2double(answer{1});
if isnan(embryoMultFactor) || embryoMultFactor <= 0
    uiwait(errordlg('The factor must be a positive scalar.'));
    return
end

% Ask the user to specify the max distance between cell centers
answer = inputdlg( ...
    ['Please specify the max allowed distance between cell centers ', ...
    'from both channels for them to be considered overlapping'], ...
    'Question', 1, {'1.0'});
if isempty(answer)
    return
end
maxCellDist = str2double(answer{1});
if isnan(maxCellDist) || maxCellDist <= 0
    uiwait(errordlg('The max distance must be a positive scalar.'));
    return
end

% Ask the user to specify the max distance to activate cell centers of
% previous time point for restoring discarded ones
answer = inputdlg( ...
    ['Please specify the max allowed distance between centers ', ...
    'of consecutive time points for recovering erroneously discarded ', ...
    'activated cells (set to 0 to disable)'], ...
    'Question', 1, {'2.0'});
if isempty(answer)
    return
end
maxDistToPreviousTimePoint = str2double(answer{1});
if isnan(maxDistToPreviousTimePoint) || maxDistToPreviousTimePoint < 0
    uiwait(errordlg( ...
        'The max distance must be a positive scalar (or 0 to disable).'));
    return
end

% Waitbar
hWaitbar = waitbar(0, 'Processing time points');

% Results
newActPos = [];
newActTime = [];
newActRadii = [];
newAllPos = [];
newAllTime = [];
newAllRadii = [];
newDiffAllPos = [];
newDiffAllTime = [];
newDiffAllRadii = [];
nTimepoints = mDataSet.GetSizeT();

numCellsPerTimepointBeforeAnyCorrection = zeros(1, nTimepoints);
numCellsPerTimepointAfterOverlappingSegmentation= zeros(1, nTimepoints);
numCellsPerTimepointAfterColocCorrection = zeros(1, nTimepoints);
numCellsPerTimepointAfterTimeCorrection = zeros(1, nTimepoints);

% Keep track of the last validated cell positions
lastActCellsPos = [];

% Extract information
allSpotPositionsXYZ = allSpots.GetPositionsXYZ();
allSpotTimeIndices = allSpots.GetIndicesT();
allSpotRadiiXYZ = allSpots.GetRadiiXYZ();
actSpotPositionsXYZ = actSpots.GetPositionsXYZ();
actSpotTimeIndices = actSpots.GetIndicesT();
actSpotRadiiXYZ = actSpots.GetRadiiXYZ();

% Process time points
for timepoint = 1 : nTimepoints

    % Find spots for current time point
    indicesAll = allSpotTimeIndices == (timepoint - 1);
    indicesAct = actSpotTimeIndices == (timepoint - 1);
    
    % Get relevant data for current timepoint
    currentAllPos = allSpotPositionsXYZ(indicesAll, :);
    currentAllTime = allSpotTimeIndices(indicesAll, :);
    currentAllRadii = allSpotRadiiXYZ(indicesAll, :);
    
    currentActPos = actSpotPositionsXYZ(indicesAct, :);
    currentActTime = actSpotTimeIndices(indicesAct, :);
    currentActRadii = actSpotRadiiXYZ(indicesAct, :);

    % Get the number of cells in this time point before correction
    numCellsPerTimepointBeforeAnyCorrection(timepoint) = size(currentActPos, 1);

    % ---------------------------------------------------------------------
    %
    % Remove overlapping cells from double-segmentations
    %
    % ---------------------------------------------------------------------
    
    D = pdist2(currentAllPos, currentAllPos);
    for i = 1 : size(D, 1)
        for j = 1 : size(D, 2)
            if j <= i
                D(i, j) = Inf;
            end
        end
    end
    [y, toDiscard] = find(D < minCellDist);
    if ~isempty(y)
        fprintf(1, ...
            'Frame %d: %d cell(s) discarded because overlapping\n', ...
            timepoint, ...
            numel(toDiscard));
        
        for i = 1 : numel(y)
           
            % We replace one of the spots with the average position of
            % the two overlapping one, and then we discard the other
            currentAllPos(y(i), :) = mean([ ...
                currentAllPos(y(i), :);
                currentAllPos(toDiscard(i), :) ...
                ], 1);
        end
        
        currentAllPos(toDiscard, :) = [];
        currentAllTime(toDiscard, :) = [];
        currentAllRadii(toDiscard, :) = [];
    end
    
    D = pdist2(currentActPos, currentActPos);
    for i = 1 : size(D, 1)
        for j = 1 : size(D, 2)
            if j <= i
                D(i, j) = Inf;
            end
        end
    end
    [y, toDiscard] = find(D < minCellDist);
    if ~isempty(y)
        fprintf(1, ...
            'Frame %d: %d cell(s) discarded because overlapping\n', ...
            timepoint, ...
            numel(toDiscard));
        
        for i = 1 : numel(y)
           
            % We replace one of the spots with the average position of
            % the two overlapping one, and then we discard the other
            currentActPos(y(i), :) = mean([ ...
                currentActPos(y(i), :);
                currentActPos(toDiscard(i), :) ...
                ], 1);
        end
        
        currentActPos(toDiscard, :) = [];
        currentActTime(toDiscard, :) = [];
        currentActRadii(toDiscard, :) = [];       
    end

    % Store the number of spots after the cell double-segmentation correction
    numCellsPerTimepointAfterOverlappingSegmentation(timepoint) = ...
        size(currentActPos, 1);
    
    % ---------------------------------------------------------------------
    %
    % Remove cells from the all cells set if they are outside the embryo
    %
    % ---------------------------------------------------------------------
    
    % Estimate the radius of the embryo as embryoMultFactor * the median
    % of the distances between cells
    allDist = double(pdist2(currentAllPos, currentAllPos));
    estDiameter = embryoMultFactor * median(allDist(:));
    
    % Discard cells that are more than estDiameter away from their center
    % of mass
    cellsToDiscard = find( ...
        pdist2(currentAllPos, median(currentAllPos, 1)) > estDiameter);
    if ~isempty(cellsToDiscard)
        fprintf(1, ...
            'Frame %d: %d cell(s) discarded because outside of embryo\n', ...
            timepoint, ...
            numel(cellsToDiscard));
    end
    currentAllPos(cellsToDiscard, :) = [];
    currentAllTime(cellsToDiscard) = [];
    currentAllRadii(cellsToDiscard, :) = [];
    
    % Append the valid ones to growing list
    newAllPos = cat(1, newAllPos, currentAllPos);
    newAllTime = cat(1, newAllTime, currentAllTime);
    newAllRadii = cat(1, newAllRadii, currentAllRadii);

    % ---------------------------------------------------------------------
    %
    % Remove cells from the converted cells set if they do not 
    % colocalize with cells from the all cells set
    %
    % ---------------------------------------------------------------------
    
    % Keep track of the cells before the filter is applied
    currentActPosBeforeColoc = currentActPos;
    currentActTimeBeforeColoc = currentActTime;
    currentActRadiiBeforeColoc = currentActRadii;
    
    % Calculate the distances
    D = double(pdist2(currentActPos, currentAllPos));
    D(D > maxCellDist) = -1;
    
    % Match the cells
    try
        L = lap(D, -1, 0, 1);
        L = L(1 :  size(currentActPos, 1));
    catch
        L = [];
    end
    validMatches = find(L <= size(currentAllPos, 1));

    if numel(validMatches) < size(currentActPos, 1)
        fprintf(1, ...
            'Frame %d: %d cell(s) discarded because not colocalizing\n', ...
            timepoint, ...
            size(currentActPos, 1) - numel(validMatches));        
    end
    
    % Extract the spots to store
    currentActPos = currentActPos(validMatches, :);
    currentActTime = currentActTime(validMatches, :);
    currentActRadii = currentActRadii(validMatches, :);
    
    % Store the number of spots after the cell colocalization correction
    numCellsPerTimepointAfterColocCorrection(timepoint) = ...
        size(currentActPos, 1);
    
    % We also store the cells from the all cells set that do not overlap
    % with cells from the converted cells set
    currentDelIndices = L(L <= size(currentAllPos, 1));
    currentDiffAllPos = currentAllPos;
    currentDiffAllPos(currentDelIndices, :) = [];
    currentDiffAllTime = currentAllTime;
    currentDiffAllTime(currentDelIndices) = [];
    currentDiffAllRadii = currentAllRadii;
    currentDiffAllRadii(currentDelIndices, :) = [];

    % Append the valid ones to growing list
    newDiffAllPos = cat(1, newDiffAllPos, currentDiffAllPos);
    newDiffAllTime = cat(1, newDiffAllTime, currentDiffAllTime);
    newDiffAllRadii = cat(1, newDiffAllRadii, currentDiffAllRadii);
    
    % ---------------------------------------------------------------------
    %
    % Recover cells from the converted cells set if they do not 
    % colocalize with cells from the all cells set BUT they are close
    % to a cell from the converted cells set in the previous time point.
    %
    % ---------------------------------------------------------------------
    
    % Compare the validated cells with previous timepoint
    if maxDistToPreviousTimePoint ~= 0
        if numel(validMatches) < size(lastActCellsPos, 1)
            
            % We lost cells in this time point compared to the
            % previous
            %
            % Compare discarded positions from this timepoint with
            % positions in the previous
            discardedCurrentActPos = currentActPosBeforeColoc;
            discardedCurrentActPos(validMatches, :) = [];
            discardedCurrentValidActTime = currentActTimeBeforeColoc;
            discardedCurrentValidActTime(validMatches, :) = [];
            discardedCurrentValidActRadii = currentActRadiiBeforeColoc;
            discardedCurrentValidActRadii(validMatches, :) = [];
            
            % Is there one that has a partner within a certain
            % distance
            D = double(pdist2(discardedCurrentActPos, lastActCellsPos));
            [y, ~] = find(D < maxDistToPreviousTimePoint);
            
            if ~isempty(y)
                fprintf(1, ...
                    'Frame %d: %d discarded cell(s) recovered\n', ...
                    timepoint, ...
                    numel(y));
                
                currentActPos = cat(1, currentActPos, ...
                    discardedCurrentActPos(y, :));
                currentActTime = cat(1, currentActTime, ...
                    discardedCurrentValidActTime(y, :));
                currentActRadii = cat(1, currentActRadii, ...
                    discardedCurrentValidActRadii(y, :));
            else
                fprintf(1, ...
                    'Frame %d: could not recover discarded cells\n', ...
                    timepoint);
            end
        end
    end
    
    % Append the valid ones to growing list
    newActPos = cat(1, newActPos, currentActPos);
    newActTime = cat(1, newActTime, currentActTime);
    newActRadii = cat(1, newActRadii, currentActRadii);
    
    % Store the number of cells after correcton
    numCellsPerTimepointAfterTimeCorrection(timepoint) = ...
        size(currentActPos, 1);
    
    % Store the positions for next round
    lastActCellsPos = currentActPos;
    
    % Inform if needed
    if numCellsPerTimepointAfterTimeCorrection(timepoint) ~= ...
            numCellsPerTimepointBeforeAnyCorrection(timepoint)
        fprintf(1, 'Frame %d: number of cells %d -> %d\n', ...
            timepoint, ...
            size(actSpotPositionsXYZ(indicesAct, :), 1), ...
            size(currentActPos, 1));
    end
    
    % Update the waitbar
    waitbar(timepoint / nTimepoints, hWaitbar);
    
end

% -------------------------------------------------------------------------
%
% Create a quality control figure
%
% -------------------------------------------------------------------------

figure;
plot(numCellsPerTimepointBeforeAnyCorrection, 'k', ...
    'DisplayName', 'Before any correction');
hold on;
plot(numCellsPerTimepointAfterOverlappingSegmentation, 'm', ...
    'DisplayName', 'After cell double segmentation');
hold on;
plot(numCellsPerTimepointAfterColocCorrection, 'b', ...
    'DisplayName', 'After coloc correction');
if maxDistToPreviousTimePoint ~= 0
    plot(numCellsPerTimepointAfterTimeCorrection, 'r', ...
        'DisplayName', 'After time correction');
end
title([char(allSpots.GetName()), ' vs. ', char(actSpots.GetName())]);
xlabel('Time point');
ylabel('# of spots');
legend('show', 'Location', 'bestoutside');
grid on

% Close the waitbar
close(hWaitbar);

% -------------------------------------------------------------------------
%
% Create the new spot objects and add them to the scene
%
% -------------------------------------------------------------------------

% Corrected all cells set
conn.createAndSetSpots(newAllPos, newAllTime, newAllRadii, ...
    ['Processed ', char(allSpots.GetName())], [rand(1, 3), 0]);

% Corrected activated cells set
conn.createAndSetSpots(newActPos, newActTime, newActRadii, ...
    ['Processed ', char(actSpots.GetName())], [rand(1, 3), 0]);

% All cells - Activated cells set
conn.createAndSetSpots(newDiffAllPos, newDiffAllTime, newDiffAllRadii, ...
    ['Processed ', char(allSpots.GetName()), ' - Processed ', ...
    char(actSpots.GetName())], [rand(1, 3), 0]);
