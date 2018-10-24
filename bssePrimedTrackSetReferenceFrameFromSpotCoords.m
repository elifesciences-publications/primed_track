function bssePrimedTrackSetReferenceFrameFromSpotCoords(aImarisId)
% Extracts a position and an orientation and sets a reference frame.
%
% Code for the paper:
%
% Welling et al. "High fidelity lineage tracing in mouse pre-implantation 
% embryos using primed conversion of photoconvertible proteins".
%
% This Imaris XTension required IceImarisConnector to run.
% See: https://github.com/aarpon/IceImarisConnector
%
%  Installation:
%
%  Copy this file into the XTensions folder in the Imaris installation 
%  directory. You will find this function in the Image Processing menu.
%
%    <CustomTools>
%      <Menu>
%       <Submenu name="BSSE">
%        <Submenu name="Primed Track">
%         <Item name="Set reference frames from cells" icon="Matlab">
%           <Command>MatlabXT::bssePrimedTrackSetReferenceFrameFromSpotCoords(%i)</Command>
%         </Item>
%        </Submenu>
%       </Submenu>
%      </Menu>
%    </CustomTools>
%
% Aaron Ponti (BSSE) 2017, 2018

% Initialize IceImarisConnector
conn = IceImarisConnector(aImarisId);

% Get the spots
spots = conn.getAllSurpassChildren(1, 'Spots');

if numel(spots) < 2
    error('At least 2 Spots objects expected.');
end

spotsNames = cell(1, numel(spots));
for i = 1 : numel(spots)
    spotsNames{i} = char(spots{i}.GetName());
end

[s, v] = listdlg('PromptString', ...
    'Select spots for estimating embryo center point:',...
    'ListSize', [400, 300], 'SelectionMode', 'single', ...
    'ListString', spotsNames);
if v == 0
    return;
end

spotsCenter = spots{s};

% Remove the selected spot
spots(s) = [];
spotsNames(s) = [];

if numel(spots) == 1
    spotsOrientation = spots{1};
else

    [s, v] = listdlg('PromptString', ...
        'Select spots for estimating the orientation of the embryo:',...
        'SelectionMode', 'single', 'ListString', spotsNames);
    if v == 0
        return;
    end

    spotsOrientation = spots{s};

end

% Now process all timepoints
sz = conn.getSizes();
nTimepoints = sz(5);

% Get all coordinates, radii and timepoints for both Spots objects
centerValues = spotsCenter.Get();
orientValues = spotsOrientation.Get();

centerCoords = zeros(nTimepoints, 3);
orientCoords = zeros(nTimepoints, 3);
orientVector = zeros(nTimepoints, 3);
quaternions = zeros(nTimepoints, 4);

extends = conn.getExtends();
centerOfDataset = 0.5 * (extends(1 : 2 : 5) + extends(2 : 2 : 6));
originOfDataset = extends(1 : 2 : 5);

refVector = [0 1 0];
for t = 0 : nTimepoints - 1

    centerIndices = centerValues.mIndicesT == t;
    orientIndices = orientValues.mIndicesT == t;

    centerPos = mean(centerValues.mPositionsXYZ(centerIndices, :), 1);
    if any(isnan(centerPos))
        centerPos = centerOfDataset;
    end
    orientPos = mean(orientValues.mPositionsXYZ(orientIndices, :), 1);
    if any(isnan(orientPos))
        orientPos = centerOfDataset;
    end

    centerCoords(t + 1, 1 : 3) = centerPos;
    orientCoords(t + 1, 1 : 3) = orientPos;
    vector = orientPos - centerPos;
    orientVector(t + 1, 1 : 3) = vector;
    if all(vector == 0)
        quaternions(t + 1, 1 : 4) = [0 0 0 1]; % No rotation
    else
        quaternions(t + 1, 1 : 4) = conn.calcRotationBetweenVectors3D(refVector, vector);
    end

end

% Prevent Imaris 8 from crashing (!)
epsilon = 1e-4;
indx = find(( ...
    abs(quaternions(:, 1) + 1) < epsilon | ...
    abs(quaternions(:, 2) + 1) < epsilon | ...
    abs(quaternions(:, 3) + 1) < epsilon ) & ...
    abs(quaternions(:, 4) < epsilon));
for i = 1 : numel(indx)
    % No rotation
    quaternions(indx(i), :) = [0 0 0 1];
end

% The coordinates of a reference frame are RELATIVE to the origin of the 
% dataset -- in contrast to all other objects
centerCoords = centerCoords - ...
    repmat(originOfDataset, size(centerCoords, 1), 1);

% Create a reference frame
referenceFrame = conn.mImarisApplication.GetFactory().CreateReferenceFrames();
referenceFrame.SetKeyFramesPositionsXYZT(0:(nTimepoints - 1), centerCoords);
referenceFrame.SetKeyFramesByQuaternionTime(0:(nTimepoints - 1), quaternions);

% Add the reference frame back to the scene
conn.mImarisApplication.GetSurpassScene().AddChild(referenceFrame, -1);
