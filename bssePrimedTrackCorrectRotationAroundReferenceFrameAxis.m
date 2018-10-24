function bssePrimedTrackCorrectRotationAroundReferenceFrameAxis(aImarisId)
% Reorient a reference frame to standard orientation [0 1 0]
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
%         <Item name="Correct rotation around reference frame axis" icon="Matlab">
%           <Command>MatlabXT::bssePrimedTrackCorrectRotationAroundReferenceFrameAxis(%i)</Command>
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

if isempty(spots)
    error('No spots objects found.');
end

if numel(spots) == 1
    spot = spots{1};
else

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

    spot = spots{s};

end

% Get the reference frames
refs = conn.getAllSurpassChildren(1, 'ReferenceFrames');

if isempty(refs)
    error('No reference frames found!');
end

if numel(refs) == 1
    ref = refs{1};

else

    refsNames = cell(1, numel(refs));
    for i = 1 : numel(refs)
        refsNames{i} = char(refs{i}.GetName());
    end

    [s, v] = listdlg('PromptString', ...
        'Select the reference frame:',...
        'ListSize', [400, 300], 'SelectionMode', 'single', ...
        'ListString', refsNames);
    if v == 0
        return;
    end

    ref = refs{s};
end

% Now process all timepoints
sz = conn.getSizes();
nTimepoints = sz(5);

% Get the coordinate system matrix 
axisAngles = ref.GetKeyFramesOrientationsAxisAngleTime();

% Rotate around the Y axis
[~, ~, rotation_axis, ~] = conn.mapAxisAngleToRotationMatrix(...
    axisAngles.mAxisAngles.mAxis(1, :), axisAngles.mAxisAngles.mAngles(1));

% Initialize all target quaternions to be the standard [0 0 0 1] and same location
% as the source one
pos = ref.GetKeyFramesPositionsXYZT();
timeIndices = pos.mIndicesT;
coords = pos.mPositions;

% Prepare the target reference frame
referenceFrame = conn.mImarisApplication.GetFactory().CreateReferenceFrames();
referenceFrame.SetKeyFramesPositionsXYZT(timeIndices, coords);
all_axes = repmat(rotation_axis, nTimepoints, 1);
all_angles = zeros(nTimepoints, 1);
referenceFrame.SetKeyFramesByAxisAngleTime(timeIndices, all_axes, all_angles);

% Estimate the rotation around the axis
% Get all coordinates, radii and timepoints for both Spots objects
spotValues = spot.Get();

% Keep track of all individual correction angles; we will need to
% accumulate them in the end
all_angles = zeros(nTimepoints, 1);

% Keep track of the missing angles
missing = zeros(nTimepoints, 1);

for t = 0 : nTimepoints - 1

    % Get spots for curent time points
    spotIndices = spotValues.mIndicesT == t;

    spotCoords = spotValues.mPositionsXYZ(spotIndices, :);
    if isempty(spotCoords)
        missing(t + 1) = 1;
        angle = 0.0;
    else
      
        if t > 0
            % Find optimal rotation
            angle = findRotationAngle(spotCoords, previousSpotCoords, rotation_axis);
        else
            angle = 0.0;
        end
    end

    % Keep track of current spot coordinates for next timepoint
    previousSpotCoords = spotCoords;

    % Store current angle
    all_angles(t + 1) = angle;

end

% Now fill in the missing angles
missing_indices = find(missing);
if ~isempty(missing_indices)

    while numel(missing_indices) > 0

        i = 1;
            
        index = missing_indices(i);
            
        % Find the closest value below and above
        low_found = false;
        c = 1;
        while ~low_found
            prev = index - c;
            if prev < 1
                low = [];
                low_found = true;
            elseif ~isempty(find(missing_indices == prev, 1))
                c = c + 1;
                continue;
            else
                low = prev;
                low_found = true;
            end
        end
        
        high_found = false;
        d = 1;
        while ~high_found
            next = index + d;
            if next >= nTimepoints
                high = [];
                high_found = true;
            elseif ~isempty(find(missing_indices == next, 1))
                d = d + 1;
                continue;
            else
                high = next;
                high_found = true;
            end
        end
        
        if ~isempty(low) && ~isempty(high)
            missing_angle = (all_angles(low) + all_angles(high)) / (1 + high - low);
        elseif ~isempty(low) && isempty(high)
            missing_angle = all_angles(low);
        else
            missing_angle = all_angles(high);
        end
        
        all_angles(index) = missing_angle;
        
        missing_indices(i) = [];
    end
end

% Accumulate all angles
all_angles = mod(cumsum(all_angles), 2*pi);

% Calculate all rotations
quaternionTimes = ref.GetKeyFramesOrientationsQuaternionTime();
restQuaternion = quaternionTimes.mQuaternions(1, :);
quaternions = zeros(nTimepoints, 4);
for i = 1 : nTimepoints
    quat = conn.mapAxisAngleToQuaternion(rotation_axis, all_angles(i));
    quaternions(i, :) = IceImarisConnector.multiplyQuaternions(quat, restQuaternion);
end

% Create a new reference frame
referenceFrame = conn.mImarisApplication.GetFactory().CreateReferenceFrames();
referenceFrame.SetKeyFramesPositionsXYZT(timeIndices, coords);
referenceFrame.SetKeyFramesByQuaternionTime(timeIndices, quaternions);

% Add the reference frame back to the scene
conn.mImarisApplication.GetSurpassScene().AddChild(referenceFrame, -1);

% =========================================================================

function angle = findRotationAngle(spotCoords, previousSpotCoords, ref_axis)

% Center the coordinates around their center of mass
centerPos = mean(spotCoords, 1);
spotCoords = spotCoords - repmat(centerPos, size(spotCoords, 1), 1);
previousCenterPos = mean(previousSpotCoords, 1);
previousSpotCoords = previousSpotCoords - repmat(previousCenterPos, size(previousSpotCoords, 1), 1);

% Make sure the vector is normalized
ref_axis = IceImarisConnector.normalize(ref_axis);

step = pi/180;
angles = 0: step : 2*pi - step;
allD = zeros(1, numel(angles));

for i = 1 : numel(angles)

    % Current angle
    a = angles(i);

    % Rotation  matrix by alpha around the normal vector
    R = IceImarisConnector.mapAxisAngleToRotationMatrix(ref_axis, a);

    % Rotate
    currSpotCoords = spotCoords * R(1:3, 1:3);

    % Calculate the distance
    D = pdist2(currSpotCoords, previousSpotCoords, 'squaredeuclidean');

    % Try to be a bit robust against outliers
    allD(i) = median(min(D, [], 2));

end

% Find the smallest distance and return the corresponding angle
angle = angles(allD == min(allD));
angle = angle(1);
