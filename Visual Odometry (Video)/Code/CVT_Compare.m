function [R,t] = CVT_Compare(BW_old, BW, points1, points2, cameraParams)

%% Find Valid Features
[features1,valid_points1] = extractFeatures(BW_old, points1, 'Upright', true);
[features2,valid_points2] = extractFeatures(BW, points2, 'Upright', true);
%% Find Matched Points
[indexPairs, ~] = matchFeatures(features1,features2);
matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);

% if at least 16 points were found to corespond...
if matchedPoints2.Count >= 16
    %% Find the best fundamental matrix
    [F_matlab,inlierIndex] = estimateFundamentalMatrix(matchedPoints1,matchedPoints2);

    %% extract inliers
    m1X = matchedPoints1.Location(:,1);
    m1Y = matchedPoints1.Location(:,2);
    inliers1 = [m1X(inlierIndex) m1Y(inlierIndex)];

    m2X = matchedPoints2.Location(:,1);
    m2Y = matchedPoints2.Location(:,2);
    inliers2 = [m2X(inlierIndex) m2Y(inlierIndex)];

    %% Get the Correct Transformation Combo
    [R,t] = relativeCameraPose(F_matlab,cameraParams,inliers1,inliers2);

else
    
    % set to "nothing"
    R = eye(3);
    t = [0,0,0];
    
end

end

