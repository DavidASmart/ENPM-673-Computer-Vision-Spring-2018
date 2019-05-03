function [matchedPoints1,matchedPoints2] = Find_Matched_Points(BW_old, BW, points1, points2)

%% Find Valid Features
[features1,valid_points1] = extractFeatures(BW_old,points1, 'Upright', true);
[features2,valid_points2] = extractFeatures(BW,points2, 'Upright', true);
%% Find Matched Points
[indexPairs, ~] = matchFeatures(features1,features2);
matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);

%% Remove Incorrect Matches
mp1 = matchedPoints1.Location; % extract possitions
mp2 = matchedPoints2.Location;
% Find Average Distanc
for i = size(mp1,1):-1:1 % for all points
    % calculate distance
    dx = mp1(i,1) - mp2(i,1);
    dy = mp1(i,2) - mp2(i,2);
    d(i) = sqrt(dx^2+dy^2);
end
avgd = mean(d);
% Remove any pair who's distance is twice the average
thresh = avgd*2;
for i = size(mp1,1):-1:1 % for all points
    % calculate distance
    dx = mp1(i,1) - mp2(i,1);
    dy = mp1(i,2) - mp2(i,2);
    d = sqrt(dx^2+dy^2);
    % evalute threshold
    if d > thresh
        % remove
        matchedPoints1(i) = [];
        matchedPoints2(i) = [];
    end
end

end

