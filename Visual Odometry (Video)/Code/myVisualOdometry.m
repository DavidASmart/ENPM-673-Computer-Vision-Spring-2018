function [R,t] = myVisualOdometry(BW_old,BW, points1, points2, K)

%% Find Valid Features, Find Matched Points, Remove Incorrect Matches
[matchedPoints1,matchedPoints2] = Find_Matched_Points(BW_old, BW, points1, points2);

% if at least 16 points were found to corespond...
if matchedPoints2.Count >= 16
    %% Find the best fundamental matrix using RANSAC
    F = RANSAC(matchedPoints1, matchedPoints2);

    %% Convert to Esential Matrix
    E = K'*F*K;
    % Perform second decompose to get SE = diag(1, 1, 0)
    [U,~,V] = svd(E);
    E = U * [1 0 0;0 1 0;0 0 0] * V';
    E = E / norm(E);

    %% Get 4 Possible Transformation Combos from Essential Matrix
    [R4, t4] = get_4Rt(E);

    %% Get the Correct Transformation Combo
    [R, t] = get_Rt(R4,t4);

else % if less than 16 points were coresponded
    R = eye(3); % set to nothing
    t = zeros(3,1);
end

end

