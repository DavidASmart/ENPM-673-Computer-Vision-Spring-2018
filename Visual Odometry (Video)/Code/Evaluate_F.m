function [score, inliers] = Evaluate_F(F, matchedpoints1, matchedpoints2)

    %% Evaluate the candidate Fundamental Matrix
    % initialize
    inliers = [];
    
    % for all the matched points
    for i = 1:size(matchedpoints1,1) 
        
       % convert points to homogeneous cordinates
       x1 = [matchedpoints1(i,:)';1];
       x2 = [matchedpoints2(i,:)';1]; 
       
       % calculate the epipolar line from image 1
       el = F * x1;
       
       % calculate distance from the epipolar line of the points in image 2 
       err = abs(x2' * F * x1) / (sqrt(el(1)^2 + el(2)^2)); 
       
       if err < 3 % if less than 3 pixels off
           inliers = [inliers; i]; % add to list of inliers
       end
       
    end
    
    % calculate the fraction of inliers
    score = size(inliers,1) / size(matchedpoints1,1);

end