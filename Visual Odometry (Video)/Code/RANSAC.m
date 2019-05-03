function F = RANSAC(matchedpoints1, matchedpoints2)

% initialize
F = zeros(3,3);
n = 0; % iteration counter
N = 2000; % max no. iterations
p = 0.95; % p-percent chance of 100 inliers
npts = matchedpoints2.Count;
old_score = 0;

while n < N % for max number of iterations set
    
    %% Select at random 8 datapoints
    ind = randsample(npts, 8);
    x1 = matchedpoints1(ind);
    x2 = matchedpoints2(ind);
    
    %% Estimate the Fundamental Matrix
    M = FundamentalMatrix(x1.Location, x2.Location);
    
    %% Evaluate the candidate Fundamental Matrix
    [new_score,~] = Evaluate_F(M, matchedpoints1.Location, matchedpoints2.Location);
    
    %% Compare against the Previous Best
    if  new_score > old_score % if better
        F = M; % update the best Fundamental Matrix
        old_score = new_score; % update the percentage of inliers
        N = min(N, log(1-p) / log(1 - new_score^8)); % update the max number of iterations probabilisticaly
    end
    
    n = n + 1; % iterate
end

end
