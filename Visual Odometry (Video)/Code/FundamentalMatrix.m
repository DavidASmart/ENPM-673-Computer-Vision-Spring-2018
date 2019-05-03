function F = FundamentalMatrix(x1, x2)

%% First Transform
cx1_x = mean(x1(:,1)); % Centroid
cx1_y = mean(x1(:,2));
nx1_x = x1(:,1) - cx1_x;% Shift origin to centroid.
nx1_y = x1(:,2) - cx1_y;
avg_d = sqrt(sum(nx1_x.^2  + nx1_y.^2)) / length(x1);
s1 = sqrt(2) / avg_d; % scalling factor

x1(:,1) = s1 * nx1_x; % normalize points
x1(:,2) = s1 * nx1_y;

T_1 = [s1,0,-s1*cx1_x;
    0,s1,-s1*cx1_y;
    0,0,1]; % transform matrix

%% Second Transform
cx2_x = mean(x2(:,1)); % Centroid
cx2_y = mean(x2(:,2));
nx2_x = x2(:,1) - cx2_x; % Shift origin to centroid.
nx2_y = x2(:,2) - cx2_y;
avg_d = sqrt(sum(nx2_x.^2  + nx2_y.^2)) / length(x2);
s2 = sqrt(2) / avg_d; % scalling factor

x2(:,1) = s2 * nx2_x; % normalize points
x2(:,2) = s2 * nx2_y;

T_2 = [s2 0 -s2*cx2_x;
       0 s2 -s2*cx2_y;
       0 0 1]; % transform matrix

%% create matrix Af = 0
A = [x1(:,1).*x2(:,1), x1(:,1).*x2(:,2), x1(:,1), x1(:,2).*x2(:,1), x1(:,2).*x2(:,2), x1(:,2), x2(:,1), x2(:,2), ones(length(x1),1)];

%% solve for f
[~,~,V] = svd(A);
f = V(:,end);
F = reshape(f,[3,3]);

% normalize F
F = F / norm(F);

%% enforce rank 2
[U,S,V] = svd(F);
S(3,3) = 0;
F = U*S*V';

%% de-normalize
F = T_2'*F*T_1;
end

