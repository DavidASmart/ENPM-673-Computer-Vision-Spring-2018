%% K_MEANS.m
function [idx, C] = K_MEANS(N, data, plot_path)
%% PART 1.
% INPUT:
%       N - number of Gaussians
%       data - Dxn matrix data (D is dimensionaility)
%       plot path - where to save the resulting plot (relative to current folder)
% OUTPUT: 
%       mu - the array NxD mean
%       covar - NxDxD covariance matrix of model parameters 
% Also plots the computed clusters and saves it in plot_path

% **BUILT IN: [idx,C] = kmeans (X,k)**
%       X : n-by-p data matrix X (p is no. of dimensions)
%       k : no. of clusters
%       idx : n-by-1 vector cluster indices of each X(i).
%       C : k-by-p matrix of cluster center locations

% NOTES FROM: http://www.cs.ubbcluj.ro/~csatol/gep_tan/Bishop-CUED-2006.pdf
%       r(n,k) = {0,1};
%       sum(r(n,k) = 1) % each point is assigned to 1 cluster
%   step E 
%       J = sum_n(sum_k(r(n,k)*||x(n)-mu(k)||^2));
%   step M
%       mu(j) = sum_n(r(n,j)*x(n))/sum_n(r(n,j))

    % set up convergance constraints
    max_i = 100; % max no. of iterations
    tol = 1; % negligable change in Centroid locations

    % initialize
    C = initCentroids(data, N);
    while max(max(isnan(C)))
        C = initCentroids(data, N);
    end
    C_old = C; Cdif = 1000;
    i = 1; 

    % Estimation Maximization
    while i < max_i && Cdif > tol % until any convergance criteria is met

        % step E : Assigning data points to clusters
        idx = getClosestCentroids(data, C);

        % step M : Re-computing centroids of clusters
        C = computeCentroids(data, idx, N);

        i = i + 1; % keep track of iterations
        Cdif = max(100*sum(abs((C - C_old)./C_old),1)); % change in centriod possitions
        if ~isfinite(Cdif) || isnan(Cdif)
            Cdif = 1000;
        end
        C_old = C; % update

    end
    
    if ~isempty(plot_path) % if a folder was specified
    
        figure; title(strcat('K-MEANS',num2str(size(data,1)),'D',num2str(N),'N')); hold on;
        % depending on the dimensionality of the data, the output should be differnt
        if size(data,1) == 1 % 1D linear
            for j = 1:N % for each cluster
                xi = data(:,idx == j); % get the data belonging to that cluster
                plot(xi,zeros(1,size(xi,2)),'.'); % plot data in each culster as a distinct color
                plot(C(1,j),0,'ko','MarkerFaceColor','k'); % plot cluster centers
            end
        
        elseif size(data,1) == 2 % 2D space
            for j = 1:N % for each cluster
                xi = data(:,idx == j); % get the data belonging to that cluster
                plot(xi(1,:),xi(2,:),'.'); % plot data in each culster as a distinct color
                plot(C(1,j),C(2,j),'ko','MarkerFaceColor','k'); % plot cluster centers
            end
            
        elseif size(data,1) == 3 % 3D volume
             for j = 1:N % for each cluster
                xi = data(:,idx == j); % get the data belonging to that cluster
                plot3(xi(1,:),xi(2,:),xi(3,:),'.'); % plot data in each culster as a distinct color
                plot3(C(1,j),C(2,j),C(3,j),'ko','MarkerFaceColor','k'); % plot cluster centers
                view(45,45);
             end
             
        elseif size(data,1) == 5 % 2D + RGB = color image
            for j = 1:N % for each cluster
                if ~isnan(C(:,j))
                    xi = data(:,idx == j); % get the data belonging to that cluster
                    scatter(xi(1,:),xi(2,:),'s','MarkerFaceColor',[C(3,j) C(4,j) C(5,j)],'MarkerEdgeAlpha',0); % plot data as its cluster color
                end
            end
            
        end
        hold off

        % save
        currentfolder = pwd; % save the current folder so that it can be returned to after saving figure
        cd(plot_path); % change to folder specified
        saveas(gcf,strcat('K_MEANS',num2str(size(data,1)),'D',num2str(N),'N.jpg'));
        cd(currentfolder); % return to the scripts folder
    end

end

%% ------------------------------------------------------------------------
function C = initCentroids(data, N)
% initialize

    % C = zeros(N,size(data,1)); % initialize
    idx = randi(N,1,size(data,2)); % random assignment
    C = computeCentroids(data, idx, N); % compute intial centroids
    
end

%% ------------------------------------------------------------------------
function idx = getClosestCentroids(data, C)
% step E : Assigning data points to clusters  

    N = size(C, 2); % number of clusters
    idx = zeros(1,size(data,2)); % initialization

    for i = 1:size(data,2) % for all data points
        k = 1; % initialize to first cluster centroid
        min_dist = sqrt(sum((data(:,i) - C(:,1)).^ 2)); % initialize minimum distance as distance to fist cluster centroid
        for j = 2:N % for the rest of the clusters
            dist = sqrt(sum((data(:,i) - C(:,j)).^ 2)); % claculate the distance to the centroid
            if(dist < min_dist) % if the distance is less
                min_dist = dist; % update the minimum distance & ...
                k = j; % .. update the centroid assignmnet
            end
        end
        idx(i) = k; % make the assignment
    end
    
end

%% ------------------------------------------------------------------------
function C = computeCentroids(data, idx, N)
% step M : Re-computing centroids of clusters

	D = size(data,1); % get dimensionality
	C = zeros(D,N); % initialize
  
    for j = 1:N % for all clusters
        xi = data(:,idx == j); % get the datapoints which belong to this cluster
        C(:,j) = sum(xi,2) / size(xi,2); % calculate centroid
    end
    
end
%% ------------------------------------------------------------------------
