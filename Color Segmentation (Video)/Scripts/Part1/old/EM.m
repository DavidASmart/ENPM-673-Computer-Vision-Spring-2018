%% EM.m
function [mu, covar] = EM(N, data, plot_path)
%% PART 1.
% INPUT: 
%       N - number of Gaussians
%       data - Dxn array of data (D is dimensionality)
%       plot path - where to save the resulting plot (relative to current folder)
% OUTPUT: 
%       mu - the array DxN mean
%       covar - DxDxN covariance matrix of model parameters 
% Also plots the computed Gaussians and saves it in plot_path 
%
% **BUILT IN: gm = fitgmdist(X,k) **
% X : n-by-p data matrix X (p is dimensionality)
% k : no. of clusters
% gm : structure which is put into ...
% **BUILT IN: idx = cluster(gm,X)
% idx : n-by-1 vector cluster indices of each X(i).

%% ------------------------------------------------------------------------
% Choose initial parameter values.

    
    tol = 1; % negligable change in means
    maxi = 1000; % max no. of iterations
    
    n = size(data, 2); % number of data points
    D = size(data, 1); % dimensionality of data
    
%     % initialze GMM cluster assignment by using K-means
    [idx, C] = K_MEANS(N, data, '');
%     for some undiscovered reason the k-means does not always find enough
%     clusters and so that caused this function to fail...
    
    while max(isnan(C))
        % initialize with random start
        idx = randi(N,1,size(data,2)); % random assignment
        for j = 1:N % for all clusters
            xi = data(:,idx == j); % get the datapoints which belong to this cluster
            C(:,j) = sum(xi) / size(xi,2); % calculate centroid
        end
    end

    mu = C; % initial means are the centroids from K-means
    for j = 1:N % for each cluster
        % Use the co-variance of the clusters from K-means
        covar{j} = cov(data(:,idx == j)); % sigma is a structure for ease of storage
        
        % current proportion of data in each cluster
        phi(j) = length(data(:,idx == j))/n;
    end

%% ------------------------------------------------------------------------
    % Expectation Maximization
    for iter = 1:maxi % until maxi interations have been completed

        %% Step-E : Expectation -------------------------------------------

        % initialize
        pdf = zeros(N,n);
        pdf_w = zeros(N,n);
        
        % For each gaussian-cluster
        for j = 1:N
            % calculate the probability of each data point fitting that gaussian distribution
            pdf(j,:) = 1/sqrt(((2*pi)^D)*det(covar{j})) * ...
                exp((-1/2)*(data - mu(:,j)).*inv(covar{j}).*(data - mu(:,j)));
            
            pdf_w(j,:) = pdf(j,:).*phi(j); % wheigted probability
        end
        
        % normalize the weighted probabilities
        W = pdf_w./sum(pdf_w);


        %% Step-M : Maximization ------------------------------------------

        mu_old = mu; % update   

        for j = 1:N % For each of the clusters...

            % prior probability (i.e. "proportion of data in each cluster")
            phi(j) = mean(W(j,:));
            
            % mean = weighted average of all data points
            mu(:,j) = sum(W(j,:).*data)./sum(W(j,:));
            
            % covariance = weighted covarnaice of all data points
            % calculate the contribution of each datapoint to the covariance matrix
            sigma_sum = 0; % initialize
            for i = 1:n
                sigma_sum = sigma_sum + (W(j,i).*((data(:,i) - mu(:,j))'*(data(:,i) - mu(:,j))));
            end
            % normalize
            covar{j} = sigma_sum./sum(W(j,:)) + eye(D)*(1e-6);
                        
        end

        % check for convergence
        mudif = max(100*abs((mu - mu_old)./mu_old));
        if ~isfinite(mudif) || isnan(mudif) % just in case
            mudif = 1000;
        end
        if mudif < tol % converged?
            break
        end
    end
    
    [~,idx] = max(W); % get the 'hard' assignments from the 'soft' assignments
    
%% ------------------------------------------------------------------------

    if ~isempty(plot_path) % if a folder was specified
    
        figure; title(strcat('EM',num2str(size(data,1)),'D',num2str(N),'N')); hold on;
        % depending on the dimensionality of the data, the output should be differnt
        if size(data,1) == 1 % 1D linear
            for j = 1:N % for each cluster
                xi = data(:,idx == j); % get the data belonging to that cluster
                plot(xi,zeros(1,size(xi,2)),'.'); % plot data in each culster as a distinct color
                plot(mu(1,j),0,'ko','MarkerFaceColor','k'); % plot cluster centers
            end
            
        elseif size(data,1) == 2 % 2D space
            for j = 1:N % for each cluster
                xi = data(:,idx == j); % get the data belonging to that cluster
                plot(xi(1,:),xi(2,:),'.'); % plot data in each culster as a distinct color
                plot(mu(1,j),mu(2,j),'ko','MarkerFaceColor','k'); % plot cluster centers
            end
            
        elseif size(data,1) == 3 % 3D volume
             for j = 1:N % for each cluster
                xi = data(:,idx == j); % get the data belonging to that cluster
                plot3(xi(1,:),xi(2,:),xi(3,:),'.'); % plot data in each culster as a distinct color
                plot3(mu(1,j),mu(2,j),mu(3,j),'ko','MarkerFaceColor','k'); % plot cluster centers
             end
             
        elseif size(data,1) == 5 % 2D + RGB = color image
            for j = 1:N % for each cluster
                xi = data(:,idx == j); % get the data belonging to that cluster
                plot(xi(1,:),xi(2,:),'s','MarkerFaceColor',[mu(3,j) mu(4,j) mu(5,j)]); % plot data as its cluster color
                plot(mu(1,j),mu(2,j),'ko','MarkerSize',10); % plot cluster centers
            end
            
        end % ...I don't know what to do if it is a diffent dimensionality...
        hold off

        % save
        currentfolder = pwd; % save the current folder so that it can be returned to after saving figure
        cd(plot_path); % change to folder specified
        saveas(gcf,strcat('EM',num2str(size(data,1)),'D',num2str(N),'N.jpg'));
        cd(currentfolder); % return to the scripts folder
    end

end

%% ------------------------------------------------------------------------
