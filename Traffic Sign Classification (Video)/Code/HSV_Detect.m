%% ************************************************************************
function [Blue_BB, Red_BB, Im_BB] = HSV_Detect(Im, debugging)
%% 2.1.1 Thresholding in HSV Color Space
% reads in image, returns bounding boxes of possible signs

    Im_BB = Im; % initialize

    %% Use Gaussian Filter to denoise the image

    sigma = 4;
    Im_filt = imgaussfilt(Im,sigma);
    
    %% Threshold in HSV color space by Hue and Value
    
    % convert to hsv
    Im_hsv = rgb2hsv(Im_filt);
    
    % divide hsv
    h = Im_hsv(:,:,1); % hue
    s = Im_hsv(:,:,2); % saturation
    v = Im_hsv(:,:,3); % value
    
    % filter based on hue -------------------------------------------------
    %C_mask = roicolor(h,160/360,180/360);  % Cyan:     120~180
    B_mask = roicolor(h,190/360,240/360);   % Blue:     180~240
    %B_mask = logical(B_mask + C_mask);
    
    M_mask = roicolor(h,270/360,300/360);   % Magenta:  240~300 (260-300) (280 - 300)
    R_mask = roicolor(h,300/360,1);         % Red:      300~360
    Y_mask = roicolor(h,0/360,30/360);      % Yellow:   0~60 (0 - 30)
    R_mask = logical(R_mask + Y_mask + M_mask);
    
    % filter based on value -----------------------------------------------
    V_mask = roicolor(v,0.25,0.80);
    B_mask = logical(B_mask.*V_mask);
    R_mask = logical(R_mask.*V_mask);
    
    %% Filter Regions Based on Convex Area & Eccentricity
    
    % filter based on convex area -----------------------------------------
    mA = (size(Im,1)/40)*(size(Im,2)/40); % minimum area % /30
    MA = (size(Im,1)/10)*(size(Im,2)/10);   % maximum area % /10
    B_mask = bwpropfilt(B_mask,'ConvexArea',[mA,MA]);
    R_mask = bwpropfilt(R_mask,'ConvexArea',[mA,MA]);
    
    % filter based on Eccentricity ----------------------------------------
    % Eccentricity is similar to Aspect Ratio...
    mE = 0;      % minimum Eccentricity
    ME = 0.9;   % maximum Eccentricity
    B_mask = bwpropfilt(B_mask,'Eccentricity',[mE,ME]);
    R_mask = bwpropfilt(R_mask,'Eccentricity',[mE,ME]);
    
    %% Extract the bounding box's
    % Make sure the bounding box covers the entire sign,
    
    % initialize
    Blue_BB = [];
    Red_BB = [];
    
    % get bounding box [x, y, w, h] & axes
    B_BB = regionprops(B_mask,'BoundingBox','MajorAxisLength','MinorAxisLength');
    R_BB = regionprops(R_mask,'BoundingBox','MajorAxisLength','MinorAxisLength');
    
    % Blue ----------------------------------------------------------------
    if ~isempty(B_BB)
        for i = 1:size(B_BB,1) % for each blue region

            % get the bouding box cordinates
            B_bb = B_BB(i).BoundingBox;
            B_l = round((B_BB(i).MajorAxisLength+B_BB(i).MinorAxisLength)*0.05);

            Blue_BB(i,:) = [max(1,B_bb(1)-B_l), max(1,B_bb(2)-B_l),...
                min(size(Im,2)-max(1,B_bb(1)-B_l),B_bb(3)+B_l*2),...
                min(size(Im,1)-max(1,B_bb(2)-B_l),B_bb(4)+B_l*2)];
        end
    end
    
    % Red -----------------------------------------------------------------
    if ~isempty(R_BB)
        for i = 1:size(R_BB,1) % for each blue region

            % get the bouding box cordinates
            R_bb = R_BB(i).BoundingBox;
            R_l = round((R_BB(i).MajorAxisLength+R_BB(i).MinorAxisLength)*0.05);

            Red_BB(i,:) = [max(1,R_bb(1)-R_l), max(1,R_bb(2)-R_l),...
                min(size(Im,2)-max(1,R_bb(1)-R_l),R_bb(3)+R_l*2),...
                min(size(Im,1)-max(1,R_bb(2)-R_l),R_bb(4)+R_l*2)];
        end
    end
    
    
    
    %% display
    if debugging
        % Blue ------------------------------------------------------------
        if ~isempty(Blue_BB)
            for i = 1:size(Blue_BB,1) % for each blue region

                % make shape from cordinates
                xiB = [Blue_BB(i,1), Blue_BB(i,1)+Blue_BB(i,3), Blue_BB(i,1)+Blue_BB(i,3), Blue_BB(i,1)];
                yiB = [Blue_BB(i,2), Blue_BB(i,2), Blue_BB(i,2)+Blue_BB(i,4), Blue_BB(i,2)+Blue_BB(i,4)];

                % make mask from shape
                B_poly = poly2mask(xiB,yiB,size(B_mask,1),size(B_mask,2));

                % get outline of shape
                B_p = bwperim(B_poly);

                % dilate outline for display purposes
                r = 3; SE = strel('sphere',r);
                B_p = imdilate(B_p,SE);

                % overlay on cropped image
                Im_BB = imoverlay(Im_BB,B_p,'b');
            end
        end
        
        % RED -------------------------------------------------------------
        if ~isempty(Red_BB)
            for i = 1:size(Red_BB,1) % for each red region

                % make shape from cordinates
                xiR = [Red_BB(i,1), Red_BB(i,1)+Red_BB(i,3), Red_BB(i,1)+Red_BB(i,3), Red_BB(i,1)];
                yiR = [Red_BB(i,2), Red_BB(i,2), Red_BB(i,2)+Red_BB(i,4), Red_BB(i,2)+Red_BB(i,4)];

                % make mask from shape
                R_poly = poly2mask(xiR,yiR,size(R_mask,1),size(R_mask,2));

                % get outline of shape
                R_p = bwperim(R_poly);

                % dilate outline for display purposes
                r = 3; SE = strel('sphere',r);
                R_p = imdilate(R_p,SE);

                % overlay on cropped image
                Im_BB = imoverlay(Im_BB,R_p,'r');

            end
        end
        
        % display ---------------------------------------------------------
        figure(1); imshow(Im_BB); title('HSV Detection');
    end
end
%% ************************************************************************
