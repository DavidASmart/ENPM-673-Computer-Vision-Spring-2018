function [Blue_BB, Red_BB, Im_BB] = MSER_Detect(Im, B_BB, R_BB, debugging)
%% 2.1.2 Maximally Stable Extremal Region (MSER) Detection

    Im_BB = Im; % initialize

    %% Denoise the image
    
%     sigma = 0.2;
%     alpha = 3;
%     Im_filt = im2double(locallapfilt(Im,sigma,alpha)); % coolest filter ever!
% 
%     dos = (0.01*diff(getrangefromclass(Im)).^2)*3;
%     sigma = 4;
%     Im_filt = im2double(imbilatfilt(Im,dos,sigma));

    sigma = 4;
    Im_filt = im2double(imgaussfilt(Im,sigma));
    
    %% Contrast Normalization
    Im_cn = imadjust(Im_filt,stretchlim(Im_filt),[]);

    %% Divide RGB Components
    R = Im_cn(:,:,1);
    G = Im_cn(:,:,2);
    B = Im_cn(:,:,3);
    
    %% Normalize Intensity
    
    % individual adjustement (simple method) ------------------------------
    R2 = R./(R+G+B);
    G2 = G./(R+G+B);
    B2 = B./(R+G+B);
    R = R2; G = G2; B = B2;
    
    % red and blue enhancement --------------------------------------------
    R2 = max(0, min(R-B, R-G)./(R+G+B));
    B2 = max(0, (B-R)./(R+G+B));
    R = R2; B = B2;
    
    %% Extract & filter MSER regions (from each regions detected by HSV_detect)
    
    % initialize incase B_BB and R_BB are empty
    ccB = [];
    ccR = [];
    
    % BLUE ----------------------------------------------------------------
    if ~isempty(B_BB)
        for i = 1:size(B_BB,1) % for each blue region
            
            % get the bouding box cordinates
            B_bb = B_BB(i,:);

            if i > 1
                % MSER detection
                [RegionsBt,ccBt] = detectMSERFeatures(B,'ThresholdDelta',4,'MaxAreaVariation',0.175,...
                    'ROI', B_bb); %[1, 1, size(Im,2), size(Im,1)]); % def: 2, 0.25; **prev: 3,0.175**
                
                % add to others
                idxS = RegionsB.Count+1;
                idxF = RegionsB.Count+RegionsBt.Count;
                RegionsB.Count = idxF;
                RegionsB.Location(idxS:idxF,:) = RegionsBt.Location;
                RegionsB.Axes(idxS:idxF,:) = RegionsBt.Axes;
                RegionsB.Orientation(idxS:idxF,:) = RegionsBt.Orientation;
                if size(RegionsBt.PixelList,2) == 2
                    RegionsB.PixelList(idxS:idxF,:) = {RegionsBt.PixelList};
                else
                    RegionsB.PixelList(idxS:idxF,:) = RegionsBt.PixelList;
                end
                
                ccB.NumObjects = ccB.NumObjects + ccBt.NumObjects;
                ccB.PixelIdxList(:,idxS:idxF) = ccBt.PixelIdxList;
                
            else
                % MSER detection
                [RegionsBt,ccBt] = detectMSERFeatures(B,'ThresholdDelta',4,'MaxAreaVariation',0.20,...
                    'ROI', B_bb);%[1, 1, size(Im,2), size(Im,1)]); % def: 2, 0.1; **prev: 3, 0.175**
                
                % create non-"MSERRegions object" so that it can be eddited
                RegionsB.Count = RegionsBt.Count;
                RegionsB.Location = RegionsBt.Location;
                RegionsB.Axes = RegionsBt.Axes;
                RegionsB.Orientation = RegionsBt.Orientation;
                if size(RegionsBt.PixelList,2) == 2
                    RegionsB.PixelList = {RegionsBt.PixelList};
                else
                    RegionsB.PixelList = RegionsBt.PixelList;
                end
                
                ccB.Connectivity = 8;
                ccB.ImageSize = size(B);
                ccB.NumObjects = ccBt.NumObjects;
                ccB.PixelIdxList = ccBt.PixelIdxList;
                
            end
        end
        % remake a true "MSERRegions object"
        RegionsB = MSERRegions(RegionsB.PixelList);
        
        % filter based on convex area -----------------------------------------
        mA = (size(Im,1)/40)*(size(Im,2)/40); % minimum area
        MA = (size(Im,1)/10)*(size(Im,2)/10);   % maximum area
        statsB = regionprops('table',ccB,'ConvexArea');
        CAIdxB = statsB.ConvexArea > mA;
        RegionsB = RegionsB(CAIdxB); 
        ccB.PixelIdxList = ccB.PixelIdxList(CAIdxB);
        ccB.NumObjects = length(ccB.PixelIdxList);
        statsB = regionprops('table',ccB,'ConvexArea');
        CAIdxB = statsB.ConvexArea < MA;
        RegionsB = RegionsB(CAIdxB);
        ccB.PixelIdxList = ccB.PixelIdxList(CAIdxB); 
        ccB.NumObjects = length(ccB.PixelIdxList);
        
        % filter based on "eccentricity"  -------------------------------------
        statsB = regionprops('table',ccB,'Eccentricity');
        eccIdxB = statsB.Eccentricity < 0.85;
        RegionsB = RegionsB(eccIdxB);
        ccB.PixelIdxList = ccB.PixelIdxList(eccIdxB); 
        ccB.NumObjects = length(ccB.PixelIdxList);
        
        % filter based on possition  ------------------------------------------
        % [top-left top-right right-top right-bottom bottom-right bottom-left left-bottom left-top]

        % limits
        %Left = size(Im,2)*0.025;
        %Right = size(Im,2)*(1-0.025);
        %Top = size(Im,1)*0.025;
        Bottom = size(Im,1)*(1-0.35); % 30
        
        statsB = regionprops('table',ccB,'Extrema'); statsB = statsB.Extrema;
        ExIdxB = ones(ccB.NumObjects,1);  % initilize list of those to be removed
        for i = 1:length(statsB) % for every MSER region
            sB = statsB{i}; % get current MSER extrema
            % check each extrema
            %if sB(1,2) < Top || sB(2,2) < Top || sB(3,2) < Top || sB(8,2) < Top % top-left or top-right or right-top or left-top = far-top
                %ExIdxB(i) = 0;
            %elseif sB(2,1) > Right || sB(3,1) > Right || sB(4,1) > Right || sB(5,1) > Right % top-right or right-top or right-bottom or bottom-right = far-right
                %ExIdxB(i) = 0;
            %   else
            if sB(4,2) > Bottom || sB(5,2) > Bottom || sB(6,2) > Bottom || sB(7,2) > Bottom % right-bottom or bottom-right or bottom-left or left-bottom = far-bottom
                ExIdxB(i) = 0;
            %elseif sB(1,1) < Left || sB(6,1) < Left || sB(7,1) < Left || sB(8,1) < Left % top-left or bottom-left or left-bottom or left-top = far-left
                %ExIdxB(i) = 0;
            end
        end
        ExIdxB = logical(ExIdxB);
        RegionsB = RegionsB(ExIdxB);
        ccB.PixelIdxList = ccB.PixelIdxList(ExIdxB); 
        ccB.NumObjects = length(ccB.PixelIdxList);
        
        % keep only the largest in a given region -----------------------------
        idxB = ones(ccB.NumObjects,1); % initilize list of those to be removed
        for r = 1:length(RegionsB) % for every MSER region
            for rr = 1:length(RegionsB) % check against every other MSER region
                if r ~= rr % don't compare against self

                    % measure distance from currect MSER center to other MSER features
                    dx = RegionsB.Location(r,1) - RegionsB.Location(rr,1);
                    dy = RegionsB.Location(r,2) - RegionsB.Location(rr,2);
                    d = sqrt(dx^2+dy^2);

                    % if within current MSER's major-axis length
                    if d < RegionsB.Axes(r,1)
                        % if current is larger
                        if RegionsB.Axes(r,1) > RegionsB.Axes(rr,1)
                            % add other MSER feature to list to be removed
                            idxB(rr) = 0;
                        end
                    end
                end
            end
        end
        idxB = logical(idxB);
        RegionsB = RegionsB(idxB);
        ccB.PixelIdxList = ccB.PixelIdxList(idxB); 
        ccB.NumObjects = length(ccB.PixelIdxList);
        
    end
    
    % RED -----------------------------------------------------------------
    if ~isempty(R_BB)
        for i = 1:size(R_BB,1) % for each blue region
            
            % get the bouding box cordinates
            R_bb = R_BB(i,:);

            if i > 1
                % MSER detection
                [RegionsRt,ccRt] = detectMSERFeatures(R,'ThresholdDelta',2,'MaxAreaVariation',0.25,...
                    'ROI', R_bb);%[1, 1, size(Im,2), size(Im,1)]); % def: 2, 0.25 
                
                % add to others
                idxS = RegionsR.Count+1;
                idxF = RegionsR.Count+RegionsRt.Count;
                RegionsR.Count = idxF;
                RegionsR.Location(idxS:idxF,:) = RegionsRt.Location;
                RegionsR.Axes(idxS:idxF,:) = RegionsRt.Axes;
                RegionsR.Orientation(idxS:idxF,:) = RegionsRt.Orientation;
                if size(RegionsRt.PixelList,2) == 2
                    RegionsR.PixelList(idxS:idxF,:) = {RegionsRt.PixelList};
                else
                    RegionsR.PixelList(idxS:idxF,:) = RegionsRt.PixelList;
                end
                
                ccR.NumObjects = ccR.NumObjects + ccRt.NumObjects;
                ccR.PixelIdxList(:,idxS:idxF) = ccRt.PixelIdxList;
                
            else
                % MSER detection
                [RegionsRt,ccRt] = detectMSERFeatures(R,'ThresholdDelta',2,'MaxAreaVariation',025,...
                    'ROI', R_bb);%[1, 1, size(Im,2), size(Im,1)]); % def: 2, 0.25 
                
                % create non-"MSERRegions object" so that it can be eddited
                RegionsR.Count = RegionsRt.Count;
                RegionsR.Location = RegionsRt.Location;
                RegionsR.Axes = RegionsRt.Axes;
                RegionsR.Orientation = RegionsRt.Orientation;
                if size(RegionsRt.PixelList,2) == 2
                    RegionsR.PixelList = {RegionsRt.PixelList};
                else
                    RegionsR.PixelList = RegionsRt.PixelList;
                end
                
                ccR.Connectivity = 8;
                ccR.ImageSize = size(R);
                ccR.NumObjects = ccRt.NumObjects;
                ccR.PixelIdxList = ccRt.PixelIdxList;
                
            end
        end
        % remake a true "MSERRegions object"
        RegionsR = MSERRegions(RegionsR.PixelList);
        
        % filter based on convex area -----------------------------------------
        mA = (size(Im,1)/40)*(size(Im,2)/40); % minimum area
        MA = (size(Im,1)/10)*(size(Im,2)/10);   % maximum area
        statsR = regionprops('table',ccR,'ConvexArea');
        CAIdxR = statsR.ConvexArea > mA;
        RegionsR = RegionsR(CAIdxR);
        ccR.PixelIdxList = ccR.PixelIdxList(CAIdxR); 
        ccR.NumObjects = length(ccR.PixelIdxList);
        statsR = regionprops('table',ccR,'ConvexArea');
        CAIdxR = statsR.ConvexArea < MA;
        RegionsR = RegionsR(CAIdxR);
        ccR.PixelIdxList = ccR.PixelIdxList(CAIdxR); 
        ccR.NumObjects = length(ccR.PixelIdxList);
        
        % filter based on "eccentricity"  -------------------------------------
        statsR = regionprops('table',ccR,'Eccentricity');
        eccIdxR = statsR.Eccentricity < 0.85;
        RegionsR = RegionsR(eccIdxR);
        ccR.PixelIdxList = ccR.PixelIdxList(eccIdxR); 
        ccR.NumObjects = length(ccR.PixelIdxList);
        
        % filter based on possition  ------------------------------------------
        % [top-left top-right right-top right-bottom bottom-right bottom-left left-bottom left-top]

        % limits
        %Left = size(Im,2)*0.05;
        %Right = size(Im,2)*(1-0.05);
        %Top = size(Im,1)*0.05;
        Bottom = size(Im,1)*(1-0.05);
        
        statsR = regionprops('table',ccR,'Extrema'); statsR = statsR.Extrema;
        ExIdxR = ones(ccR.NumObjects,1);  % initilize list of those to be removed
        for i = 1:length(statsR) % for every MSER region
            sR = statsR{i}; % get current MSER extrema
            % check each extrema
            %if sR(1,2) < Top || sR(2,2) < Top || sR(3,2) < Top || sR(8,2) < Top % top-left or top-right or right-top or left-top = far-top
                % ExIdxR(i) = 0;
            %elseif sR(2,1) > Right || sR(3,1) > Right || sR(4,1) > Right || sR(5,1) > Right % top-right or right-top or right-bottom or bottom-right = far-right
                %ExIdxR(i) = 0;
            %    else
            if sR(4,2) > Bottom || sR(5,2) > Bottom || sR(6,2) > Bottom || sR(7,2) > Bottom % right-bottom or bottom-right or bottom-left or left-bottom = far-bottom
                ExIdxR(i) = 0;
            %elseif sR(1,1) < Left || sR(6,1) < Left || sR(7,1) < Left || sR(8,1) < Left % top-left or bottom-left or left-bottom or left-top = far-left
                %ExIdxR(i) = 0;
            end
        end
        ExIdxR = logical(ExIdxR);
        RegionsR = RegionsR(ExIdxR);
        ccR.PixelIdxList = ccR.PixelIdxList(ExIdxR); 
        ccR.NumObjects = length(ccR.PixelIdxList);
        
        % keep only the largest in a given region -----------------------------
        idxR = ones(ccR.NumObjects,1); % initilize list of those to be removed
        for r = 1:length(RegionsR) % for every MSER region
            for rr = 1:length(RegionsR) % check against every other MSER region
                if r ~= rr % don't compare against self

                    % measure distance from currect MSER center to other MSER features
                    dx = RegionsR.Location(r,1) - RegionsR.Location(rr,1);
                    dy = RegionsR.Location(r,2) - RegionsR.Location(rr,2);
                    d = sqrt(dx^2+dy^2);

                    % if within current MSER's major-axis length
                    if d < RegionsR.Axes(r,1)
                        % if current is larger
                        if RegionsR.Axes(r,1) > RegionsR.Axes(rr,1)
                            % add other MSER feature to list to be removed
                            idxR(rr) = 0;
                        end
                    end
                end
            end
        end
        idxR = logical(idxR);
        RegionsR = RegionsR(idxR);
        ccR.PixelIdxList = ccR.PixelIdxList(idxR); 
        ccR.NumObjects = length(ccR.PixelIdxList);
        
    end
    
    %% Extract the new bounding box's
    % Make sure the bounding box covers the entire sign, 

    % initialize
    Blue_BB = [];
    Red_BB = [];
    
    % Blue ----------------------------------------------------------------
    if ~isempty(ccB)
        B_BB = regionprops('table',ccB,'BoundingBox','MajorAxisLength','MinorAxisLength');
        if ~isempty(B_BB)
            for i = 1:size(B_BB,1) % for each blue region

               % get the bouding box cordinates
                B_bb = B_BB.BoundingBox(i,1:4); 
                B_l = 0;%round((B_BB.MajorAxisLength(i)+B_BB.MinorAxisLength(i))*0.025);

                Blue_BB(i,:) = [max(1,B_bb(1)-B_l), max(1,B_bb(2)-B_l),...
                    min(size(Im,2)-max(1,B_bb(1)-B_l),B_bb(3)+B_l*2),...
                    min(size(Im,1)-max(1,B_bb(2)-B_l),B_bb(4)+B_l*2)];
            end
        end
    end
    
    % Red -----------------------------------------------------------------
    if ~isempty(ccR)
        R_BB = regionprops('table',ccR,'BoundingBox','MajorAxisLength','MinorAxisLength');
        if ~isempty(R_BB)
            for i = 1:size(R_BB,1) % for each blue region

                % get the bouding box cordinates
                R_bb = R_BB.BoundingBox(i,1:4); 
                R_l = 0;%round((R_BB.MajorAxisLength(i)+R_BB.MinorAxisLength(i))*0.025);

                 Red_BB(i,:) = [max(1,R_bb(1)-R_l), max(1,R_bb(2)-R_l),...
                    min(size(Im,2)-max(1,R_bb(1)-R_l),R_bb(3)+R_l*2),...
                    min(size(Im,1)-max(1,R_bb(2)-R_l),R_bb(4)+R_l*2)];
            end
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
                B_poly = poly2mask(xiB,yiB,size(Im,1),size(Im,2));

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
                R_poly = poly2mask(xiR,yiR,size(Im,1),size(Im,2));

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
        figure(2); imshow(Im_BB); title('MSER Detection');
    end
    
end

