function [Im_C] = HOG_SVM_Classification(Im, Blue_BB, Red_BB, classifierB, classifierR, imagefolder)
%% Classification
% looks at regions detected by MSER and doed a sliding window HOG feature
% extraction and performs svm classification, returns  an image with classification annotations

Im_C = Im; % initialize

%% Blue -------------------------------------------------------------------
labelsB = [{'00035'};{'00038'};{'00045'};{'Background'}];
if ~isempty(Blue_BB)
    for i = 1:size(Blue_BB,1) % for each blue region

        % get the bouding box cordinates
        B_bb = round(Blue_BB(i,:)); % [x,y,w,h]
        szbb = [B_bb(3), B_bb(4)]; % size of bounding box [w,h]

        % initialize
        prediction = 1; 
        score = -1;
        w_bb = B_bb;
        
        for scale = 1:3 % for multiple size windows

            % get window size
            if scale == 1
                szsqw = round(szbb*0.6); % 60-% 
            elseif scale == 2
                szsqw = round(szbb*0.8); % 80-%
            else % if scale == 3
                szsqw = szbb; % 100-%
            end

            % sweep the square over the bounding box 
            for y = B_bb(2):B_bb(2)+B_bb(4)-szsqw(2) % vert
                for x = B_bb(1):B_bb(1)+B_bb(3)-szsqw(1) % horiz
                    
                    xx = min(x+szsqw(1),size(Im,2));
                    yy = min(y+szsqw(2),size(Im,1));
                    
                    % extract that window from the image
                    Im_ex = Im(y:yy, x:xx, 1:3);

                    % contrast normalization
                    %Im_ex = imadjust(Im_ex,stretchlim(Im_ex),[]); 
                    
                    % resize it to 64×64
                    Im_ex = imresize(Im_ex,[64 64]); 

                    % extract HOG features and predict the sign using the trained SVM model. 
                    cellSize = [2 2]; Features = extractHOGFeatures(Im_ex, 'CellSize', cellSize); 

                    % use trained svm to predict the classification
                    [p, s] = predict(classifierB, Features);
                    
                    % if better prediction than previous best
                    if (double(p) ~= size(labelsB,1)) && max(s) > score
                        % save prediction & bounding box
                        prediction = p;
                        score = max(s);
                        w_bb = round([max(1,x-szsqw(1)*0.2), max(1,y-szsqw(2)*0.2), ...
                            min(size(Im,2)-max(1,x-szsqw(1)*0.2),szsqw(1)*1.2),...
                            min(size(Im,1)-max(1,y-szsqw(2)*0.2),szsqw(2)*1.2)]); % add 20-%
                    end

                end
            end
        end

        if score > -0.25
            % Highlight the bounding box of the traffic sign 
            Im_C = insertObjectAnnotation(Im_C,'rectangle',w_bb,labelsB{double(prediction)},...
                'Color','blue','LineWidth',6,'TextBoxOpacity',1,...
                'TextColor','white','FontSize',20);


            % Paste the appropriate sign beside it
            % get appropriate sign
            Im_s = imread(fullfile(imagefolder,strcat(labelsB{double(prediction)},'.jpg')));
            
            % resize it to be aprox. the same as the bounding box
            szs = 64;
            Im_s = imresize(Im_s,[szs szs]);
            
            if w_bb(1)+w_bb(3)+szs <= size(Im,2)
                % insert it into the picture at top right side of bounding box
                Im_C(w_bb(2):w_bb(2)+szs-1,w_bb(1)+w_bb(3)+1:w_bb(1)+w_bb(3)+szs,:) = Im_s;
            end
        end
    end
end

%% RED --------------------------------------------------------------------
labelsR = [{'00001'};{'00014'};{'00017'};{'00019'};{'00021'};{'Background'}];
if ~isempty(Red_BB)
    for i = 1:size(Red_BB,1) % for each blue region

        % get the bouding box cordinates
        R_bb = round(Red_BB(i,:)); % [x,y,w,h]
        szbb = [R_bb(3), R_bb(4)]; % size of bounding box [w,h]

        % initialize
        prediction = 1; 
        score = -1;
        w_bb = R_bb;
        
        for scale = 2:3 % for multiple size windows

            % get window size
            if scale == 1
                szsqw = round(szbb*0.6); % 60-% 
            elseif scale == 2
                szsqw = round(szbb*0.8); % 80-%
            else % if scale == 3
                szsqw = szbb; % 100-%
            end

            % sweep the square over the bounding box 
            for y = R_bb(2):R_bb(2)+R_bb(4)-szsqw(2) % vert
                for x = R_bb(1):R_bb(1)+R_bb(3)-szsqw(1) % horiz
                    
                    xx = min(x+szsqw(1),size(Im,2));
                    yy = min(y+szsqw(2),size(Im,1));
                    
                    % extract that window from the image
                    Im_ex = Im(y:yy, x:xx, 1:3);

                    % contrast normalization
                    %Im_ex = imadjust(Im_ex,stretchlim(Im_ex),[]); 
                    
                    % resize it to 64×64
                    Im_ex = imresize(Im_ex,[64 64]); 

                    % extract HOG features and predict the sign using the trained SVM model. 
                    cellSize = [2 2]; Features = extractHOGFeatures(Im_ex, 'CellSize', cellSize); 

                    % use trained svm to predict the classification
                    [p, s] = predict(classifierR, Features);

                    % if better prediction than previous best
                    if (double(p) ~= size(labelsR,1)) && max(s) > score
                        % save prediction & bounding box
                        prediction = p;
                        score = max(s);
                        w_bb = round([max(1,x-szsqw(1)*0.2), max(1,y-szsqw(2)*0.2), ...
                            min(size(Im,2)-max(1,x-szsqw(1)*0.2),szsqw(1)*1.2),...
                            min(size(Im,1)-max(1,y-szsqw(2)*0.2), szsqw(2)*1.2)]); % add 20-%
                    end
                end
            end
        end

        if score > -0.25
            % Highlight the bounding box of the traffic sign 
            Im_C = insertObjectAnnotation(Im_C,'rectangle',w_bb,labelsR{double(prediction)},...
                'Color','red','LineWidth',6,'TextBoxOpacity',1,...
                'TextColor','white','FontSize',20);

            % Paste the appropriate sign beside it
            % get appropriate sign
            Im_s = imread(fullfile(imagefolder,strcat(labelsR{double(prediction)},'.jpg')));
            
            % resize it to be aprox. the same as the bounding box
            szs = 64;
            Im_s = imresize(Im_s,[szs szs]);
            
            if w_bb(1)+w_bb(3)+szs <= size(Im,2)
                % insert it into the picture at top right side of bounding box
                Im_C(w_bb(2):w_bb(2)+szs-1,w_bb(1)+w_bb(3)+1:w_bb(1)+w_bb(3)+szs,:) = Im_s;
            end
        end
    end
end
    
end
%% ************************************************************************
