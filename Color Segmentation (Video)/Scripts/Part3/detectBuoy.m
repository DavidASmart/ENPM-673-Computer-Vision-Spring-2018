%% PART 3.
% Buoy  Detection (40) 
%
% 1.
% Compute the color-segmented binary image for each buoy 
% and compute the corresponding contours and center.
%
% 2.
% Draw the bounding contours around the detected buoy 
% (color of the contour should be the color of buoy detected). 
 
function [] = detectBuoy(Frame, muY, muR, muG, covarY, covarR, covarG, plot_path)

    %% --------------------------------------------------------------------
    
    current = pwd; % current folder path
    
    % read in images from
    TrainingSetFolder = '../../Images/TrainingSet/Frames';
    TestSetFolder = '../../Images/TestSet/Frames';
    
    % try both image folders to find the correct image
    cd(current); cd(TrainingSetFolder); % change to input folder
    ImageInfo = dir ('*.jpg'); % gather image file names from folder
    for i = 1:length(ImageInfo)
        c(i) = strcmp(ImageInfo(i).name,strcat(num2str(Frame),'.jpg'));
    end
    [found,id] = max(c); % in trainingset?

    if ~found % in testset
        cd(current); cd(TestSetFolder); % change to input folder
        ImageInfo = dir ('*.jpg'); % gather image file names from folder
        for i = 1:length(ImageInfo)
            c(i) = strcmp(ImageInfo(i).name,strcat(num2str(Frame),'.jpg'));
        end
        [~,id] = max(c); % in trainingset?
    end
    
    Im = imread(ImageInfo(id).name); % read in image file
    % **filtering**
    Im2 = imgaussfilt(Im);
    % Im2 = medfilt2(Im);
    I = double(Im2);
    
    %% --------------------------------------------------------------------
    
    % initialize masks
    my = zeros(size(I,1),size(I,2)); % yellow bouy mask
    mr = zeros(size(I,1),size(I,2)); % red bouy mask
    mg = zeros(size(I,1),size(I,2)); % green bouy mask

    % probability of fitting to 3D model which passes through the filter
    tolY = 0.4; 
    tolR = 0.4; 
    tolG = 0.5;
    
    % threshold all pixels based on 3D-gaussian models
    
    % YELLOW BUOY
    for i = 1:size(I,2) % all x
        for j = 1:size(I,1) % all y
            R = I(j,i,1); G = I(j,i,2); B = I(j,i,3); X = [R;G;B];
            for k = 1:size(muY,2) % for each gaussian model
                pm(k) = mvnpdf(muY(:,k),muY(:,k),covarY{k}); % max pdf
                p(k) = mvnpdf(X,muY(:,k),covarY{k}); % pdf
                p(k) = p(k)/pm(k); % probability of fitting
            end
            
            if (max(p) > tolY) % if it passes any of the gaussians
                my(j,i) = 1;
            else
                my(j,i) = 0;
            end

        end
    end
    
    % RED BUOY
    for i = 1:size(I,2) % all x
        for j = 1:size(I,1) % all y
            R = I(j,i,1); G = I(j,i,2); B = I(j,i,3); X = [R;G;B];
            for k = 1:size(muR,2) % for each gaussian model
                pm(k) = mvnpdf(muR(:,k),muR(:,k),covarR{k}); % max pdf
                p(k) = mvnpdf(X,muR(:,k),covarR{k}); % pdf
                p(k) = p(k)/pm(k); % probability of fitting
            end
            
            if (max(p) > tolR) % if it passes any of the gaussians
                mr(j,i) = 1;
            else
                mr(j,i) = 0;
            end

        end
    end

    % GREEN BUOY
    for i = 1:size(I,2) % all x
        for j = 1:size(I,1) % all y
            R = I(j,i,1); G = I(j,i,2); B = I(j,i,3); X = [R;G;B];
            for k = 1:size(muG,2) % for each gaussian model
                pm(k) = mvnpdf(muG(:,k),muG(:,k),covarG{k}); % max pdf
                p(k) = mvnpdf(X,muG(:,k),covarG{k}); % pdf
                p(k) = p(k)/pm(k); % probability of fitting
            end
            
            if (max(p) > tolG) % if it passes any of the gaussians
                mg(j,i) = 1;
            else
                mg(j,i) = 0;
            end

        end
    end
    
    %% --------------------------------------------------------------------
    
    % binary images
    my = imbinarize(my);
    mr = imbinarize(mr);
    mg = imbinarize(mg);
    r = 3;
    SE = strel('sphere',r);
    my = imdilate(my,SE);
    mr = imdilate(mr,SE);
    mg = imdilate(mg,SE);
%     figure(1);imshow(my);title('yellow');
%     figure(2);imshow(mr);title('red');
%     figure(3);imshow(mg);title('green');
    
    my = bwpropfilt(my,'ConvexArea',[450,7854]);
    mr = bwpropfilt(mr,'ConvexArea',[450,7854]);
    mg = bwpropfilt(mg,'ConvexArea',[450,7854]);
%     figure(1);imshow(my);title('yellow');
%     figure(2);imshow(mr);title('red');
%     figure(3);imshow(mg);title('green');
    
    my = bwpropfilt(my,'Extent',1,'largest');
    mr = bwpropfilt(mr,'Extent',1,'largest'); 
    mg = bwpropfilt(mg,'Extent',1,'largest');

%     figure(1);imshow(my);title('yellow');
%     figure(2);imshow(mr);title('red');
%     figure(3);imshow(mg);title('green');
    mb = my + mr + mg;
    
    
    %% --------------------------------------------------------------------
    % centers and contours
    
    r = 1;
    SE = strel('sphere',r);
    
    myc = regionprops(my,'Centroid');
    myp = bwperim(my); myp = imdilate(myp,SE);
%     figure; imshow(myp);title('yellow');
%     myf = imoverlay(my,myp,'y');
%     figure; imshow(myf);title('yellow');
    
    mrc = regionprops(mr,'Centroid');
    mrp = bwperim(mr); mrp = imdilate(mrp,SE);
%     figure;imshow(mrp);title('red');
%     mrf = imoverlay(mr,mrp,'r');
%     figure;imshow(mrf);title('red');
    
    mgc = regionprops(mg,'Centroid');
    mgp = bwperim(mg); mgp = imdilate(mgp,SE);
%     figure;imshow(mgp);title('green');
%     mgf = imoverlay(mg,mgp,'g');
%     figure;imshow(mgf);title('green');
    
    Im2= imoverlay(Im,myp,'y'); % plot on image
    Im2 = imoverlay(Im2,mrp,'r'); % plot on image
    Im2 = imoverlay(Im2,mgp,'g'); % plot on image
    figure(4); imshow(Im2);
    
    %% --------------------------------------------------------------------
    % Identify bouys
    rmin = 12; rmax = 50;
    [centersy,radiy] = imfindcircles(my,[rmin,rmax],'ObjectPolarity','bright');
    [centersr,radir] = imfindcircles(mr,[rmin,rmax],'ObjectPolarity','bright');
    [centersg,radig] = imfindcircles(mg,[rmin,rmax],'ObjectPolarity','bright');
    
    % take only the largest hits
    [~,ixmy] = max(radiy); centersy = centersy(ixmy,:); radiy = radiy(ixmy);
    [~,ixmr] = max(radir); centersr = centersr(ixmr,:); radir = radir(ixmr);
    [~,ixmg] = max(radig); centersg = centersg(ixmg,:); radig = radig(ixmg);

    % final segmentation plot
    t = 1:15:360; % theta
    SE = strel('sphere',1);

    if ~isempty(radiy)
        xy = centersy(1) + radiy*cosd(t); % x-cord
        yy = centersy(2) + radiy*sind(t); % y-cord
        BWy = roipoly(I,xy,yy); % make mask
        BWy = bwperim(BWy); % find outline
        BWy = imdilate(BWy,SE); % make outline thicker
        Im = imoverlay(Im,BWy,'y'); % plot on image
    end
    if ~isempty(radir)
        xr = centersr(1) + radir*cosd(t); % x-cord
        yr = centersr(2) + radir*sind(t); % y-cord
        BWr = roipoly(I,xr,yr); % make mask
        BWr = bwperim(BWr); % find outline
        BWr = imdilate(BWr,SE); % make outline thicker
        Im = imoverlay(Im,BWr,'r'); % plot on image
    end

    if ~isempty(radig)
        xg = centersg(1) + radig*cosd(t); % x-cord
        yg = centersg(2) + radig*sind(t); % y-cord
        BWg = roipoly(I,xg,yg); % make mask
        BWg = bwperim(BWg); % find outline
        BWg = imdilate(BWg,SE); % make outline thicker
        Im = imoverlay(Im,BWg,'g'); % plot on image
    end

    figure(5); imshow(Im);
    % set(gcf,'Position',[1536*0.05 864*0.1 1536*0.6 864*0.8]);

    
    %% --------------------------------------------------------------------
    % save images
    cd(current); cd(plot_path); % change to output folder
    
    imwrite(my,strcat('binary_Y_',num2str(Frame),'.jpg')); % save as binary_<FrameID>.jpg
    imwrite(mr,strcat('binary_R_',num2str(Frame),'.jpg')); % save as binary_<FrameID>.jpg
    imwrite(mg,strcat('binary_G_',num2str(Frame),'.jpg')); % save as binary_<FrameID>.jpg
    imwrite(mb,strcat('binary_',num2str(Frame),'.jpg')); % save as binary_<FrameID>.jpg
    
    imwrite(Im,strcat('out_',num2str(Frame),'.jpg')); % save as out_<FrameNo>.jpg  
    imwrite(Im2,strcat('out2_',num2str(Frame),'.jpg')); % save as out_<FrameNo>.jpg
    
    cd(current); % return
end