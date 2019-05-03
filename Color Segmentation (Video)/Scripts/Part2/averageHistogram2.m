%% averageHistogram2.m
%% PART 2.
% 1.
% Compute the average histogram for each color channel of the training set images. 
% This should help decide on the number of Gaussians to fit the color histogram.

function [Y_data, R_data, G_data] = averageHistogram2()

    % current location is ...
    ScriptsPart2Folder = pwd;
    % want to read cropped images from ...
    CroppedTrainingSetFolder = '../../Images/TrainingSet/CroppedBuoys';
    % want to saves the histogram plots in ...
    HistogramFolder = '../../Output/Part2/';

    cd(CroppedTrainingSetFolder) % change to input folder
    ImageInfo = dir ('*.jpg'); % gather image file names from folder

    % initialize histogram information
    GR = zeros(256,1);GG = zeros(256,1);GB = zeros(256,1);GC = 0;
    RR = zeros(256,1);RG = zeros(256,1);RB = zeros(256,1);RC = 0;
    YR = zeros(256,1);YG = zeros(256,1);YB = zeros(256,1);YC = 0;

    for k = 1:length(ImageInfo) % for all images

        I = imread(ImageInfo(k).name); % read in image files
        % **filtering**
        I = imgaussfilt(I);
        % I = medfilt2(I);
        Iname = ImageInfo(k).name; % image file name

        BW = imbinarize(rgb2gray(I)); % distinguish between colored & none-colored pixels
        PixelList = regionprops(BW,'PixelList'); % locations of colored & none-colored pixels
        if ~isempty(PixelList) % if colored pixels exist

            PixelList = PixelList.PixelList; % extract locations of colored pixels

            R = imhist(I(PixelList(:,2),PixelList(:,1),1)); % red component
            G = imhist(I(PixelList(:,2),PixelList(:,1),2)); % green component
            B = imhist(I(PixelList(:,2),PixelList(:,1),3)); % blue component

            Rd = I(PixelList(:,2),PixelList(:,1),1);
            Gd = I(PixelList(:,2),PixelList(:,1),2); 
            Bd = I(PixelList(:,2),PixelList(:,1),3); 
            
            % sum over all images
            if Iname(1) == 'G' % if green bouy
                GR = GR + R;
                GG = GG + G;
                GB = GB + B;
                GC = GC + 1;
                
                if GC > 1
                    % save all the actual pixel values for making the model
                    Rgdata = [Rgdata , Rd(:)'];
                    Ggdata = [Ggdata , Gd(:)'];
                    Bgdata = [Bgdata , Bd(:)'];
                else
                    % save all the actual pixel values for making the model
                    Rgdata = Rd(:)';
                    Ggdata = Gd(:)';
                    Bgdata = Bd(:)';
                end
                
            elseif Iname(1) == 'R' % if red bouy
                RR = RR + R;
                RG = RG + G;
                RB = RB + B;
                RC = RC + 1;
                
                if RC > 1
                    % save all the actual pixel values for making the model
                    Rrdata = [Rrdata , Rd(:)'];
                    Grdata = [Grdata , Gd(:)'];
                    Brdata = [Brdata , Bd(:)'];
                else
                    % save all the actual pixel values for making the model
                    Rrdata = Rd(:)';
                    Grdata = Gd(:)';
                    Brdata = Bd(:)';
                end
                
            elseif Iname(1) == 'Y' % if yellow bouy
                YR = YR + R;
                YG = YG + G;
                YB = YB + B;
                YC = YC + 1;
                
                if YC > 1
                    % save all the actual pixel values for making the model
                    Rydata = [Rydata , Rd(:)'];
                    Gydata = [Gydata , Gd(:)'];
                    Bydata = [Bydata , Bd(:)'];
                else
                    % save all the actual pixel values for making the model
                    Rydata = Rd(:)';
                    Gydata = Gd(:)';
                    Bydata = Bd(:)';
                end
                
            end

        end
    end

    % average historgrams
    GR = GR/GC; GG = GG/GC; GB = GB/GC;
    RR = RR/RC; RG = RG/RC; RB = RB/RC;
    YR = YR/YC; YG = YG/YC; YB = YB/YC;

    % save as <Color>_hist.jpg
    cd(ScriptsPart2Folder) % first switch to original folder
    cd(HistogramFolder) % then switch to output folder

    figure; plot(GR,'r'); hold on, plot(GG,'g'); plot(GB,'b'); 
    legend(' Red channel','Green channel','Blue channel'); hold off;
    title('GREEN Bouy');
    saveas(gcf,'G_hist.jpg')

    figure; plot(RR,'r'); hold on, plot(RG,'g'); plot(RB,'b'); 
    legend(' Red channel','Green channel','Blue channel'); hold off;
    title('RED Bouy');
    saveas(gcf,'R_hist.jpg')

    figure; plot(YR,'r'); hold on, plot(YG,'g'); plot(YB,'b'); 
    legend(' Red channel','Green channel','Blue channel'); hold off;
    title('YELLOW Bouy');
    saveas(gcf,'Y_hist.jpg')

    % also save numerical values...
    save('G_hist.mat','GR','GG','GB');
    save('R_hist.mat','RR','RG','RB');
    save('Y_hist.mat','YR','YG','YB');
    
    Y_data = [Rydata; Gydata; Bydata];
    R_data = [Rrdata; Grdata; Brdata];
    G_data = [Rgdata; Ggdata; Bgdata];
    
    save('Y_data.mat','Y_data');
    save('R_data.mat','R_data');
    save('G_data.mat','G_data');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd(ScriptsPart2Folder);
end