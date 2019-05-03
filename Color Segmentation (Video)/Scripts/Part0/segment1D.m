function [] = segment1D(Frame)
%% PART 0.
% 3. 
% Model the color distributions of the buoys using 1-D Gaussians. 
% Segment the buoys based on their color (red green and yellow).

% current folder is ...
ScriptsPart0Folder = pwd;
% need to read in training data from both
TrainingSetFolder = '../../Images/TrainingSet/Frames';
TestSetFolder = '../../Images/TestSet/Frames';
% need to output to ...
OutputFolder = '../../Output/Part0/';
% need to read in the average histogram data from that same folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODELS

% reterive the average histograms for each of the bouys
cd(OutputFolder);
load('Y_hist.mat'); % yellow bouy
YR = YR*100./sum(YR);YG = YG*100./sum(YG);YB = YB*100./sum(YB);
load('R_hist.mat'); % red bouy
RR = RR*100./sum(RR);RG = RG*100./sum(RG);RB = RB*100./sum(RB);
load('G_hist.mat'); % green bouy
GR = GR*100./sum(GR);GG = GG*100./sum(GG);GB = GB*100./sum(GB);


%% yellow bouy
% convert histogram back into 1D-data
tr = 1; tg = 1; tb = 1;
for i = 1:256
    dtr = round(YR(i));
    if dtr > 0
        yr(tr:tr+dtr) = i;
        tr = tr + dtr + 1;
    end
    dtg = round(YG(i));
    if dtg > 0
        yg(tg:tg+dtg) = i;
        tg = tg + dtg + 1;
    end
    dtb = round(YB(i));
    if dtb > 0
        yb(tb:tb+dtb) = i;
        tb = tb + dtb + 1;
    end
end

% initial cut
muyr = mean(yr);
sigmayr = std(yr);
muyg = mean(yg);
sigmayg = std(yg);
muyb = mean(yb);
sigmayb = std(yb);

% only take the larger cluster
for i = size(yr,2):-1:1
    if yr(i) < muyr - 2*sigmayr
        yr(i) = [];
    end
end
for i = size(yg,2):-1:1
    if yg(i) < muyg - 2*sigmayg
        yg(i) = [];
    end
end
for i = size(yb,2):-1:1
    if yb(i) < muyb - 2*sigmayb
        yb(i) = [];
    end
end

% second cut
muyr = mean(yr);
sigmayr = std(yr);
muyg = mean(yg);
sigmayg = std(yg);
muyb = mean(yb);
sigmayb = std(yb);

%% red bouy
% convert histogram back into 1D-data
tr = 1; tg = 1; tb = 1;
for i = 1:256
    dtr = round(RR(i));
    if dtr > 0
        rr(tr:tr+dtr) = i;
        tr = tr + dtr + 1;
    end
    dtg = round(RG(i));
    if dtg > 0
        rg(tg:tg+dtg) = i;
        tg = tg + dtg + 1;
    end
    dtb = round(RB(i));
    if dtb > 0
        rb(tb:tb+dtb) = i;
        tb = tb + dtb + 1;
    end
end

% initial cut
murr = mean(rr);
sigmarr = std(rr);
murg = mean(rg);
sigmarg = std(rg);
murb = mean(rb);
sigmarb = std(rb);

% only take the larger cluster
for i = size(rr,2):-1:1
    if rr(i) < murr - 2*sigmarr
        rr(i) = [];
    end
end
for i = size(rg,2):-1:1
    if rg(i) < murg - 2*sigmarg
        rg(i) = [];
    end
end
for i = size(rb,2):-1:1
    if rb(i) < murb - 2*sigmarb
        rb(i) = [];
    end
end

% second cut
murr = mean(rr);
sigmarr = std(rr);
murg = mean(rg);
sigmarg = std(rg);
murb = mean(rb);
sigmarb = std(rb);

%% green bouy
% convert histogram back into 1D-data
tr = 1; tg = 1; tb = 1;
for i = 1:256
    dtr = round(GR(i));
    if dtr > 0
        gr(tr:tr+dtr) = i;
        tr = tr + dtr + 1;
    end
    dtg = round(GG(i));
    if dtg > 0
        gg(tg:tg+dtg) = i;
        tg = tg + dtg + 1;
    end
    dtb = round(GB(i));
    if dtb > 0
        gb(tb:tb+dtb) = i;
        tb = tb + dtb + 1;
    end
end

% initial cut
mugr = mean(gr);
sigmagr = std(gr);
mugg = mean(gg);
sigmagg = std(gg);
mugb = mean(gb);
sigmagb = std(gb);

% only take the larger cluster
for i = size(gr,2):-1:1
    if gr(i) < mugr - 2*sigmagr
        gr(i) = [];
    end
end
for i = size(gg,2):-1:1
    if gg(i) < mugg - 2*sigmagg
        gg(i) = [];
    end
end
for i = size(gb,2):-1:1
    if gb(i) < mugb - 2*sigmagb
        gb(i) = [];
    end
end

% second cut
mugr = mean(gr);
sigmagr = std(gr);
mugg = mean(gg);
sigmagg = std(gg);
mugb = mean(gb);
sigmagb = std(gb);


%% plot models (only needs to be done once)--------------------------------
if Frame == 1
    
    cd(ScriptsPart0Folder);
    cd(OutputFolder);
    
    %% yellow bouy
    % plot average histograms
    figure(1); plot(YR,'r'); hold on; plot(YG,'g'); plot(YB,'b');
    % plot synthesized data
    plot(yr,zeros(length(yr),1),'r.');
    plot(yg,zeros(length(yg),1),'g.'); 
    plot(yb,zeros(length(yb),1),'b.');
    % plot the gaussian models
    x = 1:256;
    pdfyr = (1/(sqrt(2*pi)*sigmayr)) * exp((-1/(2*sigmayr^2))*(x - muyr).^2);
    plot(x, pdfyr*100,'r--');
    pdfyg = (1/(sqrt(2*pi)*sigmayg)) * exp((-1/(2*sigmayg^2))*(x - muyg).^2);
    plot(x, pdfyg*100,'g--');
    pdfyb = (1/(sqrt(2*pi)*sigmayb)) * exp((-1/(2*sigmayb^2))*(x - muyb).^2);
    plot(x, pdfyb*100,'b--'); axis([1 256 0 10]); title('yellow bouy'); hold off;
    % Save the plot of the 1-D Gaussian Model  
    saveas(gcf,'gauss1D_y.jpg');
    
    %% red bouy
    % plot average histograms
    figure(2); plot(RR,'r'); hold on; plot(RG,'g'); plot(RB,'b'); hold on;
    % plot synthesized data
    plot(rr,zeros(length(rr),1),'r.'); 
    plot(rg,zeros(length(rg),1),'g.'); 
    plot(rb,zeros(length(rb),1),'b.');
    % plot the gaussian models
    pdfrr = (1/(sqrt(2*pi)*sigmarr)) * exp((-1/(2*sigmarr^2))*(x - murr).^2);
    plot(x, pdfrr*100,'r--');
    pdfrg = (1/(sqrt(2*pi)*sigmarg)) * exp((-1/(2*sigmarg^2))*(x - murg).^2);
    plot(x, pdfrg*100,'g--');
    pdfrb = (1/(sqrt(2*pi)*sigmarb)) * exp((-1/(2*sigmarb^2))*(x - murb).^2);
    plot(x, pdfrb*100,'b--'); axis([1 256 0 10]); title('red bouy'); hold off;
    % Save the plot of the 1-D Gaussian Model
    saveas(gcf,'gauss1D_r.jpg');
    
    %% green bouy
    % plot average histograms
    figure(3); plot(GR,'r'); hold on; plot(GG,'g'); plot(GB,'b'); hold on;
    % plot synthesized data
    plot(gr,zeros(length(gr),1),'r.'); 
    plot(gg,zeros(length(gg),1),'g.'); 
    plot(gb,zeros(length(gb),1),'b.');
    % plot the gaussian models
    pdfgr = (1/(sqrt(2*pi)*sigmagr)) * exp((-1/(2*sigmagr^2))*(x - mugr).^2);
    plot(x, pdfgr*100,'r--');
    pdfgg = (1/(sqrt(2*pi)*sigmagg)) * exp((-1/(2*sigmagg^2))*(x - mugg).^2);
    plot(x, pdfgg*100,'g--');
    pdfgb = (1/(sqrt(2*pi)*sigmagb)) * exp((-1/(2*sigmagb^2))*(x - mugb).^2);
    plot(x, pdfgb*100,'b--'); axis([1 256 0 10]); title('green bouy'); hold off;
    % Save the plot of the 1-D Gaussian Model 
    saveas(gcf,'gauss1D_g.jpg');
    
    clc
    close all
    
end


%% segmentation! ----------------------------------------------------------

% try both image folders to find the correct image
cd(ScriptsPart0Folder); cd(TrainingSetFolder); % change to input folder
ImageInfo = dir ('*.jpg'); % gather image file names from folder
for i = 1:length(ImageInfo)
    c(i) = strcmp(ImageInfo(i).name,strcat(num2str(Frame),'.jpg'));
end
[found,id] = max(c); % in trainingset?

if ~found % in testset
    cd(ScriptsPart0Folder); cd(TestSetFolder); % change to input folder
    ImageInfo = dir ('*.jpg'); % gather image file names from folder
    for i = 1:length(ImageInfo)
        c(i) = strcmp(ImageInfo(i).name,strcat(num2str(Frame),'.jpg'));
    end
    [~,id] = max(c); % in trainingset?
end

I = imread(ImageInfo(id).name); % read in image file
Iname = ImageInfo(id).name; % image file name

% **filtering**
I2 = imgaussfilt(I);
% I2 = medfilt2(I);

% extract colors
R = im2double(I2(:,:,1)); % figure(5); imshow(R); title('R Image');
G = im2double(I2(:,:,2)); % figure(6); imshow(G); title('G Image');
B = im2double(I2(:,:,3)); % figure(7); imshow(B); title('B Image');

% initialize masks
my = zeros(size(I,1),size(I,2)); % yellow bouy mask
mr = zeros(size(I,1),size(I,2)); % red bouy mask
mg = zeros(size(I,1),size(I,2)); % green bouy mask

% max values of pdf functions
pyrm = (1/(sqrt(2*pi)*sigmayr));
pygm = (1/(sqrt(2*pi)*sigmayg));
pybm = (1/(sqrt(2*pi)*sigmayb));

prrm = (1/(sqrt(2*pi)*sigmarr));
prgm = (1/(sqrt(2*pi)*sigmarg));
prbm = (1/(sqrt(2*pi)*sigmarb));

pgrm = (1/(sqrt(2*pi)*sigmagr));
pggm = (1/(sqrt(2*pi)*sigmagg));
pgbm = (1/(sqrt(2*pi)*sigmagb));

tolY = 0.6;
tolR = 0.4;
tolG = 0.75;

% threshold all pixels based on 1D-gaussian models
for i = 1:size(G,2) % all x
    for j = 1:size(G,1) % all y
        
        % probability of fitting the gaussian model for the yellow bouy?
        pyr = (1/(sqrt(2*pi)*sigmayr)) * exp((-1/(2*sigmayr^2))*(R(j,i)*255 - muyr).^2);
        pyg = (1/(sqrt(2*pi)*sigmayg)) * exp((-1/(2*sigmayg^2))*(G(j,i)*255 - muyg).^2);
        pyb = (1/(sqrt(2*pi)*sigmayb)) * exp((-1/(2*sigmayb^2))*(B(j,i)*255 - muyb).^2);
        
        if (pyr/pyrm > tolY) && (pyg/pygm > tolY) && (pyb/pybm > tolY)
            my(j,i) = 1;
        else
            my(j,i) = 0;
        end
        
        % red bouy?
        prr = (1/(sqrt(2*pi)*sigmarr)) * exp((-1/(2*sigmarr^2))*(R(j,i)*255 - murr).^2);
        prg = (1/(sqrt(2*pi)*sigmarg)) * exp((-1/(2*sigmarg^2))*(G(j,i)*255 - murg).^2);
        prb = (1/(sqrt(2*pi)*sigmarb)) * exp((-1/(2*sigmarb^2))*(B(j,i)*255 - murb).^2);
        
        if (prr/prrm > tolR) && (prg/prgm > tolR) && (prb/prbm > tolR)
            mr(j,i) = 1;
        else
            mr(j,i) = 0;
        end
        
        % green bouy?
        pgr = (1/(sqrt(2*pi)*sigmagr)) * exp((-1/(2*sigmagr^2))*(R(j,i)*255 - mugr).^2);
        pgg = (1/(sqrt(2*pi)*sigmagg)) * exp((-1/(2*sigmagg^2))*(G(j,i)*255 - mugg).^2);
        pgb = (1/(sqrt(2*pi)*sigmagb)) * exp((-1/(2*sigmagb^2))*(B(j,i)*255 - mugb).^2);
        
        if (pgr/pgrm > tolG) && (pgg/pggm > tolG) && (pgb/pgbm > tolG)
            mg(j,i) = 1;
        else
            mg(j,i) = 0;
        end
        
    end
end

my = imbinarize(my);
mr = imbinarize(mr);
mg = imbinarize(mg);
r = 3;
SE = strel('sphere',r);
my = imdilate(my,SE);
mr = imdilate(mr,SE);
mg = imdilate(mg,SE);
% figure;imshow(my);title('yellow');
% figure;imshow(mr);title('red');
% figure;imshow(mg);title('green');

% Remove "Artifacts"
my = bwpropfilt(my,'ConvexArea',[400,4000]);
mr = bwpropfilt(mr,'ConvexArea',[400,4000]);
mg = bwpropfilt(mg,'ConvexArea',[400,4000]);
% figure;imshow(my);title('yellow');
% figure;imshow(mr);title('red');
% figure;imshow(mg);title('green');

% my = bwpropfilt(my,'Extent',[0.5,1]);
% mr = bwpropfilt(mr,'Extent',[0.5,1]);
% mg = bwpropfilt(mg,'Extent',[0.5,1]);
my = bwpropfilt(my,'Extent',1,'largest');
mr = bwpropfilt(mr,'Extent',1,'largest'); 
mg = bwpropfilt(mg,'Extent',1,'largest');
% figure;imshow(my);title('yellow');
% figure;imshow(mr);title('red');
% figure;imshow(mg);title('green');

% Identify bouys
rmin = 10; rmax = 60;
[centersy,radiy] = imfindcircles(my,[rmin,rmax],'ObjectPolarity','bright');
[centersr,radir] = imfindcircles(mr,[rmin,rmax],'ObjectPolarity','bright');
[centersg,radig] = imfindcircles(mg,[rmin,rmax],'ObjectPolarity','bright');

% take only the largest hits
[~,ixmy] = max(radiy); centersy = centersy(ixmy,:); radiy = radiy(ixmy);
[~,ixmr] = max(radir); centersr = centersr(ixmr,:); radir = radir(ixmr);
[~,ixmg] = max(radig); centersg = centersg(ixmg,:); radig = radig(ixmg);

% final segmentation plot
h = figure(11); hold on
t = 1:15:360; % theta
SE = strel('sphere',1);

if ~isempty(radiy)
    xy = centersy(1) + radiy*cosd(t); % x-cord
    yy = centersy(2) + radiy*sind(t); % y-cord
    BWy = roipoly(I,xy,yy); % make mask
    BWy = bwperim(BWy); % find outline
    BWy = imdilate(BWy,SE); % make outline thicker
    I = imoverlay(I,BWy,'y'); % plot on image
end
if ~isempty(radir)
    xr = centersr(1) + radir*cosd(t); % x-cord
    yr = centersr(2) + radir*sind(t); % y-cord
    BWr = roipoly(I,xr,yr); % make mask
    BWr = bwperim(BWr); % find outline
    BWr = imdilate(BWr,SE); % make outline thicker
    I = imoverlay(I,BWr,'r'); % plot on image
end

if ~isempty(radig)
    xg = centersg(1) + radig*cosd(t); % x-cord
    yg = centersg(2) + radig*sind(t); % y-cord
    BWg = roipoly(I,xg,yg); % make mask
    BWg = bwperim(BWg); % find outline
    BWg = imdilate(BWg,SE); % make outline thicker
    I = imoverlay(I,BWg,'g'); % plot on image
end

imshow(I);
set(gcf,'Position',[1536*0.05 864*0.1 1536*0.6 864*0.8]);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the segmented images
cd(ScriptsPart0Folder);
cd(OutputFolder); % change to output folder
imwrite(I,strcat('seg_',Iname)); % save as seg_<FrameNo>.jpg

end

