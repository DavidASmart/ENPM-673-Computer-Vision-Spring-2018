% PROJECT1_Lane_Detection_regular.m 
% Simple Hough Lines Approach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% David Smart
% ENPM 673 - Perception
% University of Maryland, College Park
% Project 1: Lane Detection
% 2/25/2018 (due: 2/28/2018)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem Statement
% In this project we aim to do simple Lane Detection 
% to mimic Lane Departure Warning systems used in Self Driving Cars
%
% Submit a report explaining your approach with relevant outputs after each step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear;
clc

% Video Frames
vO = VideoReader('project_video.mp4');
k = 1;
while vO.CurrentTime <= 15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original frame
frame = im2double(readFrame(vO));

% Crop to Region Of Interest 
frame2 = imcrop(frame,[350 500 740 200]);

% comnvert to B&W
frame3 = rgb2gray(frame2);

% % Denoise the image
m = 5;n = 5;frame3 = medfilt2(frame3,[m n]); % median

% mask out unnessessary stuff
x = [1,740,740,555,370,555,185,370,185,1,1];
y = [1,1,  100,1,  1,  200,200,1,  1, 100,1];
mask = roipoly(frame3,x,y);
frame4 = regionfill(frame3,mask);

% Binarize
T = 0.7;
frame4 = imbinarize(frame4,T);
frame4 = bwpropfilt(frame4,'MajorAxisLength',[30 1000]);
frame4 = bwpropfilt(frame4,'MinorAxisLength',[1 15]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edge detection
frame5 = edge(frame4,'sobel');

% Get Hough Lines
[H,theta,rho] = hough(frame5);

% Find longest hough lines
P = houghpeaks(H,4);
lines = houghlines(frame5,theta,rho,P,'FillGap',50,'MinLength',10);

% group longest hough lines by positive/negative theta and extrapolate
Ry = []; Ly = [];
for i = 1:length(lines)
   if lines(i).theta < -45 && lines(i).theta > -75 && isempty(Ry)
       Rm = (lines(i).point2(2)-lines(i).point1(2))/(lines(i).point2(1)-lines(i).point1(1));
       Ri = lines(i).point1(2)-Rm*(lines(i).point1(1)+350)+500;
       Rx = 1:vO.Width;
       Ry = Rm*Rx + Ri;
       if ~isempty(Ly)
           break
       end
   elseif lines(i).theta > 45 && lines(i).theta < 75 && isempty(Ly)
       Lm = (lines(i).point2(2)-lines(i).point1(2))/(lines(i).point2(1)-lines(i).point1(1));
       Li = lines(i).point1(2)-Lm*(lines(i).point1(1)+350)+500;
       Lx = 1:vO.Width;
       Ly = Lm*Lx + Li;
       if ~isempty(Ry)
           break
       end
   end
end
% if a line was not found, then assume nothing has changed and use the old one
R_old = 0; L_old = 0;
if isempty(Ry)
   Ry = Ry_old;
   R_old = 1;
end
if isempty(Ly)
   Ly = Ly_old;
   L_old = 1;
end

% turn direction?
Right = 0;Left = 0;
dif = 100*abs((Ly - Ry)./Ry);
vanishing_point = Lx(min(find(dif < 2)));
if vanishing_point > vO.Width/2 + vO.Width/10
    Right = 1;
elseif vanishing_point < vO.Width/2 + vO.Width/20
    Left = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY
figure(1);
set(gcf, 'Position', [1536*0.1 864*0.1 1536*0.8 864*0.8]);
ax1 = subplot(2,2,1);imshow(frame3);title('Denoised');
ax2 = subplot(2,2,2);imshow(frame4);title('Hugh Lines');
hold on
for i = 1:length(lines)
   xy = [lines(i).point1; lines(i).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','r');
end
plot(x,y,':y','LineWidth',1.5);
hold off

ax3 = subplot(2,2,[3,4]);imshow(frame);title('Lane');
hold on
if L_old == 1
    plot(Lx(1:vO.Width/2),Ly(1:vO.Width/2),'--','LineWidth',2,'Color','g');
else
    plot(Lx(1:vO.Width/2),Ly(1:vO.Width/2),'LineWidth',2,'Color','g');
end
if R_old == 1
    plot(Rx(vO.Width/2:vO.Width),Ry(vO.Width/2:vO.Width),'--','LineWidth',2,'Color','g');
else
    plot(Rx(vO.Width/2:vO.Width),Ry(vO.Width/2:vO.Width),'LineWidth',2,'Color','g');
end
if Right
    text(24,24,'Right','Color','red','FontSize',14)
elseif Left
    text(24,24,'Left','Color','red','FontSize',14)
else
    text(24,24,'Straight','Color','r','FontSize',14)
end
text(vO.Width*4/5,24,strcat(num2str(round(vO.CurrentTime,3)),'s'),'Color','r','FontSize',14)
hold off

% save frame for real-time play back
F(k) = getframe(ax3);

% save the lines in case they are temporarily lost
Ry_old = Ry; Ly_old = Ly;

k = k+1; % move on
end

%% play back in real-time
close all
imshow(F(1).cdata, 'Border', 'tight')
movie(F,3,vO.FrameRate)

%% Write to File
vW = VideoWriter('P_regular.avi');
vW.FrameRate = 25;
open(vW)
for i = 1:k-1
    writeVideo(vW,F(i).cdata)
end
close(vW)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

