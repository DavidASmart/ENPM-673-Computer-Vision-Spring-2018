function [RGB,BW,points2] = Read_Convert_Detect(LUT,ImageFolder,ImageNames,k)

%% Read Image Files
I = imread(fullfile(ImageFolder,ImageNames(k).name));

%% Convert to Correct Format
% convert from B&W Bayer format to RGB color images using "demosaic".
RGB = demosaic(I,'gbrg');
% Undistort the current frame using "UndistortImage"
RGB = UndistortImage(RGB, LUT);
% smooth image
RGB = imgaussfilt(RGB,1.2);
% a grayscale version is required for feature matching
BW = rgb2gray(RGB);
% adjust contrast for clarity
% BW = histeq(BW);

%% Find Features
points2 = detectSURFFeatures(BW);
% select a subset of points uniformly distributed throughout the image.
% points2 = selectUniform(points2, 256*4, size(BW));


end

