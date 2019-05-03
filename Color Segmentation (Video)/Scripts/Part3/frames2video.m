%% PART 3.
% Buoy  Detection (40) 
%
% 3.
%
% Write a function to generate a videos equence by reading images from 
% ColorSeg/Output/Part3/Frames/
% and save it in 
% ColorSeg/Output/Part3/Video/
 
function [] = frames2video(inpath,outpath)

    current = pwd;
    
    cd(inpath);
    % read all the frames in order
    for Frame = 1:200
        I = imread(strcat('out_',num2str(Frame),'.jpg'));
        imshow(I);
        set(gcf,'Position',[1536*0.05 864*0.1 1536*0.8 864*0.8]);
        F(Frame) = getframe(gcf); % save figure
    end

    cd(current); cd(outpath);
    % Write to File
    vW = VideoWriter('Final.avi');
    vW.FrameRate = 5;
    open(vW)
    for Frame = 1:200
        writeVideo(vW,F(Frame).cdata)
    end
    close(vW)
    
    %% --------------------------------------------------------------------
    
    cd(current); cd(inpath);
    % read all the frames in order
    for Frame = 1:200
        I = imread(strcat('out2_',num2str(Frame),'.jpg'));
        imshow(I);
        set(gcf,'Position',[1536*0.05 864*0.1 1536*0.8 864*0.8]);
        F(Frame) = getframe(gcf); % save figure
    end

    cd(current); cd(outpath);
    % Write to File
    vW = VideoWriter('Final2.avi');
    vW.FrameRate = 5;
    open(vW)
    for Frame = 1:200
        writeVideo(vW,F(Frame).cdata)
    end
    close(vW)


end