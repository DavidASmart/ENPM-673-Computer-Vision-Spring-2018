function [] = pos_plot_compare(RGB,pos_my,pos_CVT)

%% PLOT
figure(1); 
set(gcf, 'Position', [1536*0.1 864*0.1 1536*0.8 864*0.8]);

% Video
%subplot(1,2,1);
subplot(2,1,1);
imshow(RGB);

% Visual Odometry
% subplot(1,2,2);
subplot(2,1,2);
plot(-pos_my(:,1),pos_my(:,2),'b','Linewidth',1.5); % mine
hold on;
plot(pos_CVT(:,1),pos_CVT(:,2),'g','Linewidth',1.5); % Matlab's Computer Vision Toolbox

% labels
xlabel('X');
ylabel('Z');
legend('Location','southoutside');
legend('my code','MATLAB CVT');
% view(0,0);

% % update plot volume
% minX = min( min(pos_CVT(1,:)), min(pos_my(1,:)) ) - 5;
% minY = min( min(pos_CVT(2,:)), min(pos_my(2,:)) ) - 5;
% minZ = min( min(pos_CVT(3,:)), min(pos_my(3,:)) ) - 5;
% maxX = max( max(pos_CVT(1,:)), max(pos_my(1,:)) ) + 5;
% maxY = max( max(pos_CVT(2,:)), max(pos_my(2,:)) ) + 5;
% maxZ = max( max(pos_CVT(3,:)), max(pos_my(3,:)) ) + 5;
% axis([minX maxX minY maxY minZ maxZ])

hold off;
pause(0.0001);

end

