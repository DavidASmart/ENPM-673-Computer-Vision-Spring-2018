function [R, t] = get_Rt(R4, t4)

%% *** Find the "correct" R&t Combo ***

%% Assumption #1: Car is Always Moving Forward
ind = []; % initialize possible R&t Combos
for i = 1:4 % for all R&t combos
    
    if (t4(3,i) > 0) % if delta z is positive
        ind = [ind; i]; % add to list of possible R&t Combos
    end
    
end
R4 = R4(:,:,ind); % reduce to those that satisfy asumption #1
t4 = t4(:,ind);

%% Assumption #2: Car Only Rotates in "Yaw" (aka about the Y axis)
if (size(t4, 2) > 0) % if any R&t Combos satisfy Assumption #1...
    ind = []; % re-initialize possible R&t Combos
    for i = 1:size(t4, 2) % for all remaining R&t Combos
        
        % extract the rotations
        alpha = atan2d(R4(2,1,i),R4(1,1,i)); % yaw
        beta = atan2d(-R4(3,1,i),sqrt(R4(3,2,i)^2+R4(3,3,i)^2)); % pitch
        gamma = atan2d(R4(3,2,i),R4(3,3,i)); %roll

        % check the rotation values
        if abs(alpha) < 90/40 && abs(beta) <= 5 && abs(gamma) < 5
            ind = [ind; i]; % add to list of possible R&t Combos
        end   
    end
    R4 = R4(:,:,ind); % reduce to those that satisfy asumption #2
    t4 = t4(:,ind);

    %% Assumption #3: Car not moving up and down, or sliding side to side
    
    if (size(t4, 2) > 0) % if any R&t Combos satisfy Assumption #1&2...
        
        R = R4(:,:,1); % start with the first combo
        t = t4(:,1);
        min_d = sqrt(t4(1,1)^2 + t4(2,1)^2); % current delta x-y
        
        for i = 2:size(t4,2) % for the rest
            
            d = sqrt(t4(1,i)^2 + t4(2,i)^2); % check the delta x-y
            
            if d < min_d % if less
                min_d = d; % update
                R = R4(:,:,i);
                t = t4(:,i);
            end
        end

    else %% if no R&T Combos satisfy Assumption #1&2...
        % set R & t to "nothing"
        R = eye(3);
        t = [0;0;0];
    end
    
else %% if no R&T Combos satisfy Assumption #1...
    % set R & t to "nothing"
    R = eye(3);
    t = [0;0;0];
end

%% *** REFINE chosen R and t ***

% Assumption #2: Car Only Rotates in "Yaw" (aka about the Y axis)
R(1,2) = 0;
R(2,1) = 0;
R(2,3) = 0;
R(3,2) = 0;
% remove "noise"
if abs(R(1,3)) < 0.005
    R(1,3) = 0;
end
if abs(R(3,1)) < 0.005
    R(3,1) = 0;
end  

% Assumption #1: Car is Always Moving Forward
% Assumption #3: Car not moving up and down, or sliding side to side
if abs(t(1)) < 0 || R(1,1) > 0.90
        t = [0;0;t(3)];
end

end