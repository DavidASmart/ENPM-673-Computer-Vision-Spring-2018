function [R4,t4] = get_4Rt(E)

    %% Single Value Decomposition
    [U,~,V] = svd(E);
    W = [0 -1 0; 1 0 0; 0 0 1];
    
    %% 2 Possible Rotation Matricies
    R1 = U*W*V';
    R2 = U*W'*V';
    
    %% 2 Possible Translation Matricies
    c1 = U(:,3);
    c2 = -U(:,3);

    %% 4 Possible R & T Combos
    % #1 - R1,t1
    R = R1; c = c1;
    if det(R) < 0 % Force det(R) = 1
        R = -R;
        c = -c;
    end
    R4(:,:,1) = R; t4(:,1) = -R'*c;
    
    % #2 - R1,t2
    R = R1; c = c2;
    if det(R) < 0 % Force det(R) = 1
        R = -R;
        c = -c;
    end
    R4(:,:,2) = R; t4(:,2) = -R'*c;
    
    % #3 - R2,t1
    R = R2; c = c1;
    if det(R) < 0 % Force det(R) = 1
        R = -R;
        c = -c;
    end
    R4(:,:,3) = R; t4(:,3) = -R'*c;
    
    % #4 - R2,t2
    R = R2; c = c2;
    if det(R) < 0 % Force det(R) = 1
        R = -R;
        c = -c;
    end
    R4(:,:,4) = R; t4(:,4) = -R'*c;
    
end
