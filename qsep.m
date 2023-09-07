function [Qb] = qsep(Q,beta)
%% Lyne-Hollick recursive filtering method for separating direct runoff and baseflow

    % Direct runoff calculations
    Qd = zeros(size(Q));
    % forward
    for i = 2:length(Q)
        Qd(i,1) = beta*Qd(i-1,1) + (1+beta)/2*(Q(i,1)-Q(i-1,1)); 
        if Qd(i,1) < 0
            Qd(i,1) = 0;
        end
    end
    % backward
    for i = length(Q)-1:-1:1
        Qd(i,1) = beta*Qd(i+1,1) + (1+beta)/2*(Qd(i,1)-Qd(i+1,1)); 
        if Qd(i,1) < 0
            Qd(i,1) = 0;
        elseif Qd(i,1) > Q(i,1)
            Qd(i,1) = Q(i,1);    
        end
    end
    % forward
    for i = 2:length(Q)
        Qd(i,1) = beta*Qd(i-1,1) + (1+beta)/2*(Qd(i,1)-Qd(i-1,1)); 
        if Qd(i,1) < 0
            Qd(i,1) = 0;
        end
    end
    Qb = Q-Qd;
end