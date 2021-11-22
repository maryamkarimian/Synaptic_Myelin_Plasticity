function [K,Phase,OMG] = kuramoto_Timm(N,Phase0,w0,K0,dt,allstp,fstp,LR,LE,tow)

Phase       = zeros(N,allstp); % phases
Phase(:,1)  = Phase0; 
dPhase      = zeros(N);        % all pairs of phase offsets
OMG         = zeros(N,allstp);
K           = ones(N,N,allstp);

for t = 1 : fstp-1
    OMG(:,t)        = w0;
    Phase(:,t+1)    = w0 * dt + Phase(:,t);
    K(:,:,t+1) = K0;
end

%------ Dynamics dependent & delayed coupled Kuramoto oscillatores --------

for t = fstp : allstp-1
    for i = 1:N
        for j = 1:N
            dPhase(i,j) = Phase(j,t - tow(i,j))-Phase(i,t);
        end   
    end
    OMG(:,t)     = sum(K(:,:,t) .* sin(dPhase),2) / N + w0;
    Phase(:,t+1) = OMG(:,t) * dt + Phase(:,t);
    
    % ------------------ adaptive connection weights ----------------------
    K(:,:,t+1)   = LR * ( (LE * cos(dPhase) ) - K(:,:,t) ) * dt + K(:,:,t);
end

K       = single(K);
OMG     = single(OMG);
Phase   = single(Phase);