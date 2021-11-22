
function [K,Cond,Phase,OMG] = kuramoto_Cond_Hebb(N,Phase0,w0,K0,dt,allstp,fstp,LR,LE,CR,CE,Cond0,towbase,L)


Phase       = zeros(N,allstp); 
Phase(:,1)  = Phase0; 
dPhase      = zeros(N);        
OMG         = zeros(N,allstp);
K           = zeros(N,N,allstp);
tow = floor(( L ./ (Cond0+1)).*towbase);
for t = 1 : fstp-1
    OMG(:,t)        = w0;
    Phase(:,t+1)    = w0 * dt + Phase(:,t);
    K(:,:,t+1) = K0;
    Cond(:,:,t+1) = Cond0;
end
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
    Cond(:,:,t+1) = CR * ( (CE * cos(dPhase) ) - Cond(:,:,t) ) * dt + Cond(:,:,t);
%     Cond(:,:,t+1) = max(0,Cond(:,:,t+1)); % this term is added in order to keep the values positive.
    tow = floor(( L ./ (Cond(:,:,t+1)+1)).*towbase);
end

K       = single(K);
OMG     = single(OMG);
Phase   = single(Phase);

% Phase = zeros(N,allstp); %phases
% Phase(:,1) = Phase0;
% OMG = zeros(N,allstp);
% K = zeros(N,N,allstp);
% tow = floor(( L ./ (Cond0+1)).*towbase);
% for t = 1 : fstp-1
%     for i = 1 : N
%         OMG(i,t) = w0(i);
%         Phase(i,t+1) = w0(i) * dt + Phase(i,t);
%     end
%     K(:,:,t+1) = K0;
%     Cond(:,:,t+1) = Cond0;
%     
% end
% 
% %------ Dynamics dependent & delayed coupled Kuramoto oscillatores --------
% 
% for t = fstp : allstp-1
%     for i = 1:N
%         for j = 1:N
%             Phase(i,t+1) = Phase(i,t+1) + K(i,j,t) * sin(Phase(j,t - tow(i,j))-Phase(i,t));
%         end
%         OMG(i,t) = Phase(i,t+1)/N + w0(i);
%         Phase(i,t+1) = OMG(i,t) * dt + Phase(i,t);    
%     end
%     
%     %             --------------- - adoptive conduction ---------------
% 
%     for i = 1:N
%         for j = 1:N
%              K(i,j,t+1) = LR * ( (LE * cos(Phase(i,t)-Phase(j,t-tow(i,j))) ) - K(i,j,t) ) * dt + K(i,j,t);
%              Cond(i,j,t+1) = CR * ( (CE * cos(Phase(i,t)-Phase(j,t-tow(i,j))) ) - K(i,j,t) ) * dt + K(i,j,t);
%              Cond(i,j,t+1) = max(0,Cond(i,j,t+1)); % this term is added in order to keep the values positive.
%         end
%     end
%     tow = floor(( L ./ (Cond(:,:,t+1)+1)).*towbase);
% end
% 
% K = single(K);
% OMG = single(OMG);
% Phase = single(Phase);