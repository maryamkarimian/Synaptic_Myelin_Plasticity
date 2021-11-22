clear all
close all
clc

nC = 11;
nL = 11;
N = 100;
k = 1;
time = 200;
dt =0.01;
threshold = 0.01;
LE = 1;
p      = 0;
sigma   = 0.01; % The standard deviation of the intrinsic frequency normal distanceribution in Hz.
mu      = 1; % The mean value of the intrinsic frequensies in Hz.
%v       = 0.1 ; % The velosity of signal propagation in terms of cm per time step(10^(-5))s.
fstp     = 600; % steps are considered to be passed without time delay and Dynamics coupling,...
%               so that I have memory in my data to be used for the next steps that I have time delay.
%               The number of these steps should be the max amount of time delay(max distanceance/v)...
%               in terms of time step.
LRall = linspace(0,0.1,nL);
LEall = linspace(0,1,nL);
CRall = linspace(0,0.1,nC);
CEall = linspace(0,1,nC);
allstp = floor(time/dt);
load InitialCondition.mat;
Phase0 = InitialCondition.Phase0;
K0   = InitialCondition.K0;
Cond0   = InitialCondition.Cond0;
Cond0     = Cond0*0.04; % such that T0 = 7
towbase  = InitialCondition.towbase;
w0   = InitialCondition.w0;

tic;
dir1 = 'P:\FSE_MACSBIO\maryam.karimian\data\Timm_CondPlasticity_lowCR_InitT7_1st/';
L = 1; %so, the max(tow) would be 5000 and the min(tow) would be 10. (T(cond<0.1)=L/0.1)
% pool = parpool;
LE = 1;
%for iLR = 1:2
    %LR = (iLR - 1)*0.1;
    LR = 0.1;
    for iCR = 1:numel(CRall)
        CR = CRall(iCR);
        for iCE = 1:numel(CEall)
%         parfor iCE = 1:numel(CEall)
            CE = CEall(iCE);
            [K,Cond,Phase,OMG] = kuramoto_Cond_Hebb(N,Phase0,w0,K0,dt,allstp,fstp,LR,LE,CR,CE,Cond0,towbase,L);
            fname = [dir1,'LR',num2str(LR),'_LE',num2str(LE),'_CR',num2str(CR),'_CE',num2str(CE),'.mat'];
            parsave_cond_Hebb(fname,K,Cond,OMG,Phase,allstp,threshold,N);
        end
    end
%end

toc;
% delete(pool);
disp ('Kuramoto finished')
%%
pool = parpool;
tic;
% R1 = zeros(numel(CRall),numel(CEall),numel(LRall),numel(LEall));
% R2 = zeros(numel(CRall),numel(CEall),numel(LRall),numel(LEall));
R1 = zeros(numel(CRall),numel(CEall),numel(LRall));
R2 = zeros(numel(CRall),numel(CEall),numel(LRall));
for iLR = 1:numel(LRall)
    LR = LRall(iLR);
%     for iLE = 1:numel(LEall)
%         LE = LEall(iLE);
        for iCR = 1:numel(CRall)
            CR = CRall(iCR);
            parfor iCE = 1:numel(CEall)
                CE = CEall(iCE);
%                 fname = [dir1,'LR',num2str(LR),'_LE',num2str(LE),'_CR',num2str(CR),'_CE',num2str(CE),'.mat'];
                fname = [dir1,'LR',num2str(LR),'_CR',num2str(CR),'_CE',num2str(CE),'.mat'];
                [Phase,K2] = load_func_cond_Hebb(fname);
                [R1_t,R2_t,tt] = OrderParameter_mod(N,dt,allstp,Phase,K2);
%                 R1(iCR,iCE,iLR,iLE) = mean(R1_t(1,(allstp/100 + 1)-10:(allstp/100 + 1)));
%                 R2(iCR,iCE,iLR,iLE) = mean(R2_t(1,(allstp/100 + 1)-10:(allstp/100 + 1)));
                R1(iCR,iCE,iLR) = mean(R1_t(1,(allstp/100 + 1)-10:(allstp/100 + 1)));
                R2(iCR,iCE,iLR) = mean(R2_t(1,(allstp/100 + 1)-10:(allstp/100 + 1)));
            end
        end
    %end
end
% OrderParameter_CRCELRLE_td.R1 = R1;
% OrderParameter_CRCELRLE_td.R2 = R2;
% save ( 'OrderParameter_CRCELRLE_td.mat', 'OrderParameter_CRCELRLE_td', '-v7.3' )
OrderParameter_CRCELR_td.R1 = R1;
OrderParameter_CRCELR_td.R2 = R2;
save ( 'OrderParameter_CRCELR_td.mat', 'OrderParameter_CRCELR_td', '-v7.3' )
delete(pool);

toc;
