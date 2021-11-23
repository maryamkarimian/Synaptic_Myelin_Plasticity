
%% Connection matrix
N = 100;
time = 200;
dt = 0.01;
LR = 0.1;
LE = 1;
allstp = floor(time/dt);
nL = 11;
nC = 11;
%LRall = linspace(0,0.1,nL);
CRall = linspace(0,0.1,nC);
CEall = linspace(0,1,nC);

tLim = 100;
dir1 = 'P:\FSE_MACSBIO\maryam.karimian\data\Timm_CondPlasticity_lowCR_InitT7_1st/';
dir2 = 'P:\FSE_MACSBIO\maryam.karimian\data\output\Timm_CondPlasticity_lowCR_InitT7_1st\Connections/';
fig1 = figure;
for iCR = 1:numel(CRall)
    CR = CRall(iCR);
    for iCE = 1:numel(CEall)
        CE = CEall(iCE);
        load([dir1,'LR',num2str(LR),'_LE',num2str(LE),'_CR',num2str(CR),'_CE',num2str(CE),'.mat']);

        loc = 11*11-( ( (iCR-1)*11 )+(11-iCE) );
        subplot(11,11,loc)
        imagesc(K2(:,:,end));
        colormap Jet;
%         colorbar;
        caxis([-1 1]);
        x = gca;
        xlim([1 N])
        ylim([1 N])
        
        x.XDir = 'normal';
        x.YDir = 'normal';
        x.XTick = [1:49:N];
        x.YTick = [1:49:N];
        x.XTickLabel = {'1', '50', '100'};
        x.YTickLabel = {'1', '50', '100'};
        pbaspect([1 1 1])

    end
end
tightfig;
fig1 = gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 20 10];
fnametif = [dir2,'Connection_LR',num2str(LR),'.tif'];
saveas(fig1, fnametif);
fnamefig = [dir2,'Connection_LR',num2str(LR),'.fig'];
saveas(fig1, fnamefig);
hold off;
%% Single connection matrix

LR = 0;
LE = 1;
CR = 0.06;
CE = 0;
nC = 11;
nL = 11;
time = 200;
dt = 0.01;
tLim = 100;
allstp = floor(time/dt);
LRall = linspace(0,0.1,nL);
dir1 = 'P:\FSE_MACSBIO\maryam.karimian\data\Timm_CondPlasticity_lowCR_InitT7_1st/';
dir2 = 'P:\FSE_MACSBIO\maryam.karimian\data\output\Timm_CondPlasticity_lowCR_InitT7_1st\Connections/';

load([dir1,'LR',num2str(LR),'_LE',num2str(LE),'_CR',num2str(CR),'_CE',num2str(CE),'.mat']);
fig1 = figure;
imagesc(K2(:,:,end));
colormap Jet;
caxis([-1 1]);
x = gca;
xlim([1 N])
ylim([1 N])

x.XDir = 'normal';
x.YDir = 'normal';

x.XTick = [1:19:N];
x.YTick = [1:19:N];
pbaspect([1 1 1])
x.XTickLabel = {'1' , '20', '40' , '60' , '80', '100'};
x.YTickLabel = {'1' , '20', '40' , '60' , '80', '100'};
set(gca,'FontSize',16);
titname = ['\eta=',num2str(CR),', \beta=',num2str(CE)];
title(titname);
fnametif = [dir2,'LR',num2str(LR),'CR',num2str(CR),'_CE',num2str(CE),'.tif'];
saveas(fig1, fnametif);
fnamefig = [dir2,'LR',num2str(LR),'CR',num2str(CR),'_CE',num2str(CE),'.fig'];
saveas(fig1, fnamefig);

%% Correlation Matrix


LR = 0.1;
LE = 1;
tLim = 100;
dir1 = 'P:\FSE_MACSBIO\maryam.karimian\data\Timm_CondPlasticity_lowCR_InitT7_22d/';
dir2 = 'P:\FSE_MACSBIO\maryam.karimian\data\output\Timm_CondPlasticity_lowCR_InitT7_22d/Correlations\';
fig2 = figure;
for iCR = 1:numel(CRall)
    CR = CRall(iCR);
    for iCE = 1:numel(CEall)
        CE = CEall(iCE);
        load([dir1,'LR',num2str(LR),'_LE',num2str(LE),'_CR',num2str(CR),'_CE',num2str(CE),'.mat']);
        %[D] = correlation(tLim,Phase,N);
        %-------------------------------
        D = zeros(N,N);
        l = size(Phase,2)-tLim;
        %l = 9900-tLim;
        for k = 1:tLim
            for i = 1:N
                for j = i:N
                    D(i,j) = D(i,j) + cos(Phase(i,l)-Phase(j,l)) / tLim;
                    D(j,i) = D(i,j);
                end
            end
            l = l + 1;
        end
        %-------------------------------
%         fig2 = figure;
%         %-------Reordering
%         W = K2(:,:,1);
%         [M,Q] = community_louvain(W,gamma);
%         M0 = M;
%         Ph = zeros(N,1);
%         W = K2(:,:,end);
%         [M,Q] = community_louvain(W,gamma,M0,'negative_asym');
%         [On,Wr] = reorder_mod(W,M);
%         for i = 1:N
%             Ph(i) = Phase(On(i),end);
%         end
%         D = zeros(N,N);
%         for i = 1:N
%             for j = i:N
%                 D(i,j) = cos(Ph(i,end)-Ph(j,end));
%                 D(j,i) = D(i,j);
%             end
%         end
        %---------------------------------
%         loc = (iCR-1)*11 + iCE;
        loc = 11*11-( ( (iCR-1)*11 )+(11-iCE) );
        subplot(11,11,loc)
        imagesc(D);
        colormap Jet;
%         colorbar;
        caxis([-1 1]);
        x = gca;
        xlim([1 N])
        ylim([1 N])
        
        x.XDir = 'normal';
        x.YDir = 'normal';
        x.XTick = [1:49:N];
        x.YTick = [1:49:N];
        x.XTickLabel = {'1', '50', '100'};
        x.YTickLabel = {'1', '50', '100'};
        pbaspect([1 1 1])
        
    end
end
tightfig;
fig1 = gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 20 10];
fnametif = [dir2,'Correlation_LR',num2str(LR),'.tif'];
saveas(fig1, fnametif);
fnamefig = [dir2,'Correlation_LR',num2str(LR),'.fig'];
saveas(fig1, fnamefig);
hold off;
%% Single correlation matrix

LR = 0.1;
LE = 1;
CR = 0.06;
CE = 0; 
nC = 11;
nL = 11;
time = 200;
dt = 0.01;
tLim = 100;
allstp = floor(time/dt);
LRall = linspace(0,0.1,nL);
dir1 = 'P:\FSE_MACSBIO\maryam.karimian\data\Timm_CondPlasticity_lowCR_InitT7_5th/';
dir2 = 'P:\FSE_MACSBIO\maryam.karimian\data\output\Timm_CondPlasticity_lowCR_InitT7_5th/Correlations\';

load([dir1,'LR',num2str(LR),'_LE',num2str(LE),'_CR',num2str(CR),'_CE',num2str(CE),'.mat']);
D = zeros(N,N);
l = size(Phase,2)-tLim;
for k = 1:tLim
    for i = 1:N
        for j = i:N
            D(i,j) = D(i,j) + cos(Phase(i,l)-Phase(j,l)) / tLim;
            D(j,i) = D(i,j);
        end
    end
    l = l + 1;
end
%-------------------------------
fig1 = figure;
imagesc(D);
colormap Jet;
caxis([-1 1]);
x = gca;
xlim([1 N])
ylim([1 N])

x.XDir = 'normal';
x.YDir = 'normal';
pbaspect([1 1 1])
x.XTick = [1:19:N];
x.YTick = [1:19:N];
x.XTickLabel = {'1' , '20', '40' , '60' , '80', '100'};
x.YTickLabel = {'1' , '20', '40' , '60' , '80', '100'};
set(gca,'FontSize',16);
titname = ['\eta=',num2str(CR),', \beta=',num2str(CE)];
title(titname);
fnametif = [dir2,'LR',num2str(LR),'CR',num2str(CR),'_CE',num2str(CE),'.tif'];
saveas(fig1, fnametif);
fnamefig = [dir2,'LR',num2str(LR),'CR',num2str(CR),'_CE',num2str(CE),'.fig'];
saveas(fig1, fnamefig);
%% Conduction matrix

LR = 0.1;
LE = 1;
N = 100;
time = 200;
dt = 0.01;
allstp = floor(time/dt);
nL = 11;
nC = 11;
CRall = linspace(0,0.1,nC);
CEall = linspace(0,1,nC);
tLim = 100;
dir1 = 'P:\FSE_MACSBIO\maryam.karimian\data\Timm_CondPlasticity_lowCR_InitT7_22d/';
dir2 = 'P:\FSE_MACSBIO\maryam.karimian\data\output\Timm_CondPlasticity_lowCR_InitT7_22d/Conductions\';
fig1 = figure;
for iCR = 1:numel(CRall)
    CR = CRall(iCR);
    for iCE = 1:numel(CEall)
        CE = CEall(iCE);
        load([dir1,'LR',num2str(LR),'_LE',num2str(LE),'_CR',num2str(CR),'_CE',num2str(CE),'.mat']);
%         loc = (iCR-1)*11 + iCE;
        loc = 11*11-( ( (iCR-1)*11 )+(11-iCE) );
        subplot(11,11,loc)
        imagesc(Cond2(:,:,end));
        colormap Jet;
%         colorbar;
        caxis([0 1]);
        x = gca;
        xlim([1 N])
        ylim([1 N])
        
        x.XDir = 'normal';
        x.YDir = 'normal';
        x.XTick = [1:49:N];
        x.YTick = [1:49:N];
        x.XTickLabel = {'1', '50', '100'};
        x.YTickLabel = {'1', '50', '100'};
        pbaspect([1 1 1])
    end
end
tightfig;
fig1 = gcf;
fig1.PaperUnits = 'inches';
fig1.PaperPosition = [0 0 20 10];
fnametif = [dir2,'Conduction_LR',num2str(LR),'.tif'];
saveas(fig1, fnametif);
fnamefig = [dir2,'Conduction_LR',num2str(LR),'.fig'];
saveas(fig1, fnamefig);
hold off;
%% Single conduction matrix

LR = 0.1;
LE = 1;
CR = 0.06;
CE = 0; 
nC = 11;
time = 200;
dt = 0.01;
tLim = 100;
allstp = floor(time/dt);
LRall = linspace(0,0.1,nL);
dir1 = 'P:\FSE_MACSBIO\maryam.karimian\data\Timm_CondPlasticity_lowCR_InitT7_5th/';
dir2 = 'P:\FSE_MACSBIO\maryam.karimian\data\output\Timm_CondPlasticity_lowCR_InitT7_5th/Conductions\';

load([dir1,'LR',num2str(LR),'_LE',num2str(LE),'_CR',num2str(CR),'_CE',num2str(CE),'.mat']);
fig1 = figure;
imagesc(Cond2(:,:,end));
colormap Jet;
caxis([0 1]);
x = gca;
xlim([1 N])
ylim([1 N])

x.XDir = 'normal';
x.YDir = 'normal';
pbaspect([1 1 1])
x.XTick = [1:19:N];
x.YTick = [1:19:N];
x.XTickLabel = {'1' , '20', '40' , '60' , '80', '100'};
x.YTickLabel = {'1' , '20', '40' , '60' , '80', '100'};
set(gca,'FontSize',16);
titname = ['\eta=',num2str(CR),', \beta=',num2str(CE)];
title(titname);
fnametif = [dir2,'LR',num2str(LR),'CR',num2str(CR),'_CE',num2str(CE),'.tif'];
saveas(fig1, fnametif);
fnamefig = [dir2,'LR',num2str(LR),'CR',num2str(CR),'_CE',num2str(CE),'.fig'];
saveas(fig1, fnamefig);
