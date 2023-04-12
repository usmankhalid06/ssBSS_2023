clear;
close all; 
clc;


%% generic parameters
% sp = simtb_create_sP('exp_params_aod');
% sp.SM_spread =4+0.0001*randn(240,150*150);
% SM = simtb_makeSM(sp,1);   % Create spatial maps
% TC = zscore(simtb_makeTC(sp,1));  % Create TCs
load sources
N = size(TC,1);
nV = sqrt(size(SM,2));
nSRCS = size(TC,2);
K = 16; %dimensionality reduction sources
nIter = 30; %algorithm iterations
tstd  = sqrt(0.6); 
sstd  = sqrt(0.01); 
Dp = dctbases(N,N); %dct basis dictionary
nAlgos = 3;
rng('default'); rng('shuffle') % random number generator
Y= (TC+tstd*randn(N,nSRCS))*(SM+sstd*randn(nSRCS,nV^2));
Y= Y-repmat(mean(Y),size(Y,1),1);
nA = 3;

%% ssBSSP
tic;
params1.K = 16;
params1.P = 8; 
params1.lam1 = 0.01; 
params1.lam2 = 0.01; 
params1.lam3 = 12; 
params1.zeta = 60;
params1.Kp = 120;
params1.nIter = nIter;
params1.algo = 'PC';
params1.upd = 'blk';
params1.alpha = 10^-8;
[D{1},X{1},E(:,1)]=ssBSS(Y,Dp,params1);
D{1} = D{1}';
toc

%% ssBSSA
tic;
params2.K = 250;
params2.P = 8;
params2.lam1 = 2; 
params2.lam2 = 2; 
params2.lam3 = 12; 
params2.zeta = 60;
params2.Kp = 120;
params2.nIter = nIter;
params2.algo = 'AE';
params2.upd = 'blk';
params2.alpha = 10^-8; 
[D{2},X{2},E(:,2)]=ssBSS(Y,Dp,params2);
D{2}= D{2}';
toc

%% ssBSSAS
tic;
params3.K = 250; 
params3.P = 8;
params3.lam1 = 2; 
params3.lam2 = 2; 
params3.lam3 = 12; 
params3.zeta = 60;
params3.Kp = 120;
params3.nIter = nIter;
params3.algo = 'AE';
params3.upd = 'seq';
params3.alpha = 10^-8;
[D{3},X{3},E(:,3)]=ssBSS(Y,Dp,params3);
D{3}= D{3}';
toc

%% plots
for jj =1:nAlgos
    [sD{jj},sX{jj},ind]=sort_TSandSM_spatial(TC,SM,D{jj},X{jj},nSRCS);
    for ii =1:nSRCS
        TCcorr(jj,ii) =abs(corr(TC(:,ii),D{jj}(:,ind(ii))));
        SMcorr(jj,ii) =abs(corr(SM(ii,:)',X{jj}(ind(ii),:)'));
    end
end
cTC = sum(TCcorr')
cSM = sum(SMcorr')

my_subplots_just3(nA,nSRCS,nV,nV,TCcorr,SMcorr,sD,sX)
