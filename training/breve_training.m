% Ligature recognition
% Modified from acfDemoCal.m
% See also acfReadme.m
%
% Piotr's Computer Vision Matlab Toolbox      Version 3.40
% Copyright 2014 Piotr Dollar.  [pdollar-at-gmail.com]
% Licensed under the Simplified BSD License [see external/bsd.txt]

%DEFINE where Piotr's toolbox installed
piotr = '/path/to/piotr/toolbox/'
addpath(genpath(piotr))
%% extract training and testing images and ground truth
cd(fileparts(which('acfDemoCal.m'))); 
dataDir='/path/to/training/note_breve/';
%{
for s=1:2
  if(s==1), type='test'; skip=[]; else type='train'; skip=4; end
  dbInfo(['Usa' type]); if(s==2), type=['train' int2str2(skip,2)]; end
%  if(exist([dataDir type '/annotations'],'dir')), continue; end
 % dbExtract([dataDir type],1,skip);
end
%}
%% set up opts for training detector (see acfTrain)
opts=acfTrain(); 
opts.modelDs=[25 25]; opts.modelDsPad=[30 30];
opts.pPyramid.pChns.pColor.smooth=0; opts.nWeak=[256];
opts.pBoost.pTree.maxDepth=5; opts.pBoost.discrete=0;
opts.pBoost.pTree.fracFtrs=1/16; opts.nNeg=900; opts.nAccNeg=900;
opts.pPyramid.pChns.pGradHist.softBin=1; opts.pJitter=struct('flip',1);
%opts.posGtDir=[dataDir 'train' int2str2(skip,2) '/annotations'];
%opts.posImgDir=[dataDir 'train' int2str2(skip,2) '/images'];
opts.posWinDir=[dataDir 'resized30'];
opts.negWinDir=[dataDir 'negative-cr30_color'];

opts.pPyramid.pChns.shrink=2; opts.name='models/breved6';
%breved5 --negative-cr30 (does not include negative samples of color-breve)
%pLoad={'lbls',{'person'},'ilbls',{'people'},'squarify',{3,.41}};
%opts.pLoad = [pLoad 'hRng',[50 inf], 'vRng',[1 1] ];

%% optionally switch to LDCF version of detector (see acfTrain)
if( 0 ), opts.filters=[5 4]; opts.name='models/LdcfCaltech'; end

%% train detector (see acfTrain)
detector = acfTrain( opts );

%% modify detector (see acfModify)
pModify=struct('cascThr',-1,'cascCal',.025);
detector=acfModify(detector,pModify);

%% run detector on a sample image (see acfDetect)
imgNms=bbGt('getFiles',{[dataDir 'test/images']});
I=imread(imgNms{3}); tic, bbs=acfDetect(I,detector); toc
figure(1); im(I); bbApply('draw',bbs); pause(.1);

%% test detector and plot roc (see acfTest)
%[~,~,gt,dt]=acfTest('name',opts.name,'imgDir',[dataDir 'test/images'],...
%  'gtDir',[dataDir 'test/annotations'],'pLoad',[pLoad, 'hRng',[50 inf],...
%  'vRng',[.65 1],'xRng',[5 635],'yRng',[5 475]],...
 % 'pModify',pModify,'reapply',0,'show',2);
