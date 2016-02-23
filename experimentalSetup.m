clear all
% Position in Space of the experimental setup
angStep = 2*pi/6;
elNum = 0:5;
r = .06;

% arrays
array1 = [-(.4 + 3*.8 + .75); 5*.8 + 2*.75];
ctrArray1 = array1 + [.03; -.1];
arrayEl1 = array1 + [0; -.05];

% array1 Channels and Positions
arrayElem1 = getArrayElemPos(arrayEl1,ctrArray1);
arrayChans1 = [35 36 37 38 39 40];

array2 = [-(.4 + 2*.8); .75 + 3*.8];
ctrArray2 = array2 + [-.08; -.03];
arrayEl2 = array2 + [-.02; 0];

% array2 Channels and Positions
arrayElem2 = getArrayElemPos(arrayEl2, ctrArray2);
arrayChans2 = [3 4 5 6 7 8];

array3 = [-0; .75];
ctrArray3 = array3 + [-.1; -.03];
arrayEl3 = array3 + [-0.05; 0];

% array3 Channels and Positions
arrayElem3 = getArrayElemPos(arrayEl3, ctrArray3);
arrayChans3 = [46 45 44 43 42 41];

arrays = [array1 array2 array3];
ctrArray = [ctrArray1 ctrArray2 ctrArray3];
arrayElems = [arrayElem1 arrayElem2 arrayElem3];

% references
singleRef = [-(.4 + 3*.8); .75 + 3*.8] + [-0.15; 0]; % correction dist
% singleRef = [-(.4 + 3*.8); .75 + 3*.8] + [-0.08; 0]; % correction dist
singleChan = 48;

% arrayRef = [-(.4 + 2*.8); .75] + [0.14; 0]; %correction of distance
arrayRef = [-(.4 + 2*.8); .75] + [0.09; 0]; %correction of distance
ctrArrayRef = arrayRef + [-0.03; 0.1];
tmpElem = getArrayElemPos(arrayRef + [0; .05],...
    ctrArrayRef);
pairRef = tmpElem(:,[2 5]);
pairChans = [1 2];
referenceChans = [singleChan pairChans];
% arrayRefElems

references = [singleRef arrayRef];

% loudspeakers
frontLoud = [-(.4 + 3*.8 + .75 + .5*.8); .5*.75];
ladLoud = [-(.4 + 3*.8 + .5*.75); .5*.75] + [.1; 0] + [0.03; 0];
loudSpeakers = [frontLoud ladLoud];

%%
% plot(arrays(1,:),arrays(2,:),'or')
% hold on, 
% plot(arrayElems(1,:), arrayElems(2,:), '.g')
% plot(ctrArray(1,:),ctrArray(2,:),'om')
% plot(references(1,:), references(2,:), 'xb')
% plot(singleRef(1),singleRef(2),'.m')
% plot(loudSpeakers(1,:),loudSpeakers(2,:), 'dg')
% plot(pairRef(1,:), pairRef(2,:), '.m')
% grid on