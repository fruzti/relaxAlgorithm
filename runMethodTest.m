function srcList = runMethodTest(micFreqData,srcTimeData,srcFreqData,...
    micPos,fs,N,K,xGrid,yGrid)

% srcTestPos = [0.955 0.3750];
srcTestPos = [0.9859 0.3984];

[~, ~, ~, Jtotal] = nfEstML_TOA_DOA(micFreqData, srcFreqData,...
    micPos(:,1:2), fs);

L = 15;

Jp = cell(L,1);

sJp = cell(L,1);

smoothSize = [10,10];
sJtotal = medfilt2(Jtotal, smoothSize);
srcList = zeros(L,2);

sJpSum = 0;


% feasible region
[r2Src,~,~] = getMic2SrcParams(micPos(:,1:2)',srcTestPos');
phi = 0:0.1:2*pi;
for it = 1:size(r2Src,1)
    pntsCirc(:,it,:) = r2Src(it) * [cos(phi)' sin(phi)'] + ...
        repmat(micPos(it,1:2),length(phi),1);
end

pntsCirc = [vec(pntsCirc(:,:,1)),vec(pntsCirc(:,:,2))];
[xM, yM] = meshgrid(xGrid,yGrid); xM = xM(:); yM = yM(:);
feasiblePoints = inpolygon(xM,yM,pntsCirc(:,1), pntsCirc(:,2));
maskFeasible = reshape(~feasiblePoints, length(yGrid), length(xGrid))';
%
for l = 1:L
    
    y_p = nfGenDelayData(srcTimeData, micPos(:,1:2)', srcTestPos', fs);
    yFreqData = getFreqMicData(y_p, N, K);
    [~, ~, ~, Jp{l}] = nfEstML_TOA_DOA(yFreqData, srcFreqData,...
        micPos(:,1:2), fs);
    sJp{l} = medfilt2(Jp{l},smoothSize);
    
    sJpSum = sJpSum + sJp{l};
    
    [~,tmpIndx] = findMaximizer(sJtotal./sJpSum,maskFeasible);
    srcList(l,:) = [xGrid(tmpIndx(1)) yGrid(tmpIndx(2))];
    
    % line search
    myCost = @(x) -nfEvalCost(micFreqData,srcFreqData, micPos(:,1:2)',...
        x,floor(N/2),N,fs);
    x0 = [srcList(l,1); srcList(l,2)];
    options = optimset('Display', 'off') ;
    xStar = fmincon(myCost,x0,[],[],[],[],x0 - 1, x0 + 1,[],options);
    srcList(l,1) = xStar(1); srcList(l,2) = xStar(2);
    
    srcTestPos = srcList(l,:);
    
end

end