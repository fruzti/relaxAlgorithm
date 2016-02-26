function srcList = runMethodTest2(micFreqData,srcTimeData,srcFreqData,...
    micPos,fs,N,K,xGrid,yGrid)

% srcTestPos = [0.955 0.3750];
srcTestPos = [0.9859 0.3984];

[~, ~, ~, Jtotal] = nfEstML_TOA_DOA(micFreqData, srcFreqData,...
    micPos(:,1:2), fs);

L = 4;

Jp = cell(L,1);

sJp = cell(L,1);

smoothSize = [10,10];
sJtotal = medfilt2(Jtotal, smoothSize);
srcList = zeros(L,2);

sJpSum = 0;

y_p = nfGenDelayData(srcTimeData, micPos(:,1:2)', srcTestPos', fs);
yFreqData = getFreqMicData(y_p, N, K);
[~, ~, ~, JBase] = nfEstML_TOA_DOA(yFreqData, srcFreqData,...
    micPos(:,1:2), fs);
sJBase = medfilt2(JBase,smoothSize);

for l = 1:L
    
    j0 = l;
    itErr = 1;
    
    while(itErr > 1e-3)
        
        for j = j0:l
    
            indVec = 1:l;
            indVec = indVec(indVec ~= j);
            
            sJpSum = JBase;
            
            for indx = indVec
                sJpSum = sJpSum + Jp{indx};
            end
            
            sJpSum = medfilt2(sJpSum,smoothSize);
            
            [~,tmpIndx] = findMaximizer(sJtotal./sJpSum);
            srcList(j,:) = [xGrid(tmpIndx(1)) yGrid(tmpIndx(2))];
    
            % line search
            myCost = @(x) -nfEvalCost(micFreqData,srcFreqData, micPos(:,1:2)',...
                    x,floor(N/2),N,fs);
            x0 = [srcList(j,1); srcList(j,2)];
            options = optimset('Display', 'off') ;
            xStar = fmincon(myCost,x0,[],[],[],[],x0 - 1, x0 + 1,[],options);
            srcList(j,1) = xStar(1); srcList(j,2) = xStar(2);
    
            srcTestPos = srcList(j,:);
            
            y_p = nfGenDelayData(srcTimeData, micPos(:,1:2)', srcTestPos', fs);
            yFreqData = getFreqMicData(y_p, N, K);
            [~, ~, ~, Jp{j}] = nfEstML_TOA_DOA(yFreqData, srcFreqData,...
                micPos(:,1:2), fs);
            sJp{j} = medfilt2(Jp{j},smoothSize);
            
        end
    
   j0 = 1;
        
   currCost = nfGetCost(applyIFFT(micFreqData,N), applyIFFT(srcFreqData,N),...
       micPos(:,1:2)',...
       [0.9859; srcList(:,1)], [0.3984; srcList(:,2)], ones(l+1,1), l+1, fs);
        
   if (l == 1)
       prevCost = currCost;
   end
   
   % RELAX Mode
   itErr = abs( currCost - prevCost );

   % SEQUENTIAL Mode
   % itErr = 0;
        
   prevCost = currCost;
   
   end
    
end

end