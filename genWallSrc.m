function wallSrcs = genWallSrc(srcPos, walls)

    numWalls = size(walls,1);
    wallSrcs = nan(numWalls,2);
    
    for w = 1:numWalls
        switch (walls(w,1))
            case 0
                wallSrcs(w,1) = srcPos(1);
                wallSrcs(w,2) = 2*walls(w,2) - srcPos(2);
            case 1
                wallSrcs(w,2) = srcPos(2);
                wallSrcs(w,1) = 2*walls(w,2) - srcPos(1);
        end
    end

end