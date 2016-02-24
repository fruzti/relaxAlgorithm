function vrtSrcPos = genSrcsFromWalls(srcPos, rightCorner)

    wallOne = [-srcPos(1) ; srcPos(2)];
    wallTwo = [srcPos(1) ; -srcPos(2)];
    wallThree = [srcPos(1) + 2*(rightCorner(1) - srcPos(1)) ; srcPos(2)];
    wallFour = [srcPos(1) ; srcPos(2) + 2*(rightCorner(2) - srcPos(2))];
    
    vrtSrcPos = [wallOne wallTwo wallThree wallFour];

end