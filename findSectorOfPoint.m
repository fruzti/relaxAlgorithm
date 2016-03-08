function sector = findSectorOfPoint(refPoint, newPoint)

    sector = 0;
    
    phi = pi + atan2(newPoint(2)-refPoint(2),newPoint(1)-refPoint(1));
    
    if (phi < 5*pi/4 && phi > 3*pi/4)
%         maskMax(xGrid > 1.5,:) = 0;
        sector = 1;
    elseif (phi > 5*pi/4 && phi < 7*pi/4)
%         maskMax(:,yGrid > 0.7) = 0;
        sector = 2;
    elseif (phi > 7*pi/4 || phi < pi/4)
%         maskMax(xGrid < 0.4,:) = 0;
        sector = 3;
    elseif (phi > pi/4 && phi < 3*pi/4)
%         maskMax(:,yGrid < 0) = 0;
        sector = 4;
    end

end