function error = findSmallestAngularError(expAngle, actAngle)
    candidateErrors = actAngle-expAngle-[-2*pi;0;2*pi];
    [~,idx] = min(candidateErrors.^2);
    error = candidateErrors(idx);
end

