function [maxPeakError, samePeakLocation] = ...
    DetermineSecondaryErrorMetricsP1(numericSol, analyticSol)

numericMax = max(max(numericSol));
analyticMax = max(max(analyticSol));

maxPeakError = abs(analyticMax - numericMax) / abs(analyticMax);

samePeakLocation = find(numericSol == (numericMax)) ...
    == find(numericSol == (numericMax));

disp(['Max peak (rel) error: ' num2str(maxPeakError)]);
disp(['Is same peak location: ' num2str(samePeakLocation)]);
end