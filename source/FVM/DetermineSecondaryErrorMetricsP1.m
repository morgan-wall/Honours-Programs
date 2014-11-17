function [maxPeakError, samePeakLocation] = ...
    DetermineSecondaryErrorMetricsP1(numericSol, analyticSol)

numericSol = numericSol(:, end);
analyticSol = analyticSol(:, end);

numericMax = max(numericSol);
analyticMax = max(analyticSol);

samePeakLocation = find(numericSol == (numericMax)) ...
    == find(analyticSol == (analyticMax));

disp(['Is same peak location: ' num2str(samePeakLocation)]);
end