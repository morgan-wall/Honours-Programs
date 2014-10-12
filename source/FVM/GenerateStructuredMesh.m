function [xNodes, yNodes] = GenerateStructuredMesh(xLower, xUpper, xCount, ...
    yLower, yUpper, yCount, xGeometricParameters, yGeometricParameters)

xNodes = Generate1DMesh(xCount, xLower, xUpper, xGeometricParameters);
yNodes = Generate1DMesh(yCount, yLower, yUpper, yGeometricParameters);
end

function nodes = Generate1DMesh(count, lowerBound, upperBound, ...
    geometricParameters)

midPoint = (upperBound - lowerBound) / 2 + lowerBound;
halfCount = ceil(count / 2);

lowerHalfNodes = ...
    GenerateNodeLocations(halfCount, lowerBound, midPoint, ...
    geometricParameters.lowerIsGeometric, geometricParameters);
upperHalfNodes = ...
    GenerateNodeLocations(halfCount, upperBound, midPoint, ...
    geometricParameters.upperIsGeometric, geometricParameters);

lowerHalfNodes = lowerHalfNodes(:);
upperHalfNodes = upperHalfNodes(:);

nodes = [lowerHalfNodes; upperHalfNodes(2:end)];
end

function nodeLocations = GenerateNodeLocations(count, lowerBound, ...
    upperBound, isGeometric, geometricParameters)

if (isGeometric)
    nodeLocations = GenerateGeometricProgression(count - 1, lowerBound, ...
        upperBound, geometricParameters.commonRatio);
    nodeLocations = [lowerBound; nodeLocations + lowerBound];
    nodeLocations = sort(nodeLocations);
else
    nodeLocations = linspace(lowerBound, upperBound, count);
    nodeLocations = sort(nodeLocations);
end
end

function geometricProg = GenerateGeometricProgression(count, lowerBound, ...
    upperBound, commonRatio)

scaleFactor = (upperBound - lowerBound) ...
    * (1 - commonRatio) / (1 - commonRatio^(count));

geometricProg = [1; cumprod( commonRatio .* ones(count - 1, 1))];
geometricProg = cumsum(geometricProg .* scaleFactor);
end