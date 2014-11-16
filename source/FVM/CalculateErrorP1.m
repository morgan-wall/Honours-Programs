function [error] = CalculateErrorP1(numericDiag, analyticDiag)
error = norm(numericDiag - analyticDiag, inf);
end