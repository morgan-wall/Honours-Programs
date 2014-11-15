function [tout, yout, gmresIterations, nonlinearFnCalls] = Solver(dt, tFinal, ...
    Dxx, Dyy, Vx, Vy, source, theta, advectionHandling, nodesX, nodesY, ...
    northBC, eastBC, southBC, westBC, initialCondition, storedTimeSteps, ...
    newtonParameters, gmresParameters, forcingTermParameters, ...
    safeguardParameters, chordSteps, isGlobalised, linesearchParam, ...
    minLambda, maxLambda, maxBacktracks)
%% Solver: solve non-linear two-dimesional advection-diffusion equation.
% Determine a numerical solution for the non-linear, two-dimensional
% advection-diffusion equation given by:
%
%   $\frac{\partial \phi}{\partial t} + \Delta (\mathbf{V} \phi - \mathbf{D}
%   \Delta \phi) = S_u$
%
% where $\Delta(\phi) = diag(D_xx, D_yy)$, $\mathbf{V}(\phi) = (v_x,
% v_y)^T$, and $S_u(\phi)$.
%
% A vertex-centred finite volume strategy is used on a structured, regular
% mesh. A right preconditioned inexact Newton-GMRES solver is used to
% obtain solutions for the nonlinear system generated at each time step.
% The supported preconditioners include Jacobi, SOR, and ILU.
%
% The inexact Newton-GMRES solver uses the modified Eisenstat-Walker
% formula for determining the forcing term. The formula can be found in:
%
%   Eisenstat, S. C., & Walker, H. F. (1996). Choosing the forcing terms in
%   an inexact Newton Method. SIAM J. Sci. Comput, 17(1), 16-32.
%   doi: 10.1137/09170
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation notes:
%
%   This solver requires boundary conditions for the four faces defining
%   the 'surface' of the 2D region of interest. The form of the boundary
%   conditions are given by:
%
%       A_i \phi + B_i \frac{\partial \phi}{\partial t} = C_i.
%
%   The notation adopted for coefficients A, B, and C is used throughout
%   this solver.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%
%   dt:
%       The delta between time steps.
%
%   tFinal:
%       The final time for solutions to be determined. The PDE is solved
%       from 0 < t < tFinal using an adaptive time-stepping strategy.
%
%   Dxx:
%       A function handle that takes one parameter: the solution value at
%       a point. The function must return the horizontal diffusion rate.
%
%   Dyy:
%       A function handle that takes one parameter: the solution value at
%       a point. The function must return the vertical diffusion rate.
%
%   Vx:
%       A function handle that takes one parameter: the solution value at
%       a point. The function must return the horizontal advection rate.
%
%   Vy:
%       A function handle that takes one parameter: the solution value at
%       a point. The function must return the vertical advection rate.
%
%   source:
%       A function handle that takes one parameter: the solution value at
%       a point. The function must return the source (rate).
%
%   theta:
%       The temporal weighting (0 for forward Euler, 1 for backward Euler,
%       and 1/2 for Crank-Nicolson).
%
%   advectionHandling:
%       A string denoting the method used for approximating the advection
%       at a control volume face. Valid arguments include 'averaging' and
%       'upwinding'.
%
%   nodesX:
%       An array containing the unique x-coordinates of nodes in the 2D 
%       mesh (in ascending order).
%
%   nodesY:
%       An array containig the unique y-coordinates of nodes in the 2D 
%       mesh (in descending order).
%
%   northBC:
%       A struct specifying the coefficients for the boundary condition
%       defined on the north (i.e. upper) face of the region. The fields
%       are A, B, and C, which correspond to the coefficients in the
%       generalised boundary condition (see implementation notes). The
%       value corresponding to B must be positive.
%
%   eastBC:
%       A struct specifying the coefficients for the boundary condition
%       defined on the east (i.e. right) face of the region. The fields
%       are A, B, and C, which correspond to the coefficients in the
%       generalised boundary condition (see implementation notes). The
%       value corresponding to B must be positive.
%
%   southBC:
%       A struct specifying the coefficients for the boundary condition
%       defined on the south (i.e. lower) face of the region. The fields
%       are A, B, and C, which correspond to the coefficients in the
%       generalised boundary condition (see implementation notes). The
%       value corresponding to B must be positive.
%
%   westBC:
%       A struct specifying the coefficients for the boundary condition
%       defined on the west (i.e. left) face of the region. The fields
%       are A, B, and C, which correspond to the coefficients in the
%       generalised boundary condition (see implementation notes). The
%       value corresponding to B must be positive.
%
%   initialCondition:
%       A matrix containing the initial solution value at each node in 
%       the mesh.
%
%   storedTimeSteps:
%       Solutions at time steps that are multiples of this (integer) 
%       parameter are returned from this function.
%
%   newtonParameters:
%       A struct specifying parameters for the inexact Newton solver used
%       for solving the non-linear system of equations generated by the
%       FVM at each time step. The fields include:
%           maxIterations:
%               The maximum number of iterations of the inexact Newton
%               method.
%           tolUpdate:
%               The tolerance applied to the Newton update to determine
%               the stopping point of the inexact Newton method.
%           tolResidual:
%               The tolerance applied to the nonlinear residual to 
%               determine the stopping point of the inexact Newton method.
%
%   gmresParameters:
%       A struct specifying parameters for the GMRES solver used for
%       solving the linear systems in the Newton step. The fields include:
%           maxIterations:
%               The maximum number of iterations to perform of the GMRES
%               method.
%           restartValue:
%               The iteration count at which a restart of the GMRES method
%               occurs.
%           errorTol:
%               The relative error tolerance applied to the GMRES method.
%           preconditioningType:
%               A string denoting the type of preconditioning used in the
%               GMRES method. Only right preconditioning is supported. The
%               valid inputs include 'jacobi', 'SOR', and 'ilu'.
%           omega:
%               The value of the omega parameter used in SOR
%               preconditioning. This value is only used if SOR
%               preconditioning is used.
%
%   forcingTermParameters:
%       A struct specifying parameters for the forcing term used in the 
%       inexact Newton solver. The fields include:
%           maxForcingTerm:
%               The upper bound on the value of the forcing term $n_k$.
%               This parameter must have a value between zero and one,
%               inclusive.
%           type:
%               A string denoting the type of forcing term used. The valid
%               inputs include 'choice1' and 'choice2', which correspond to
%               the two forcing terms proposed by Eisenstat & Walker
%               (1996).
%           gamma:
%               A double in [0, 1] used in calculating the forcing term.
%               This parameter is only required if type is 'choice2'.
%           alpha:
%               A double in (1, 2] used in calculating the forcing term.
%               This parameter is only required if type is 'choice2'.
%
%   safeguardParameters:
%       A struct containing parameters for the safeguards imposed on the 
%       forcing terms in the inexact Newton method. The fields include:
%           threshold:
%               The minimum value of the forcing term to which safeguarding
%               is applied. If a proposed forcing term is less than or
%               equal to this value, then safeguarding is not imposed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
%   
%   tout:
%       An array containing the times corresponding to the solutions 
%       contained in yout.
%
%   yout:
%       A 3D matrix containing the solution of the system at each node in
%       the mesh for each time specified in tout. In terms of indexing
%       notation, the first and second indices correspond to the node in
%       the mesh. The third index corresponds to the tout variable
%       specifying solution times. For example: yout(:, :, 1) retrieves the
%       matrix of solution values at the first time step.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Errors:
%   
%   1. Errors are thrown if the B field in any boundary condition struct 
%      is zero.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unimplemented features:
%
%   1. Jacobian-free Newton-GMRES solver.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialise constants

FIRST_TIME_STEP_INDEX = 1;

%% Initialise mesh parameters

nodesY = nodesY(:);
nodesX = nodesX(:);

rows = length(nodesY);
columns = length(nodesX);
nodeCount = rows * columns;

% Absolute distances between adjacent nodes
xNodeDeltas = diff(nodesX);
yNodeDeltas = abs(diff(nodesY));

% Absolute distances between nodes and CV faces
xFaceDeltas = xNodeDeltas ./ 2;
yFaceDeltas = yNodeDeltas ./ 2;

% Control volumne dimensions
nodeWidths = zeros(columns, 1);
nodeWidths(1) = xFaceDeltas(1);
nodeWidths(2:end-1) = xFaceDeltas(2:end) + xFaceDeltas(1:end-1);
nodeWidths(end) = xFaceDeltas(end);

nodeHeights = zeros(rows, 1);
nodeHeights(1) = yFaceDeltas(1);
nodeHeights(2:end-1) = yFaceDeltas(2:end) + yFaceDeltas(1:end-1);
nodeHeights(end) = yFaceDeltas(end);

%% Modify mesh parameters for optimal vectorised code
% N.B. This improves efficiency at the expense of increased memory requirements.

indices = 1:nodeCount;

% Store the row and column for each index
rowForIndex = mod(indices - 1, rows) + 1;
columnForIndex = floor((indices - 1) ./ rows) + 1;

nodesXPos = nodesX(columnForIndex(indices));
nodesYPos = nodesY(rowForIndex(indices));
nodeWidthsAll = nodeWidths(columnForIndex(indices));
nodeHeightsAll = nodeHeights(rowForIndex(indices));

% Store indices that lie on each boundary
isNBoundaryIndex = mod(indices - 1, rows) == 0;
isEBoundaryIndex = indices > nodeCount - rows;
isSBoundaryIndex = mod(indices, rows) == 0;
isWBoundaryIndex = indices <= rows;

%% Initialise solution parameters

tStart = 0;
times = tStart:dt:tFinal;
timeSteps = length(times) - 1;

numStoredSolutions = floor(timeSteps / storedTimeSteps) + 1;
manuallyStoreFinalTimeSol = rem(timeSteps, storedTimeSteps) ~= 0;
if (manuallyStoreFinalTimeSol)
    numStoredSolutions = numStoredSolutions + 1;
end

yout = zeros(nodeCount, numStoredSolutions);
tout = zeros(numStoredSolutions, 1);

yout(:, 1) = initialCondition(:);
previousSolution = initialCondition(:);

%% Initialise performance metrics

nonlinearFnCalls = 0;
gmresIterations = 0;

%% Output conservation metrics

total = sum(yout(:, 1) ...
    .* (nodeHeights(rowForIndex(indices)) .* nodeWidths(columnForIndex(indices))));
disp(['Begin "mass": ' num2str(total)]);

%% Initialise solver parameters

isUpwinding = strcmp(advectionHandling, 'upwinding');

tauA = 1e-6;
tauR = 1e-6;

%% Iteratively solve advection-diffusion equation (time marching strategy)

for i = 1:timeSteps
    
    % Formulate the Forward Euler component of F(u) = 0
    F_forwardEuler = zeros(nodeCount, 1);
    if (theta ~= 1)
        F_forwardEuler = dt .* (1 - theta) .* GenerateFluxVec(...
            previousSolution, indices, rows, Vx, Vy, Dxx, Dyy, nodeWidths, ...
            nodeHeights, rowForIndex, columnForIndex, xNodeDeltas, ...
            yNodeDeltas, northBC, eastBC, southBC, westBC, isUpwinding, ...
            isNBoundaryIndex, isEBoundaryIndex, isSBoundaryIndex, ...
            isWBoundaryIndex, nodesX, nodesY, times(i+1), nodesXPos, ...
            nodesYPos, nodeWidthsAll, nodeHeightsAll);
        F_forwardEuler = F_forwardEuler ...
            - dt .* (1 - theta) .* source(previousSolution, nodesXPos, ...
                nodesYPos, times(i+1));
    end
    F_forwardEuler = F_forwardEuler - previousSolution;
    
    % Formulate the Backward Euler component of F(u) = 0 
    if (i == FIRST_TIME_STEP_INDEX)
        F_current_backwardEuler = dt .* theta .* GenerateFluxVec(...
            previousSolution, indices, rows, Vx, Vy, Dxx, Dyy, nodeWidths, ...
            nodeHeights, rowForIndex, columnForIndex, xNodeDeltas, ...
            yNodeDeltas, northBC, eastBC, southBC, westBC, isUpwinding, ...
            isNBoundaryIndex, isEBoundaryIndex, isSBoundaryIndex, ...
            isWBoundaryIndex, nodesX, nodesY, times(i+1), nodesXPos, ...
            nodesYPos, nodeWidthsAll, nodeHeightsAll);
        F_current_backwardEuler = F_current_backwardEuler ...
            - dt .* theta .* source(previousSolution, nodesXPos, ...
                nodesYPos, times(i+1));
        F_current_backwardEuler = F_current_backwardEuler + previousSolution;
        nonlinearFnCalls = nonlinearFnCalls + 1;
    end
    
    Fx = F_forwardEuler + F_current_backwardEuler;
    Fx_previous = Fx;
    Fx_initial = Fx;
    
    % Pre-emptively generate the Jacobian
    if (i == FIRST_TIME_STEP_INDEX)        
        [jacobian, jacobianNonlinearFnCalls] = GenerateJacobian(indices, ...
            nodeCount, previousSolution, dt, theta, rows, columns, ...
            rowForIndex, columnForIndex, Vx, Vy, Dxx, Dyy, xNodeDeltas, ...
            yNodeDeltas, nodeWidths, nodeHeights, northBC, eastBC, southBC, ...
            westBC, isUpwinding, source, F_forwardEuler, Fx, ...
            isNBoundaryIndex, isEBoundaryIndex, isSBoundaryIndex, ...
            isWBoundaryIndex, nodesX, nodesY, times(i+1), previousSolution, ...
            nodesXPos, nodesYPos, nodeWidthsAll, nodeHeightsAll);
        nonlinearFnCalls = nonlinearFnCalls + jacobianNonlinearFnCalls;
        previousJacobian = jacobian;
        
        [LPrecond, UPrecond] = GenerateMatrixPreconditioner(jacobian, ...
            gmresParameters.preconditioningType, gmresParameters.omega);
    end
    
    % Solve the non-linear system, F(x) = 0, for the next time step
    currentSolution = previousSolution;
    forcingTerm = 1/2;
    delta_x = realmax;
    current_iteration = 0;
    while (current_iteration <= newtonParameters.maxIterations ...
            && norm(delta_x) >= newtonParameters.tolUpdate * norm(currentSolution))
        
        if (current_iteration ~= 0 && mod(current_iteration, chordSteps) == 0)
            previousJacobian = jacobian;
            
            [jacobian, jacobianNonlinearFnCalls] = GenerateJacobian(indices, ...
                nodeCount, currentSolution, dt, theta, rows, columns, ...
                rowForIndex, columnForIndex, Vx, Vy, Dxx, Dyy, xNodeDeltas, ...
                yNodeDeltas, nodeWidths, nodeHeights, northBC, eastBC, ...
                southBC, westBC, isUpwinding, source, F_forwardEuler, Fx, ...
                isNBoundaryIndex, isEBoundaryIndex, isSBoundaryIndex, ...
                isWBoundaryIndex, nodesX, nodesY, times(i+1), ...
                previousSolution, nodesXPos, nodesYPos, nodeWidthsAll, ...
                nodeHeightsAll);
            nonlinearFnCalls = nonlinearFnCalls + jacobianNonlinearFnCalls;
            
            % Update preconditioner for Jacobian
            [LPrecond, UPrecond] = GenerateMatrixPreconditioner(jacobian, ...
                gmresParameters.preconditioningType, gmresParameters.omega);
        end
        
        % Determine the error tolerance based on the forcing term
        if (strcmp(forcingTermParameters.type, 'none'))
            residualError = gmresParameters.errorTol;
        else
            isInitialIteration = current_iteration == 0;
            residualError = DetermineErrorToleranceFromForcingTerm( ...
                forcingTerm, forcingTermParameters, safeguardParameters, Fx, ...
                Fx_previous, Fx_initial, previousJacobian, previousSolution, ...
                tauA, tauR, isInitialIteration);
        end
        
        % solve the linear system using GMRES
        [delta_x, iterations] = gmres_general(jacobian, Fx, currentSolution, ...
            LPrecond, UPrecond, gmresParameters.maxIterations, ...
            gmresParameters.restartValue, residualError);
        
        gmresIterations = gmresIterations + iterations;
        
        preGlobalisationSolution = currentSolution - delta_x;
        
        % Evaluate the nonlinear system for updated iterate
        F_current_backwardEuler = dt .* theta .* GenerateFluxVec(...
            preGlobalisationSolution, indices, rows, Vx, Vy, Dxx, Dyy, ...
            nodeWidths, nodeHeights, rowForIndex, columnForIndex, ...
            xNodeDeltas, yNodeDeltas, northBC, eastBC, southBC, westBC, ...
            isUpwinding, isNBoundaryIndex, isEBoundaryIndex, ...
            isSBoundaryIndex,isWBoundaryIndex, nodesX, nodesY, times(i+1), ...
            nodesXPos, nodesYPos, nodeWidthsAll, nodeHeightsAll);
        F_current_backwardEuler = F_current_backwardEuler ...
            - dt .* theta .* source(preGlobalisationSolution, nodesXPos, ...
                nodesYPos, times(i+1));
        F_current_backwardEuler = ...
            F_current_backwardEuler + preGlobalisationSolution;
        nonlinearFnCalls = nonlinearFnCalls + 1;
        
        Fx_previous = Fx;
        Fx = F_current_backwardEuler + F_forwardEuler;
        
        % Apply backtracking globalisation strategy
        if (isGlobalised)
            [Fx, delta_x, residualError, fCalls] = GlobaliseSolution(...
                currentSolution, residualError, delta_x, linesearchParam, ...
                minLambda, maxLambda, maxBacktracks, Fx, Fx_previous, dt, ...
                theta, source, F_forwardEuler, indices, rows, Vx, Vy, Dxx, ...
                Dyy, nodeWidths, nodeHeights, rowForIndex, columnForIndex, ...
                xNodeDeltas, yNodeDeltas, northBC, eastBC, southBC, westBC, ...
                isUpwinding, isNBoundaryIndex, isEBoundaryIndex, ...
                isSBoundaryIndex, isWBoundaryIndex, nodesX, nodesY, times(i+1));
            currentSolution = currentSolution - delta_x;
            nonlinearFnCalls = nonlinearFnCalls + fCalls;
        else
            currentSolution = preGlobalisationSolution;
        end
        
        current_iteration = current_iteration + 1;
    end
    
    % Ensure a solution
    if (current_iteration > newtonParameters.maxIterations ...
            && norm(delta_x) >= newtonParameters.tolUpdate * norm(currentSolution))
        error(['Method Failure: the non-linear system generated at '...
            't = ' num2str(i * dt) ' was not solved using ' ...
            num2str(newtonParameters.maxIterations) ' iterations of inexact '...
            'Newton method.']);
    end
    
    % Optionally store the solution
    if (mod(i, storedTimeSteps) == 0)
        lastStoredSolutionIndex = i / storedTimeSteps + 1;
        yout(:, lastStoredSolutionIndex) = currentSolution;
        tout(lastStoredSolutionIndex) = i * dt;
    elseif (i == timeSteps && manuallyStoreFinalTimeSol)
        yout(:, end) = currentSolution;
        tout(end) = i * dt;
    end
    
    previousSolution = currentSolution;
end

%% Output conservation metrics 

total = sum(yout(:, end) ...
    .* (nodeHeights(rowForIndex(indices)) .* nodeWidths(columnForIndex(indices))));
disp(['End "mass": ' num2str(total)]);

end

%
%   Helper Functions
%

function [Fx, deltaX, etaK, nonlinearFnCalls] = GlobaliseSolution(...
    previousSolution, etaK, deltaX, t, minLambda, maxLambda, maxIterations, ...
    Fx, Fx_previous, dt, theta, source, F_forwardEuler, indices, rows, Vx, ...
    Vy, Dxx, Dyy, nodeWidths, nodeHeights, rowForIndex, columnForIndex, ...
    xNodeDeltas, yNodeDeltas, northBC, eastBC, southBC, westBC, isUpwinding, ...
    isNBoundaryIndex, isEBoundaryIndex, isSBoundaryIndex, isWBoundaryIndex, ...
    nodesX, nodesY, time)

simpleLinesearchCoef = 2;

nonlinearFnCalls = 0;
updatedDeltaX = deltaX;
lambda = maxLambda * simpleLinesearchCoef;
iteration = 1;

while (norm(Fx) > (1 - t * (1 - etaK)) * norm(Fx_previous) ...
        && iteration < maxIterations)
    
    % Backtrack the delta
    lambda = lambda / simpleLinesearchCoef;
    if (lambda < minLambda) 
        lambda = minLambda;
    end
    
    updatedDeltaX = lambda .* updatedDeltaX;
    etaK = 1 - lambda * (1 - etaK);
    testSolution = previousSolution - updatedDeltaX;
    
    % Re-evaluate the non-linear system for the new delta
    F_current_backwardEuler = dt .* theta .* GenerateFluxVec(testSolution, ...
        indices, rows, Vx, Vy, Dxx, Dyy, nodeWidths, nodeHeights, ...
        rowForIndex, columnForIndex, xNodeDeltas, yNodeDeltas, northBC, ...
        eastBC, southBC, westBC, isUpwinding, isNBoundaryIndex, ...
        isEBoundaryIndex, isSBoundaryIndex, isWBoundaryIndex, nodesX, nodesY, ...
        time, nodesXPos, nodesYPos, nodeWidthsAll, nodeHeightsAll);
    F_current_backwardEuler = F_current_backwardEuler ...
        - dt .* theta .* source(testSolution, nodesXPos, nodesYPos, time);
    F_current_backwardEuler = F_current_backwardEuler + testSolution;
    
    Fx = F_current_backwardEuler + F_forwardEuler;
    
    nonlinearFnCalls = nonlinearFnCalls + 1;
    iteration = iteration + 1;
end

% Ensure a valid iterate was found
if (norm(Fx) > (1 - t * (1 - etaK)) * norm(Fx_previous) && lambda == minLambda ...
        && iteration >= maxIterations)
    error(['Method Failure: backtracking was unable to find an iterate to ' ...
        'satisfy the specified requirements on the residual norm (t = ' ...
        num2str(time) '). Consider increasing the maximum number of iterations' ...
        ' used to determine a new iterate.']);
end
end

function [jacobian, nonlinearFnCalls] = GenerateJacobian(indices, nodeCount, ...
    currentSolution, dt, theta, rows, columns, rowForIndex, columnForIndex, ...
    Vx, Vy, Dxx, Dyy, xNodeDeltas, yNodeDeltas, nodeWidths, nodeHeights, ...
    northBC, eastBC, southBC, westBC, isUpwinding, source, ...
    forwardEulerComponent, Fx, isNBoundaryIndex, isEBoundaryIndex, ...
    isSBoundaryIndex, isWBoundaryIndex, nodesX, nodesY, t, previousSolution, ...
    nodesXPos, nodesYPos, nodeWidthsAll, ...
    nodeHeightsAll)

nonZerosPerInternalNode = 5;
jacobianNonZeroCount = nodeCount * nonZerosPerInternalNode ...
    - (2 * rows) - (2 * columns);
jacobian = sparse([], [], [], nodeCount, nodeCount, jacobianNonZeroCount);

h = determine_newton_step_delta(currentSolution);

nonlinearFnCalls = 0;

initialOffset = 2;
offSet = rows + initialOffset;
for i = 0:rows+initialOffset
    
    % Perturb solutions at selected nodes
    initialIndex = 1 + i;
    perturbedNodeIndices = initialIndex:offSet:nodeCount;
    xStepped = currentSolution;
    xStepped(perturbedNodeIndices) = xStepped(perturbedNodeIndices) + h;
    
    % Determine impacted nodes
    impactedNodePerturbationNodeIndex = ...
        [perturbedNodeIndices perturbedNodeIndices perturbedNodeIndices ...
        perturbedNodeIndices perturbedNodeIndices];
    impactedNodes = ...
        [perturbedNodeIndices perturbedNodeIndices+1 perturbedNodeIndices-1 ...
        perturbedNodeIndices+rows perturbedNodeIndices-rows];
    isInvalidNode = impactedNodes <= 0 | impactedNodes > nodeCount;
    impactedNodePerturbationNodeIndex(isInvalidNode) = [];
    impactedNodes(isInvalidNode) = [];
    
    % Evaluate F(u + e) = 0
    F_backwardEuler_stepped = dt * theta * GenerateFluxVec(xStepped, ...
        impactedNodes, rows, Vx, Vy, Dxx, Dyy, nodeWidths, nodeHeights, ...
        rowForIndex, columnForIndex, xNodeDeltas, yNodeDeltas, northBC, ...
        eastBC, southBC, westBC, isUpwinding, isNBoundaryIndex, ...
        isEBoundaryIndex, isSBoundaryIndex, isWBoundaryIndex, ...
        nodesX, nodesY, t, nodesXPos, nodesYPos, nodeWidthsAll, ...
    nodeHeightsAll);
    F_backwardEuler_stepped = F_backwardEuler_stepped ...
        - dt * theta * source(previousSolution, nodesXPos(impactedNodes), ...
            nodesYPos(impactedNodes), t);
    F_backwardEuler_stepped = F_backwardEuler_stepped + xStepped(impactedNodes);
    nonlinearFnCalls = nonlinearFnCalls + 1;
    
    F_stepped = F_backwardEuler_stepped + forwardEulerComponent(impactedNodes);
    
    for j = 1:length(impactedNodes)
        jacobian(impactedNodes(j), impactedNodePerturbationNodeIndex(j)) = ...
            (F_stepped(j) - Fx(impactedNodes(j))) / h;
    end
end
end

function [L, U] = GenerateMatrixPreconditioner(A, type, omega)

% convert system parameters to sparse storage
if (~issparse(A))
    A = sparse(A); 
end

% construct matrix preconditioner
diagonal_index = 0;
matrixSize = size(A);
rows = matrixSize(1);

switch(type)
    case 'ilu'
        setup.type = 'nofill';
        setup.milu = 'row';
        setup.droptol = 1e-4;
        [L,U] = ilu(A, setup);
    case 'jacobi'
        L = spdiags(diag(A), diagonal_index, speye(rows));
        U = speye(rows);
    case 'sor'
        L = spdiags(diag(A) / omega, diagonal_index, tril(A));
        U = speye(rows);
    case 'none'
        L = speye(rows);
        U = L;
end
end

function residualError = DetermineErrorToleranceFromForcingTerm( ...
    forcingTerm, forcingTermParameters, safeguardParameters, Fx, ...
    Fx_previous, Fx_initial, previousJacobian, previousSolution, ...
    tauA, tauR, isInitialIteration)

switch(forcingTermParameters.type)
    case 'choice1'
        safeguardedTerm = forcingTerm^((1 + sqrt(5))/2);
        forcingTerm = ...
            norm(Fx - Fx_previous - previousJacobian * previousSolution) ...
            / norm(Fx_previous);
        if (safeguardedTerm > safeguardParameters.threshold)
            forcingTerm = max(forcingTerm, safeguardedTerm);
            forcingTerm = ...
                min(forcingTerm, forcingTermParameters.maxForcingTerm);
        end
        residualError = forcingTerm * norm(Fx);
    case 'choice2'
        safeguardedTerm = forcingTermParameters.gamma ...
            * forcingTerm^forcingTermParameters.alpha;
        forcingTerm = forcingTermParameters.gamma ...
            * (norm(Fx) / norm(Fx_previous))^forcingTermParameters.alpha;
        if (safeguardedTerm > safeguardParameters.threshold)
            forcingTerm = max(forcingTerm, safeguardedTerm);
            forcingTerm = ...
                min(forcingTerm, forcingTermParameters.maxForcingTerm);
        end
        residualError = forcingTerm * norm(Fx);
    case 'assignment'
        safeguardedTerm = forcingTermParameters.gamma ...
            * forcingTerm^forcingTermParameters.alpha;
        eta_k_R = forcingTermParameters.gamma ...
            * (norm(Fx, Inf) / norm(Fx_previous, Inf))...
            ^forcingTermParameters.alpha;
        if (isInitialIteration)
            eta_k_S = forcingTermParameters.maxForcingTerm;
        elseif (safeguardedTerm <= safeguardParameters.threshold)
            eta_k_S = min(forcingTermParameters.maxForcingTerm, eta_k_R);
        else
            eta_k_S = min(forcingTermParameters.maxForcingTerm, ...
                max(forcingTermParameters.maxForcingTerm, eta_k_R));
        end
        forcingTerm = min(forcingTermParameters.maxForcingTerm, ...
            max(eta_k_S, (0.5 * (tauA + tauR * norm(Fx_initial, Inf)) ...
            / norm(Fx, Inf))));
        residualError = forcingTerm * norm(Fx);
end
end

function flux = GenerateFluxVec(phi, indices, rows, Vx, Vy, Dxx, Dyy, ...
    rowForIndex, columnForIndex, xNodeDeltas, yNodeDeltas, northBC, eastBC, ...
    southBC, westBC, isUpwinding, isNBoundaryIndex, isEBoundaryIndex, ...
    isSBoundaryIndex, isWBoundaryIndex, time, nodesXPos, nodesYPos, ...
    nodeWidthsAll, nodeHeightsAll)

%% Initialise solution parameters

flux = zeros(length(phi), 1);

nBoundaryIndices = indices(isNBoundaryIndex(indices));
nIndices = indices(~isNBoundaryIndex(indices));

eBoundaryIndices = indices(isEBoundaryIndex(indices));
eIndices = indices(~isEBoundaryIndex(indices));

sBoundaryIndices = indices(isSBoundaryIndex(indices));
sIndices = indices(~isSBoundaryIndex(indices));

wBoundaryIndices = indices(isWBoundaryIndex(indices));
wIndices = indices(~isWBoundaryIndex(indices));

Vx_eval = Vx(phi, nodesXPos, nodesYPos, time);
Vy_eval = Vy(phi, nodesXPos, nodesYPos, time);

Dxx_eval = Dxx(phi, nodesXPos, nodesYPos, time);
Dyy_eval = Dyy(phi, nodesXPos, nodesYPos, time);

%% Formulate F(x) = 0 by considering the flux at each boundary

% North nodes
if (~isempty(nIndices))
    
    nPhi = phi(nIndices);
    northNeighbourPhi = phi(nIndices - 1);
    
    % Determine advection at face
    advectionVelocityY = Vy_eval(nIndices);
    if (isUpwinding)
        advectionAtFace = zeros(length(nIndices), 1);
        positiveAdvection = advectionVelocityY > 0;
        advectionAtFace(positiveAdvection) = phi(nIndices(positiveAdvection));
        advectionAtFace(~positiveAdvection) = ...
            phi(nIndices(~positiveAdvection) - 1);
    else
        advectionAtFace = (northNeighbourPhi + nPhi) ./ 2;
    end
    
    flux(nIndices) = flux(nIndices) + nodeWidthsAll(nIndices) ...
        .* (((Vy_eval(nIndices) + Vy_eval(nIndices - 1)) ./ 2) ...
            .* advectionAtFace ...
        - ((Dyy_eval(nIndices) + Dyy_eval(nIndices - 1)) ./ 2) ...
            .* (northNeighbourPhi - nPhi) ...
            ./ yNodeDeltas(rowForIndex(nIndices) - 1));
end
    
if (~isempty(nBoundaryIndices));    
    nBoundaryPhi = phi(nBoundaryIndices);
    diffusion = Dyy_eval(nBoundaryIndices);
    
    flux(nBoundaryIndices) = flux(nBoundaryIndices) ...
        + nodeWidthsAll(nBoundaryIndices) ...
        .* ( (Vy_eval(nBoundaryIndices) + diffusion ...
            .* northBC.A(nodesXPos(nBoundaryIndices), time) ...
            ./ northBC.B(nodesXPos(nBoundaryIndices), time)) ...
        .* nBoundaryPhi - diffusion ...
            .* northBC.C(nodesXPos(nBoundaryIndices), time) ...
            ./ northBC.B(nodesXPos(nBoundaryIndices), time) );
end

% East nodes
if (~isempty(eIndices))
    
    ePhi = phi(eIndices);
    eastNeighbourPhi = phi(eIndices + rows);
    
    % Determine advection at face
    advectionVelocityX = Vx_eval(eIndices);
    if (isUpwinding)
        advectionAtFace = zeros(length(eIndices), 1);
        positiveAdvection = advectionVelocityX > 0;
        advectionAtFace(positiveAdvection) = phi(eIndices(positiveAdvection));
        advectionAtFace(~positiveAdvection) = ...
            phi(eIndices(~positiveAdvection) + rows);
    else
        advectionAtFace = (eastNeighbourPhi + ePhi) ./ 2;
    end
    
    flux(eIndices) = flux(eIndices) + nodeHeightsAll(eIndices) ...
        .* (((Vx_eval(eIndices) + Vx_eval(eIndices + rows)) ./ 2) ...
            .* advectionAtFace ...
        - ((Dxx_eval(eIndices) + Dxx_eval(eIndices + rows)) ./ 2) ...
            .* (eastNeighbourPhi - ePhi) ...
            ./ xNodeDeltas(columnForIndex(eIndices))); 
end

if (~isempty(eBoundaryIndices))    
    eBoundaryPhi = phi(eBoundaryIndices);
    diffusion = Dxx_eval(eBoundaryIndices);
    
    flux(eBoundaryIndices) = flux(eBoundaryIndices) ...
        + nodeHeightsAll(eBoundaryIndices) ...
        .* ( (Vx_eval(eBoundaryIndices) + diffusion ...
            .* eastBC.A(nodesYPos(eBoundaryIndices), time) ...
            ./ eastBC.B(nodesYPos(eBoundaryIndices), time)) ...
        .* eBoundaryPhi - diffusion ...
            .* eastBC.C(nodesYPos(eBoundaryIndices), time) ...
            ./ eastBC.B(nodesYPos(eBoundaryIndices), time) );
end
    
% South nodes
if (~isempty(sIndices))
    
    sPhi = phi(sIndices);
    southNeighbourPhi = phi(sIndices + 1);
    
    % Determine advection at face
    advectionVelocityY = Vy_eval(sIndices);
    if (isUpwinding)
        advectionAtFace = zeros(length(sIndices), 1);
        positiveAdvection = advectionVelocityY > 0;
        advectionAtFace(~positiveAdvection) = phi(sIndices(~positiveAdvection));
        advectionAtFace(positiveAdvection) = ...
            phi(sIndices(positiveAdvection) + 1);
    else
        advectionAtFace = (southNeighbourPhi + sPhi) ./ 2;
    end
    
    flux(sIndices) = flux(sIndices) - nodeWidthsAll(sIndices) ...
        .* (((Vy_eval(sIndices) + Vy_eval(sIndices + 1)) ./ 2) ...
            .* advectionAtFace ...
        - ((Dyy_eval(sIndices) + Dyy_eval(sIndices + 1)) ./ 2) ...
            .* (sPhi - southNeighbourPhi) ...
            ./ yNodeDeltas(rowForIndex(sIndices)));
end
   
if (~isempty(sBoundaryIndices))
    sBoundaryPhi = phi(sBoundaryIndices);
    diffusion = Dyy_eval(sBoundaryIndices);

    flux(sBoundaryIndices) = flux(sBoundaryIndices) ...
        - nodeWidthsAll(sBoundaryIndices) ...
        .* ( (Vy_eval(sBoundaryIndices) - diffusion ...
            .* southBC.A(nodesXPos(sBoundaryIndices), time) ...
            ./ southBC.B(nodesXPos(sBoundaryIndices), time)) ...
        .* sBoundaryPhi + diffusion ...
            .* southBC.C(nodesXPos(sBoundaryIndices), time) ...
            ./ southBC.B(nodesXPos(sBoundaryIndices), time) );
end
    
% West nodes
if (~isempty(wIndices))
    
    wPhi = phi(wIndices);
    westNeighbourPhi = phi(wIndices - rows);
    
    % Determine advection at face
    advectionVelocityX = Vx_eval(wIndices);
    if (isUpwinding)
        advectionAtFace = zeros(length(wIndices), 1);
        positiveAdvection = advectionVelocityX > 0;
        advectionAtFace(~positiveAdvection) = phi(wIndices(~positiveAdvection));
        advectionAtFace(positiveAdvection) = ...
            phi(wIndices(positiveAdvection) - rows);
    else
        advectionAtFace = (westNeighbourPhi + wPhi) ./ 2;
    end
    
    flux(wIndices) = flux(wIndices) - nodeHeightsAll(wIndices) ...
        .* (((Vx_eval(wIndices) + Vx_eval(wIndices - rows)) ./ 2) ...
            .* advectionAtFace ...
        - ((Dxx_eval(wIndices) + Dxx_eval(wIndices - rows)) ./ 2) ...
            .* (wPhi - westNeighbourPhi) ...
            ./ xNodeDeltas(columnForIndex(wIndices) - 1));
end

if (~isempty(wBoundaryIndices))    
    wBoundaryPhi = phi(wBoundaryIndices);
    diffusion = Dxx_eval(wBoundaryIndices);
    
    flux(wBoundaryIndices) = flux(wBoundaryIndices) ...
        - nodeHeightsAll(wBoundaryIndices) ...
        .* ( (Vx_eval(wBoundaryIndices) - diffusion ...
            .* westBC.A(nodesYPos(wBoundaryIndices), time) ...
            ./ westBC.B(nodesYPos(wBoundaryIndices), time)) ...
        .* wBoundaryPhi + diffusion ...
            .* westBC.C(nodesYPos(wBoundaryIndices), time) ...
            ./ westBC.B(nodesYPos(wBoundaryIndices), time));
end

flux = flux(indices) ./ (nodeWidthsAll(indices) .* nodeHeightsAll(indices));

end

function [delta] = determine_newton_step_delta(x)
if (norm(x) == 0)
    delta = sqrt(eps);
else
    delta = sqrt(eps) * norm(x);
end
end
