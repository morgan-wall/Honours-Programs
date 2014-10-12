function [tout, yout] = Solver(dt, tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters, ...
    forcingTermParameters, safeguardParameters)
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
%   There are no errors thrown by this function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unimplemented features:
%
%   1. Selective construction of the preconditioner in the GMRES method.
%   2. Jacobian-free Newton-GMRES solver.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Validate input parameters

if (northBC.B == 0)
    error(['Invalid Boundary Condition (North): The B coefficient ' ...
            'for a boundary condition must be positive.']);
end

if (eastBC.B == 0)
    error(['Invalid Boundary Condition (East): The B coefficient ' ...
            'for a boundary condition must be positive.']);
end

if (southBC.B == 0)
    error(['Invalid Boundary Condition (South): The B coefficient ' ...
            'for a boundary condition must be positive.']);
end

if (westBC.B == 0)
    error(['Invalid Boundary Condition (West): The B coefficient ' ...
            'for a boundary condition must be positive.']);
end

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

% Store indices that lie on each boundary
isNBoundaryIndex = mod(indices - 1, rows) == 0;
isEBoundaryIndex = indices > nodeCount - rows;
isSBoundaryIndex = mod(indices, rows) == 0;
isWBoundaryIndex = indices <= rows;

%% Initialise solution parameters

tStart = 0;
timeSteps = length(tStart:dt:tFinal) - 1;

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

newtonIterations = 0;
gmresCalls = 0;

%% Initialise solver parameters

isUpwinding = strcmp(advectionHandling, 'upwinding');

tauA = 1e-6;
tauR = 1e-6;

%% Iteratively solve advection-diffusion equation (time marching strategy)

for i = 1:timeSteps
    
    % Formulate the Forward Euler component of F(u) = 0
    F_forwardEuler = zeros(nodeCount, 1);
    if (theta ~= 1)
        F_forwardEuler = dt .* (1 - theta) ...
            .* GenerateFluxVec(previousSolution, indices, rows, Vx, Vy, ...
            Dxx, Dyy, nodeWidths, nodeHeights, rowForIndex, columnForIndex, ...
            xNodeDeltas, yNodeDeltas, northBC, eastBC, southBC, westBC, ...
            isUpwinding, isNBoundaryIndex, isEBoundaryIndex, ...
            isSBoundaryIndex, isWBoundaryIndex);
    end
    F_forwardEuler = F_forwardEuler ...
        - dt .* (1 - theta) .* source(previousSolution);
    F_forwardEuler = F_forwardEuler - previousSolution;

    % Initialise variables for Newton-GMRES solver
    current_iteration = 0;
    
    % Formulate the Backward Euler component of F(u) = 0 
    if (i == FIRST_TIME_STEP_INDEX)
        F_current_backwardEuler = dt .* theta ...
            .* GenerateFluxVec(previousSolution, indices, rows, Vx, Vy, ...
            Dxx, Dyy, nodeWidths, nodeHeights, rowForIndex, columnForIndex, ...
            xNodeDeltas, yNodeDeltas, northBC, eastBC, southBC, westBC, ...
            isUpwinding, isNBoundaryIndex, isEBoundaryIndex, ...
            isSBoundaryIndex, isWBoundaryIndex);

        F_current_backwardEuler = F_current_backwardEuler ...
            - dt .* theta .* source(previousSolution);
        F_current_backwardEuler = F_current_backwardEuler + previousSolution;
    end
    
    Fx = F_forwardEuler + F_current_backwardEuler;
    Fx_previous = Fx;
    Fx_initial_norm = norm(Fx);
    Fx_initial = Fx;
    
    % Solve the non-linear system, F(x) = 0, for the next time step
    jacobian = zeros(nodeCount);
    currentSolution = previousSolution;
    forcingTerm = 1/2;
    delta_x = realmax;
    while (current_iteration <= newtonParameters.maxIterations ...
            && norm(delta_x / currentSolution) >= newtonParameters.tolUpdate ...
            && norm(Fx) >= newtonParameters.tolResidual * Fx_initial_norm)
        
        % determine finite difference approx of Jacobian
        previousJacobian = jacobian;
        
        h = determine_newton_step_delta(currentSolution);
        
        % self-dependent Jacobian components
        for j = 1:nodeCount
            nodeIndex = j;
            
            xStepped = currentSolution;
            xStepped(nodeIndex) = xStepped(nodeIndex) + h;
            
            F_backwardEuler_stepped = dt * theta * GenerateFlux(j, ...
                    rows, columns, xStepped, Vx, Vy, Dxx, Dyy, ...
                    xNodeDeltas, yNodeDeltas, nodeWidths, nodeHeights, ...
                    northBC, eastBC, southBC, westBC, advectionHandling);
                
            F_backwardEuler_stepped = F_backwardEuler_stepped ...
                - dt * theta * source(xStepped(j));
            F_backwardEuler_stepped = F_backwardEuler_stepped + xStepped(j);
            
            F_stepped = F_backwardEuler_stepped + F_forwardEuler(j);
            
            jacobian(j, nodeIndex) = (F_stepped - Fx(j)) / h;
        end
        
        % north-dependent Jacobian components
        for j = 1:nodeCount
            nodeIndex = j - 1;
            
            if (nodeIndex >= 1 && nodeIndex <= nodeCount)
                xStepped = currentSolution;
                xStepped(nodeIndex) = xStepped(nodeIndex) + h;
                
                F_backwardEuler_stepped = dt * theta * GenerateFlux(j, ...
                        rows, columns, xStepped, Vx, Vy, Dxx, Dyy, ...
                        xNodeDeltas, yNodeDeltas, nodeWidths, nodeHeights, ...
                        northBC, eastBC, southBC, westBC, advectionHandling);

                F_backwardEuler_stepped = F_backwardEuler_stepped ...
                    - dt * theta * source(xStepped(j));
                F_backwardEuler_stepped = F_backwardEuler_stepped + xStepped(j);

                F_stepped = F_backwardEuler_stepped + F_forwardEuler(j);
                
                jacobian(j, nodeIndex) = (F_stepped - Fx(j)) / h;
            end
        end
        
        % south-dependent Jacobian components
        for j = 1:nodeCount
            nodeIndex = j + 1;
            
            if (nodeIndex >= 1 && nodeIndex <= nodeCount)
                xStepped = currentSolution;
                xStepped(nodeIndex) = xStepped(nodeIndex) + h;
                
                F_backwardEuler_stepped = dt * theta * GenerateFlux(j, ...
                        rows, columns, xStepped, Vx, Vy, Dxx, Dyy, ...
                        xNodeDeltas, yNodeDeltas, nodeWidths, nodeHeights, ...
                        northBC, eastBC, southBC, westBC, advectionHandling);

                F_backwardEuler_stepped = F_backwardEuler_stepped ...
                    - dt * theta * source(xStepped(j));
                F_backwardEuler_stepped = F_backwardEuler_stepped + xStepped(j);

                F_stepped = F_backwardEuler_stepped + F_forwardEuler(j);
                
                jacobian(j, nodeIndex) = (F_stepped - Fx(j)) / h;
            end
        end
        
        % east-dependent Jacobian components
        for j = 1:nodeCount
            nodeIndex = j + rows;
            
            if (nodeIndex >= 1 && nodeIndex <= nodeCount)
                xStepped = currentSolution;
                xStepped(nodeIndex) = xStepped(nodeIndex) + h;
                
                F_backwardEuler_stepped = dt * theta * GenerateFlux(j, ...
                        rows, columns, xStepped, Vx, Vy, Dxx, Dyy, ...
                        xNodeDeltas, yNodeDeltas, nodeWidths, nodeHeights, ...
                        northBC, eastBC, southBC, westBC, advectionHandling);

                F_backwardEuler_stepped = F_backwardEuler_stepped ...
                    - dt * theta * source(xStepped(j));
                F_backwardEuler_stepped = F_backwardEuler_stepped + xStepped(j);

                F_stepped = F_backwardEuler_stepped + F_forwardEuler(j);
                
                jacobian(j, nodeIndex) = (F_stepped - Fx(j)) / h;
            end
        end
        
        % west-dependent Jacobian components
        for j = 1:nodeCount
            nodeIndex = j - rows;
            
            if (nodeIndex >= 1 && nodeIndex <= nodeCount)
                xStepped = currentSolution;
                xStepped(nodeIndex) = xStepped(nodeIndex) + h;
                
                F_backwardEuler_stepped = dt * theta * GenerateFlux(j, ...
                        rows, columns, xStepped, Vx, Vy, Dxx, Dyy, ...
                        xNodeDeltas, yNodeDeltas, nodeWidths, nodeHeights, ...
                        northBC, eastBC, southBC, westBC, advectionHandling);
                    
                F_backwardEuler_stepped = F_backwardEuler_stepped ...
                    - dt * theta * source(xStepped(j));
                F_backwardEuler_stepped = F_backwardEuler_stepped + xStepped(j);
                
                F_stepped = F_backwardEuler_stepped + F_forwardEuler(j);
                
                jacobian(j, nodeIndex) = (F_stepped - Fx(j)) / h;
            end
        end
        
        if (current_iteration == 0)
            previousJacobian = jacobian;
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
            gmresParameters.maxIterations, gmresParameters.restartValue, ...
            residualError, gmresParameters.preconditioningType, ...
            gmresParameters.omega);
        
        gmresCalls = gmresCalls + iterations;
        
        currentSolution = currentSolution - delta_x;
        
        % Evaluate the nonlinear system for updated iterate
        F_current_backwardEuler = dt .* theta ...
            .* GenerateFluxVec(currentSolution, indices, rows, Vx, Vy, ...
            Dxx, Dyy, nodeWidths, nodeHeights, rowForIndex, columnForIndex, ...
            xNodeDeltas, yNodeDeltas, northBC, eastBC, southBC, westBC, ...
            isUpwinding, isNBoundaryIndex, isEBoundaryIndex, ...
            isSBoundaryIndex, isWBoundaryIndex);

        F_current_backwardEuler = F_current_backwardEuler ...
            - dt .* theta .* source(currentSolution);
        F_current_backwardEuler = F_current_backwardEuler + currentSolution;

        Fx_previous = Fx;
        Fx = F_current_backwardEuler + F_forwardEuler;
        
        current_iteration = current_iteration + 1;
    end
    
    newtonIterations = newtonIterations + current_iteration;
    
    % Ensure an accurate solution to F(u) = 0 was determined
    if (current_iteration > newtonParameters.maxIterations ...
            && norm(delta_x / currentSolution) < newtonParameters.tolUpdate ...
            && norm(Fx) < newtonParameters.tolResidual * Fx_initial_norm)
        error(['Method Failure: the non-linear system generated at '...
            't = ' num2str(i * dt) ' was not solved using ' ...
            num2str(max_iterations) ' iterations of inexact Newton ' ...
            'method.']);
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

end

%
%   Helper Functions
%

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
    nodeWidths, nodeHeights, rowForIndex, columnForIndex, ...
    xNodeDeltas, yNodeDeltas, northBC, eastBC, southBC, westBC, ...
    isUpwinding, isNBoundaryIndex, isEBoundaryIndex, isSBoundaryIndex, ...
    isWBoundaryIndex)

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

%% Formulate F(x) = 0 by considering the flux at each boundary

% North nodes
if (~isempty(nIndices))
    
    nPhi = phi(nIndices);
    
    % Determine advection at face
    advectionAtFace = zeros(length(flux));
    advectionVelocityY = Vy(phi(nIndices));
    if (isUpwinding)
        positiveAdvection = advectionVelocityY > 0;
        advectionAtFace(positiveAdvection) = phi(nIndices(positiveAdvection));
        advectionAtFace(~positiveAdvection) = ...
            phi(nIndices(~positiveAdvection) - 1);
    else
        advectionAtFace = (phi(nIndices - 1) + nPhi) ./ 2;
    end
    
    flux(nIndices) = flux(nIndices) + nodeWidths(columnForIndex(nIndices)) ...
        .* (Vy(nPhi) .* advectionAtFace - Dyy(nPhi) ...
        .* (phi(nIndices - 1) - nPhi) ...
        ./ yNodeDeltas(rowForIndex(nIndices) - 1));
end
    
if (~isempty(nBoundaryIndices));
    nBoundaryPhi = phi(nBoundaryIndices);
    diffusion = Dyy(phi(nBoundaryIndices));
    
    flux(nBoundaryIndices) = flux(nBoundaryIndices) ...
        + nodeWidths(columnForIndex(nBoundaryIndices)) ...
        .* ( (Vy(nBoundaryPhi) + diffusion .* northBC.A ./ northBC.B) ...
        .* nBoundaryPhi - diffusion .* northBC.C ./ northBC.B );
end

% East nodes
if (~isempty(eIndices))
    
    ePhi = phi(eIndices);
    
    % Determine advection at face
    advectionAtFace = zeros(length(flux));
    advectionVelocityX = Vx(phi(eIndices));
    if (isUpwinding)
        positiveAdvection = advectionVelocityX > 0;
        advectionAtFace(positiveAdvection) = phi(eIndices(positiveAdvection));
        advectionAtFace(~positiveAdvection) = ...
            phi(eIndices(~positiveAdvection) + rows);
    else
        advectionAtFace = (phi(eIndices + rows) + ePhi) ./ 2;
    end
    
    flux(eIndices) = flux(eIndices) + nodeHeights(rowForIndex(eIndices)) ...
        .* (Vx(ePhi) .* advectionAtFace - Dxx(ePhi) ...
        .* (phi(eIndices + rows) - ePhi) ...
        ./ xNodeDeltas(columnForIndex(eIndices))); 
end

if (~isempty(eBoundaryIndices))
    eBoundaryPhi = phi(eBoundaryIndices);
    diffusion = Dxx(phi(eBoundaryIndices));
    
    flux(eBoundaryIndices) = flux(eBoundaryIndices) ...
        + nodeHeights(rowForIndex(eBoundaryIndices)) ...
        .* ( (Vx(eBoundaryPhi) + diffusion .* eastBC.A ./ eastBC.B) ...
        .* eBoundaryPhi - diffusion .* eastBC.C ./ eastBC.B );
end
    
% South nodes
if (~isempty(sIndices))
    
    sPhi = phi(sIndices);
    
    % Determine advection at face
    advectionAtFace = zeros(length(flux));
    advectionVelocityY = Vy(phi(sIndices));
    if (isUpwinding)
        positiveAdvection = advectionVelocityY > 0;
        advectionAtFace(~positiveAdvection) = phi(sIndices(~positiveAdvection));
        advectionAtFace(positiveAdvection) = ...
            phi(sIndices(positiveAdvection) + 1);
    else
        advectionAtFace = (phi(sIndices + 1) + sPhi) ./ 2;
    end
    
    flux(sIndices) = flux(sIndices) - nodeWidths(columnForIndex(sIndices)) ...
        .* (Vy(sPhi) .* advectionAtFace - Dyy(sPhi) ...
        .* (sPhi - phi(sIndices + 1)) ...
        ./ yNodeDeltas(rowForIndex(sIndices)));
end
   
if (~isempty(sBoundaryIndices))
    sBoundaryPhi = phi(sBoundaryIndices);
    diffusion = Dyy(phi(sBoundaryIndices));
    
    flux(sBoundaryIndices) = flux(sBoundaryIndices) ...
        - nodeWidths(columnForIndex(sBoundaryIndices)) ...
        .* ( (Vy(sBoundaryPhi) - diffusion .* southBC.A ./ southBC.B) ...
        .* sBoundaryPhi - diffusion .* southBC.C ./ southBC.B );
end
    
% West nodes
if (~isempty(wIndices))
    
    wPhi = phi(wIndices);
    
    % Determine advection at face
    advectionAtFace = zeros(length(flux));
    advectionVelocityX = Vx(phi(wIndices));
    if (isUpwinding)
        positiveAdvection = advectionVelocityX > 0;
        advectionAtFace(~positiveAdvection) = phi(wIndices(~positiveAdvection));
        advectionAtFace(positiveAdvection) = ...
            phi(wIndices(positiveAdvection) - rows);
    else
        advectionAtFace = (phi(wIndices - rows) + wPhi) ./ 2;
    end
    
    flux(wIndices) = flux(wIndices) - nodeHeights(rowForIndex(wIndices)) ...
        .* (Vx(wPhi) .* advectionAtFace - Dxx(wPhi) ...
        .* (wPhi - phi(wIndices - rows)) ...
        ./ xNodeDeltas(columnForIndex(wIndices) - 1));
end

if (~isempty(wBoundaryIndices))
    wBoundaryPhi = phi(wBoundaryIndices);
    diffusion = Dxx(phi(wBoundaryIndices));
    
    flux(wBoundaryIndices) = flux(wBoundaryIndices) ...
        - nodeHeights(rowForIndex(wBoundaryIndices)) ...
        .* ( (Vx(wBoundaryPhi) - diffusion .* westBC.A ./ westBC.B) ...
        .* wBoundaryPhi + diffusion .* westBC.C ./ westBC.B);
end

flux = flux(indices);
flux = flux ./ ...
    (nodeWidths(columnForIndex(indices)) .* nodeHeights(rowForIndex(indices)));

end

function flux = GenerateFlux(nodeCount, rows, columns, ...
    previousSolution, Vx, Vy, Dxx, Dyy, xNodeDeltas, yNodeDeltas, ...
    nodeWidths, nodeHeights, northBC, eastBC, southBC, westBC, ...
    advectionHandling)
%% GenerateFlux: ...
%

%% Intiailise known constants and solution variables

MIN_INDEX = 1;

row = mod(nodeCount - 1, rows) + 1;
column = floor((nodeCount - 1) / rows) + 1;

cv_prevSolution = previousSolution(nodeCount); 
cv_Dxx = Dxx(cv_prevSolution);
cv_Dyy = Dyy(cv_prevSolution);
cv_Vx = Vx(cv_prevSolution);
cv_Vy = Vy(cv_prevSolution);
cv_width = nodeWidths(column);
cv_height = nodeHeights(row);

flux = 0;

isUpwinding = strcmp(advectionHandling, 'upwinding');

%% Determine flux for control volume

% North face
if (row == MIN_INDEX)
    flux = flux + cv_width ...
        * ( (cv_Vy + cv_Dyy * northBC.A / northBC.B) * cv_prevSolution ...
        - cv_Dyy * northBC.C / northBC.B );
else
    northPrevSolution = previousSolution(nodeCount - 1);
    
    if (isUpwinding)
        if (cv_Vy > 0)
            advectionAtFace = cv_prevSolution;
        else
            advectionAtFace = northPrevSolution;
        end
    else
        advectionAtFace = (northPrevSolution + cv_prevSolution) / 2;
    end
    
    flux = flux + cv_width ...
        * (cv_Vy * advectionAtFace ...
        - cv_Dyy * (northPrevSolution - cv_prevSolution) / yNodeDeltas(row - 1));
end

% East face
if (column == columns)
    flux = flux + cv_height ...
        * ((cv_Vx + cv_Dxx * eastBC.A / eastBC.B) * cv_prevSolution ...
        - cv_Dxx * eastBC.C / eastBC.B);
else
    eastPrevSolution = previousSolution(nodeCount + rows);
    
    if (isUpwinding)
        if (cv_Vx > 0)
            advectionAtFace = cv_prevSolution;
        else
            advectionAtFace = eastPrevSolution;
        end
    else
        advectionAtFace = (eastPrevSolution + cv_prevSolution) / 2;
    end
    
    flux = flux + cv_height ...
        * (cv_Vx * advectionAtFace ...
        - cv_Dxx * (eastPrevSolution - cv_prevSolution) / xNodeDeltas(column));
end

% South face
if (row == rows)
    flux = flux - cv_width ...
        * ((cv_Vy - cv_Dyy * southBC.A / southBC.B) * cv_prevSolution ...
        + cv_Dyy * southBC.C / southBC.B);
else
    southPrevSolution = previousSolution(nodeCount + 1);
    
    if (isUpwinding)
        if (cv_Vy > 0)
            advectionAtFace = southPrevSolution;
        else
            advectionAtFace = cv_prevSolution;
        end
    else
        advectionAtFace = (southPrevSolution + cv_prevSolution) / 2;
    end
    
    flux = flux - cv_width ...
        * (cv_Vy * advectionAtFace ...
        - cv_Dyy * (cv_prevSolution - southPrevSolution) / yNodeDeltas(row));
end

% West face
if (column == MIN_INDEX)
    flux = flux - cv_height ...
        * ( (cv_Vx - cv_Dxx * westBC.A / westBC.B) * cv_prevSolution ...
        + cv_Dxx * westBC.C / westBC.B );
else
    westPrevSolution = previousSolution(nodeCount - rows);
    
    if (isUpwinding)
        if (cv_Vx > 0)
            advectionAtFace = westPrevSolution;
        else
            advectionAtFace = cv_prevSolution;
        end
    else
        advectionAtFace = (westPrevSolution + cv_prevSolution) / 2;
    end
    
    flux = flux - cv_height ...
        * (cv_Vx * advectionAtFace ...
        - cv_Dxx * (cv_prevSolution - westPrevSolution) / xNodeDeltas(column - 1));
end

flux = flux / (cv_width * cv_height);

end

function [delta] = determine_newton_step_delta(x)
if (norm(x) == 0)
    delta = sqrt(eps);
else
    delta = sqrt(eps) * norm(x);
end
end
