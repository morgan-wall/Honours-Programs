function [tout, yout] = Solver(tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps, newtonParameters, gmresParameters)
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
%           relErrorTol:
%               The relative error tolerance applied to the inexact Newton
%               method.
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
%   1. Adaptive time-stepping.
%   2. Inexact Newton method (incl. forcing terms proposed by Walker).
%   3. Selective construction of the preconditioner in the GMRES method.
%   4. Jacobian-free Newton-GMRES solver.
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

%% Initialise solution parameters

tStart = 0;
dt = 0.001;

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

%% Iteratively solve advection-diffusion equation (time marching strategy)

for i = 1:timeSteps
    
    % Formulate the Forward Euler component of F(u) = 0
    F_forwardEuler = zeros(nodeCount, 1);
    for j = 1:nodeCount
        if (theta ~= 1)
            F_forwardEuler(j) = dt * (1 - theta) * GenerateFlux(j, ...
                    rows, columns, previousSolution, Vx, Vy, Dxx, Dyy, ...
                    xNodeDeltas, yNodeDeltas, nodeWidths, nodeHeights, ...
                    northBC, eastBC, southBC, westBC, advectionHandling);
            F_forwardEuler(j) = F_forwardEuler(j) ...
                - dt * (1 - theta) * source(previousSolution(j));
        end
        
        F_forwardEuler(j) = F_forwardEuler(j) - previousSolution(j);
    end
    
    % Initialise variables for Newton-GMRES solver
    current_iteration = 0;
    
    % Formulate the Backward Euler component of F(u) = 0 
    if (i == FIRST_TIME_STEP_INDEX)
        F_current_backwardEuler = zeros(nodeCount, 1);
        
        for j = 1:nodeCount
            F_current_backwardEuler(j) = dt * theta * GenerateFlux(j, ...
                rows, columns, previousSolution, Vx, Vy, Dxx, Dyy, ...
                xNodeDeltas, yNodeDeltas, nodeWidths, nodeHeights, ...
                northBC, eastBC, southBC, westBC, advectionHandling);
            F_current_backwardEuler(j) = F_current_backwardEuler(j) ...
                - dt * theta * source(previousSolution(j));
            F_current_backwardEuler(j) = ...
                F_current_backwardEuler(j) + previousSolution(j);
        end
    end
    
    Fx = F_forwardEuler + F_current_backwardEuler;
    
    % Solve the non-linear system, F(x) = 0, for the next time step
    currentSolution = previousSolution;
    while (current_iteration <= newtonParameters.maxIterations ...
            && norm(Fx) > newtonParameters.relErrorTol)
        
        % determine finite difference approx of Jacobian
        jacobian = zeros(nodeCount);
        
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
        
        % solve the linear system using GMRES
        delta_x = gmres_general(jacobian, Fx, currentSolution, ...
            gmresParameters.maxIterations, gmresParameters.restartValue, ...
            gmresParameters.errorTol, gmresParameters.preconditioningType, ...
            gmresParameters.omega);

        currentSolution = currentSolution - delta_x;
        
        % Evaluate the nonlinear system for updated iterate
        F_current_backwardEuler = zeros(nodeCount, 1);
        for k = 1:nodeCount
            F_current_backwardEuler(k) = dt * theta * GenerateFlux(k, ...
                    rows, columns, currentSolution, Vx, Vy, Dxx, Dyy, ...
                    xNodeDeltas, yNodeDeltas, nodeWidths, nodeHeights, ...
                    northBC, eastBC, southBC, westBC, advectionHandling);
            F_current_backwardEuler(k) = F_current_backwardEuler(k) ...
                - dt * theta * source(currentSolution(k));
            F_current_backwardEuler(k) = ...
                F_current_backwardEuler(k) + currentSolution(k);
        end

        Fx = F_current_backwardEuler + F_forwardEuler;
        
        current_iteration = current_iteration + 1;
    end
    
    % Ensure an accurate solution to F(u) = 0 was determined
    if (current_iteration > newtonParameters.maxIterations ...
            && norm(Fx) > newtonParameters.relErrorTol)
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

%% Determine flux for control volume

% North face
if (row == MIN_INDEX)
    if (northBC.B ~= 0)
        flux = flux + cv_width ...
            * ( (cv_Vy + cv_Dyy * northBC.A / northBC.B) * cv_prevSolution ...
            - cv_Dyy * northBC.C / northBC.B );
    else
        error(['Invalid Boundary Condition (North): The B coefficient ' ...
            'for a boundary condition must be positive.']);
    end
else
    northPrevSolution = previousSolution(nodeCount - 1);
    
    advectionAtFace = 0;
    if (strcmp(advectionHandling, 'upwinding'))
        if (cv_Vx > 0)
            advectionAtFace = cv_prevSolution;
        else
            advectionAtFace = northPrevSolution;
        end
    elseif (strcmp(advectionHandling, 'averaging'))
        advectionAtFace = (northPrevSolution + cv_prevSolution) / 2;
    end
    
    flux = flux + cv_width ...
        * ( cv_Vy * advectionAtFace ...
        - cv_Dyy * (northPrevSolution - cv_prevSolution) / yNodeDeltas(row - 1) );
end

% East face
if (column == columns)
    if (eastBC.B ~= 0)
        flux = flux + cv_height ...
            * ( (cv_Vx + cv_Dxx * eastBC.A / eastBC.B) * cv_prevSolution ...
            - cv_Dxx * eastBC.C / eastBC.B );
    else
        error(['Invalid Boundary Condition (East): The B coefficient ' ...
            'for a boundary condition must be positive.']);
    end
else
    eastPrevSolution = previousSolution(nodeCount + rows);
    
    advectionAtFace = 0;
    if (strcmp(advectionHandling, 'upwinding'))
        if (cv_Vx > 0)
            advectionAtFace = cv_prevSolution;
        else
            advectionAtFace = eastPrevSolution;
        end
    elseif (strcmp(advectionHandling, 'averaging'))
        advectionAtFace = (eastPrevSolution + cv_prevSolution) / 2;
    end
    
    flux = flux + cv_height ...
        * ( cv_Vx * advectionAtFace ...
        - cv_Dxx * (eastPrevSolution - cv_prevSolution) / xNodeDeltas(column) );
end

% South face
if (row == rows)
    if (southBC.B ~= 0)
        flux = flux - cv_width ...
            * ( (cv_Vy - cv_Dyy * southBC.A / southBC.B) * cv_prevSolution ...
            + cv_Dyy * southBC.C / southBC.B );
    else
        error(['Invalid Boundary Condition (South): The B coefficient ' ...
            'for a boundary condition must be positive.']);
    end
else
    southPrevSolution = previousSolution(nodeCount + 1);
    
    advectionAtFace = 0;
    if (strcmp(advectionHandling, 'upwinding'))
        if (cv_Vx > 0)
            advectionAtFace = southPrevSolution;
        else
            advectionAtFace = cv_prevSolution;
        end
    elseif (strcmp(advectionHandling, 'averaging'))
        advectionAtFace = (southPrevSolution + cv_prevSolution) / 2;
    end
    
    flux = flux - cv_width ...
        * ( cv_Vy * advectionAtFace ...
        - cv_Dyy * (cv_prevSolution - southPrevSolution) / yNodeDeltas(row) );
end

% West face
if (column == MIN_INDEX)
    if (westBC.B ~= 0)
        flux = flux - cv_height ...
            * ( (cv_Vx - cv_Dxx * westBC.A / westBC.B) * cv_prevSolution ...
            + cv_Dxx * westBC.C / westBC.B );
    else
        error(['Invalid Boundary Condition (West): The B coefficient ' ...
            'for a boundary condition must be positive.']);
    end
else
    westPrevSolution = previousSolution(nodeCount - rows);
    
    advectionAtFace = 0;
    if (strcmp(advectionHandling, 'upwinding'))
        if (cv_Vx > 0)
            advectionAtFace = westPrevSolution;
        else
            advectionAtFace = cv_prevSolution;
        end
    elseif (strcmp(advectionHandling, 'averaging'))
        advectionAtFace = (westPrevSolution + cv_prevSolution) / 2;
    end
    
    flux = flux - cv_height ...
        * ( cv_Vx * advectionAtFace ...
        - cv_Dxx * (cv_prevSolution - westPrevSolution) / xNodeDeltas(column - 1) );
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
