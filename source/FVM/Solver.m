function [tout, yout] = Solver(tFinal, Dxx, Dyy, Vx, Vy, source, theta, ...
    advectionHandling, nodesX, nodesY, northBC, eastBC, southBC, westBC, ...
    initialCondition, storedTimeSteps)
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialise mesh parameters

nodesY = nodesY(:);
nodesX = nodesX(:);

rows = length(nodesY);
columns = length(nodesX);

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

% TODO: remove stub code once adaptive time-stepping is implemented.
tStart = 0;
dt = 0.0001;

tout = tStart:dt:tFinal;
timeSteps = length(tout);
yout = zeros(rows, columns, floor(timeSteps / storedTimeSteps) + 1);
yout(:, :, 1) = initialCondition;
previousSolution = initialCondition;

%% Iteratively solve advection-diffusion equation (time marching strategy)

for i = 1:timeSteps
    
    %% Pseudocode:
    
    % Formulate the nonlinear system F(u) = 0 for the next time step
    %   N.B. This formulation must support different temporal weightings 
    %   (see IMEX method). It must also support generalised boundary
    %   conditions.
    %
    %   Clearly only components of this nonlinear system can be formulated
    %   by this point. Namely, the components not involving phi_{n+1}. The
    %   reason for this is that the value of phi_{n+1} is determined during
    %   the next step of the algorithm (see below). As such, we can
    %   formulate a constant vector corresponding to all the components
    %   solely relying upon phi_n. This includes -phi_n and the forward 
    %   Euler component of the IMEX formulation (see notes).
    %
    
    % Formulate the Forward Euler component of F(u) = 0
    F_forwardEuler = zeros(rows, columns);
    if (theta == 0 || theta == 1/2)
        
        for j = 1:rows
            for k = 1:columns
                F_forwardEuler(j, k) = dt * (1 - theta) * GenerateFlux(j, k, ...
                    rows, columns, previousSolution, Vx, Vy, Dxx, Dyy, ...
                    xFaceDeltas, yFaceDeltas, nodeWidths, nodeHeights, ...
                    northBC, eastBC, southBC, westBC);
                F_forwardEuler(j, k) = F_forwardEuler(j, k) ...
                    - dt * (1 - theta) * source(previousSolution(j, k));
                F_forwardEuler(j, k) = F_forwardEuler(j, k) - previousSolution(j, k);
            end
        end
    end
    
    previousSolution = -1 * F_forwardEuler(:, :);
    
    % Iteratively solve F(u) = 0 for phi using an inexact Newton-GMRES solver
    %   N.B. This solver must use the formulates for the forcing term
    %   prescribed by Walker (see function description). 
    %   
    %   Preconditioning must be included. Additional handling should be 
    %   included to ensure that the preconditioner is not continually 
    %   rebuilt for each system. The preconditioner should only be 
    %   constructed for two reasons: (1) there does not exist a 
    %   preconditioner yet and (2) the current preconditioner does not 
    %   result in fast convergence within the GMRES-arnoldi iterations 
    %   (be careful of small restart values).
    %
    %   Do we need to "construct" F(u) for each iteration of Newton's
    %   method? Ans: between iterations phi_n does not change, however
    %   phi_{n+1} does. This means that any component involving phi_{n+1}
    %   need to be updated. This is equivalent to evaluating F(u) for a
    %   given value of u. More specifically, we need to construct all the
    %   components solely relying upon phi_{n+1} in this step. This
    %   includes phi_{n+1} itself and the backward Euler component of the
    %   IMEX method (see notes). We then add this result onto the existing
    %   vector of values based directly upon phi_{n} (see previous step).
    %
    %   Error handling should be included in the event that a sufficiently
    %   accurate solution is not found. This should terminate the solver
    %   entirely.
    %   
    
    % Optionally store the solution
    if (mod(i, storedTimeSteps) == 0)
        lastStoredSolutionIndex = i / storedTimeSteps + 1;
        yout(:, :, lastStoredSolutionIndex) = previousSolution;
    end
end

end

%
%   Helper Functions
%

function flux = GenerateFlux(row, column, rows, columns, ...
    previousSolution, Vx, Vy, Dxx, Dyy, xFaceDeltas, yFaceDeltas, ...
    nodeWidths, nodeHeights, northBC, eastBC, southBC, westBC)
%% GenerateFlux: ...
%

%% Intiailise known constants and solution variables

MIN_INDEX = 1;

cv_prevSolution = previousSolution(row, column); 
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
            'for a boundary condition cannot be negative.']);
    end
else
    northPrevSolution = previousSolution(row - 1, column);
    flux = flux + cv_width ...
        * ( cv_Vy * (northPrevSolution + cv_prevSolution) / 2 ...
        - cv_Dyy * (northPrevSolution - cv_prevSolution) / yFaceDeltas(row - 1) );
end

% East face
if (column == columns)
    if (eastBC.B ~= 0)
        flux = flux + cv_height ...
            * ( (cv_Vx + cv_Dxx * eastBC.A / eastBC.B) * cv_prevSolution ...
            - cv_Dxx * eastBC.C / eastBC.B );
    else
        error(['Invalid Boundary Condition (East): The B coefficient ' ...
            'for a boundary condition cannot be negative.']);
    end
else
    eastPrevSolution = previousSolution(row, column + 1);
    flux = flux + cv_height ...
        * ( cv_Vx * (eastPrevSolution + cv_prevSolution) / 2 ...
        - cv_Dxx * (eastPrevSolution - cv_prevSolution) / xFaceDeltas(column) );
end

% South face
if (row == rows)
    if (southBC.B ~= 0)
        flux = flux - cv_width ...
            * ( (cv_Vy - cv_Dyy * southBC.A / southBC.B) * cv_prevSolution ...
            + cv_Dyy * southBC.C / southBC.B );
    else
        error(['Invalid Boundary Condition (South): The B coefficient ' ...
            'for a boundary condition cannot be negative.']);
    end
else
    southPrevSolution = previousSolution(row + 1, column);
    flux = flux - cv_width ...
        * ( cv_Vy * (southPrevSolution + cv_prevSolution) / 2 ...
        - cv_Dyy * (cv_prevSolution - southPrevSolution) / yFaceDeltas(row) );
end

% West face
if (column == MIN_INDEX)
    if (westBC.B ~= 0)
        flux = flux - cv_height ...
            * ( (cv_Vx - cv_Dxx * westBC.A / westBC.B) * cv_prevSolution ...
            + cv_Dxx * westBC.C / westBC.B );
    else
        error(['Invalid Boundary Condition (West): The B coefficient ' ...
            'for a boundary condition cannot be negative.']);
    end
else
    westPrevSolution = previousSolution(row, column - 1);
    flux = flux - cv_height ...
        * ( cv_Vx * (westPrevSolution + cv_prevSolution) / 2 ...
        - cv_Dxx * (cv_prevSolution - westPrevSolution) / xFaceDeltas(column - 1) );
end

flux = flux / (cv_width * cv_height);

end
