  function [sol] = ddef2(F, phi, tau, j, left, right, hIn)
% WHAT'S NEW IN THIS VERSION 
% * SWITCHED MESH GETTING FUNCTION TO GETMESH_V4

%   A program to solve the initial value problem y'(t) = F(y(t), D(t)) where
%   D is a distributed delay. F must be entered as an anonymous function F =
%   @(x, D). phi is the history function for the IVP, also an anonymous
%   function. left and right are the endpoints of the solution interval for
%   the IVP, and hIn is the user requested Runge-Kutta method mesh size. 
%   numSteps is the number of subintervals the solver will use
%   to find the solution y(t) on [left, right].
    
%   Defining some constants to be used throughout the script
%--------------------------------------------------------------------------
%    INPUT VALIDATION

if ~isa(F,'function_handle')
    ME = MException('InputException:TypeError',...
        'F is not defined as a function handle.');
    throw(ME)
end
if ~isa(phi,'function_handle')
    ME = MException('InputException:TypeError',...
        'phi is not defined as a function handle.');
    throw(ME)
end
if j <= 0
    ME = MException('InputException:OutOfBounds', ...
        'Variable j is <=0');
    throw(ME)
end
if tau <= 0
    ME = MException('InputException:OutOfBounds', ...
        'Variable tau is <=0');
    throw(ME)
end
if right <= left
    ME = MException('InputException:OutOfBounds', ...
        'Interval [right, left] not defined or degenerate. Ensure right > left.');
    throw(ME)
end
if hIn <= 0
    ME = MException('InputException:OutOfBounds', ...
        'Enter a step size h > 0');
    throw(ME)
end
if hIn > right - left
    ME = MException('InputException:OutOfBounds', ...
        'Stepsize h > right - left; stepsize larger than interval');
    throw(ME)
end
%--------------------------------------------------------------------------
    numSteps = floor((right - left)./hIn) + 1; %determine uniform step size for the IVP based on the user input step size.

    a = j/tau;
    h = (right - left) /numSteps; %define the step size for IVP solution
    k = 1 + 1/j;
    beta = 1 + 5/j;
    A = (j+1)/a^(1/beta); % to be passed to the simpson's rule function
    g= @(t)((beta.*a.^j.* (k.*j./a.^(1./beta)).^(beta.*j))./gamma(j)).* exp(-a.*(-k.*j.*log(t)./a.^(1./beta)).^beta).*(-log(t)).^(beta .*j-1).*(1./t);
    
%-------------------------------------------------------------------------

%   global variables pertaining to numerical quadrature
    intCoeffs = [2, -1, 2, 2, -1, 2, 2, -1, 2]; %this will eventually be a vector of [2, -1, 2] growing as necessary
    K = [0, 0]; %this is a global variable containing the stage values K_i
    hQuadrature = 0.1 .* sqrt(h); %this is the step size which will be used for numerical integration 
    
%-----------------------------------------------------------------------

    %create space to store the solution
    SOL = zeros(numSteps+1, 4); %preallocate solution for speed. 
    SOL(1, 1) = left;
    SOL(1, 2) = phi(left); %set the value of history function at t_0
    y_n = SOL(1, 2);
    
%-----------------------------------------------------------------------
    
%   find K(1) for the first time
    t_n = left;
    quadMesh = getMesh_v4(hQuadrature, 1);
    transformedMesh = t_n - (-A.*log(quadMesh.mesh)).^(beta);
    intMeshSize = size(transformedMesh, 1);
    gMesh = g(quadMesh.mesh); 
    gMeshH = multiplyByH(gMesh, quadMesh.h, quadMesh.hTemp1, quadMesh.hTemp2, quadMesh.index, intMeshSize);
    [evalFirstPart, reuseIndex] = assignFirstPart(SOL, transformedMesh, intMeshSize, left, t_n, phi);
    I_num = integrate(1, gMeshH, evalFirstPart);
    K(1) = F(y_n, I_num);

%   main loop
    for i = 1 : numSteps
        %disp(num2str(t_n));
%       set t_right
        t_right = t_n + h;
        
%       get a new mesh
        specialPoint = exp((-1/A).*(t_right - left).^(1/beta));
        quadMesh = getMesh_v4(hQuadrature, specialPoint); %generate the mesh which will be used for the integration up to t_right
        gMesh = g(quadMesh.mesh); 
        
%       transform the mesh to be over (-\infty, t_right)
        transformedMesh = t_right - (-A.*log(quadMesh.mesh)).^(beta);

%       Assign function values to the mesh points which fall in (-\infty, t_n)
%       first call 'assignFirstPart' which will assign all function
%       values on (-\infty, t_n)
        intMeshSize = size(transformedMesh, 1);
        [evalFirstPart, reuseIndex] = assignFirstPart(SOL, transformedMesh, intMeshSize, left, t_n, phi, h); %reuseIndex is the index of the first point which must be evaluated at time > t_n
%       once reuseIndex is known, we can pre-multiply each mesh point by 
%       its corresponding hQuadrature (will only be different for the
%       interval affected by the specialPoint
        gMeshH = multiplyByH(gMesh, quadMesh.h, quadMesh.hTemp1, quadMesh.hTemp2, quadMesh.index, intMeshSize);
        
%       then call 'assignSecondPart' for subsequent stages
        stageEval1 = assignSecondPart (transformedMesh, intMeshSize, reuseIndex, t_n, y_n, 1, K, h);
        stageEval2 = assignSecondPart(transformedMesh, intMeshSize, reuseIndex, t_n, y_n, 2, K, h);
        
%       After the assignments are done, perform the numerical integration
%       on the 'assigned' meshes. Add up the first and second parts to get
%       the full integral to use later
        I_numPartial = integrate(1, gMeshH, evalFirstPart);
        I_numStage1 = integrate(reuseIndex, gMeshH, stageEval1);
        I_numStage2 = integrate(reuseIndex, gMeshH, stageEval2);
        
        K(2) = F(y_n + h.*K(1), I_numPartial + I_numStage1);
        y_next = y_n + h.*(1/2).*(K(1) + K(2));
        SOL(i+1, 1) = t_right;
        SOL(i+1, 2) = y_next;
        SOL(i, 3:4) = linTerp(t_n, t_right, y_n, y_next); %straight line between  y_n and y_next
        
%       evaluate K(1) of the endpoint t_next, but now with the final stage interpolant given by
%       the RK method as the function definition
        K(1) = F(y_next, I_numPartial + I_numStage2);
        
%       advance t and y
        t_n = t_right;
        y_n = y_next;
    end

    function I_num = integrate(leastIndex, gMeshH, assignmentValues)
%   integrates on the mesh that is supplied to it from the least index up to the highest
%   index. 
        maxIndex = leastIndex + size(assignmentValues, 1)-1;
        while size(intCoeffs, 2) < intMeshSize
                intCoeffs = [intCoeffs, 2, -1, 2]; %grow the intCoeffs row vector as necessary
        end
%       idea here is to perform a dot product, but only for the mesh points
%       needed. For example, to find the integral from -\infty to t_0, only
%       perform the dot product of the first reuseIndex - 1 coefficients,
%       and the corresponding function evaluations found in
%       'assignmentValues'
        I_num = intCoeffs(1, leastIndex: maxIndex) * (assignmentValues.*gMeshH(leastIndex: maxIndex));
    end
    SOL(numSteps + 1, 1) = right; %do this or 'deval' will shit itself because of rounding errors
    sol.x = SOL(:, 1);
    sol.y = SOL(:, 2:4);
    sol.hist = phi;
    sol.interpolantOrder = 2;
end

function [evalFirstPart, reuseIndex] = assignFirstPart(SOL, transformedMesh, intMeshSize, t_0, t_n, phi, h)
%   take as input a mesh from (-\infty, t_right) and assign the correct evaluation to each point of
%   the transformed mesh falling between -\infty and t_n. It also outputs the
%   index in this mesh of the first point which is > t_n and < t_right

%   fnDefs has format [t_n, yn, ypn, a_0, a_1] where a_i is a
%   coefficient of the linear interpolant on the interval starting with t_n
%   t_left is the left endpoint of current interval we are solving in.
       
       evalFirstPart = zeros(intMeshSize, 1);
       reuseIndex = 0;
       for loopIndex = 1 : intMeshSize
            currPoint = transformedMesh(loopIndex);
            if currPoint <= t_0 %if the integration mesh point is in the history, evaluate it as such
                evalFirstPart(loopIndex) = phi(currPoint);
            elseif currPoint > t_0 && currPoint <= t_n 
%               i.e. if the transformed point falls somewhere in the domain of
%               the already computed solution
                interval = floor((currPoint - t_0)/h)+1;
                evalFirstPart(loopIndex) =  SOL(interval, 3 : 4 ) * [1; currPoint];
            elseif currPoint > t_n %apply the stage interpolant for the current interval
                   reuseIndex = loopIndex; 
                break
            end
       end
       if reuseIndex ~= 0 %only return the relevant part of EvalFirstPart
            evalFirstPart = evalFirstPart(1: reuseIndex - 1);
       end
end

function stageEval = assignSecondPart(transformedMesh, intMeshSize, reuseIndex, t_n, y_n, stageNumber, K, h)
%   This function takes as input a mesh from (-\infty, t_right) and the index
%   of the first point which is > t_n and < t_right. Depending on the stage,
%   use a different function to evaluate these points.
    stageEval = zeros(intMeshSize - reuseIndex +1, 1); %allocate space for the stageEval vector
    switch stageNumber
        case 1
            stageEval = y_n + (K(1)).*(transformedMesh(reuseIndex : intMeshSize) - t_n);
        case 2
            b2 = (1/(2.*h)) .* ((transformedMesh(reuseIndex : intMeshSize) - t_n).^2);
            b1 = (transformedMesh(reuseIndex : intMeshSize) - t_n) - b2;
            stageEval= y_n + ( b1.*K(1) + b2.* K(2));
    end 
end

function gMeshH = multiplyByH(gMesh, hDefault, hTemp1, hTemp2, spIndex, gMeshSize)
%   This function multiplies each entry in the gMesh by its appropriate
%   step size. This is slightly complicated since the presence of the
%   special point creates six mesh points with a different spacing
    gMeshH = zeros(gMeshSize, 1);
    if spIndex < 0 %there is no special point, it's already in the mesh, and so every mesh point gets the defeault h
        gMeshH = hDefault/3 .* gMesh;
    else
        gMeshH(1 : spIndex - 1) = hDefault./3.*gMesh(1 : spIndex-1);
        gMeshH(spIndex : spIndex + 2) = hTemp1./3 .*gMesh(spIndex : spIndex + 2);
        gMeshH(spIndex+3 : spIndex + 5) = hTemp2./3 .*gMesh(spIndex+3 : spIndex + 5);
        gMeshH(spIndex + 6 : gMeshSize) = hDefault./3.*gMesh(spIndex + 6 : gMeshSize);
    end
end

function newMesh = getMesh_v4(hQuadrature, specialPoint)
    numPoints = floor((1/hQuadrature))+1; %determines a uniform mesh size which is slightly less than the one pre-calculated in the preamble of the script
    h = 1 ./ numPoints;
    %divide h by 4 once and save these values to reduce rounding error
    h_1 = h./4;
    h_2 = h./2;
    h_3 = 3*h_1;
    
    test = numPoints * specialPoint;
    spInterval = 0; %this will contain the interval of the mesh which contains the special point
    index = 0;
    if  test == floor(test) %check to see if the special point is already part of the mesh, If so, make the mesh without the special point
        mesh = zeros(3 * (numPoints), 1);
        index = -1;
        spInterval = -1;
        hTemp1 = -1;
        hTemp2 = -1;
        for i = 0 : 3: 3*(numPoints)-1
            mesh(i+1, 1) = h*(floor(i/3)) + h_1;
            mesh(i+2, 1) = h*(floor(i/3)) + h_2;
            mesh(i+3, 1) = h*(floor(i/3)) + h_3;
        end
    else 
        mesh = zeros(3*numPoints +3, 1);
        spInterval = floor(test); %this places the special point in its proper interval
        counter = 0;
        interval = 0;
        while interval < spInterval && interval < numPoints-1
            mesh(counter+1, 1) = h * interval + h_1;
            mesh(counter+2, 1) = h * interval + h_2;
            mesh(counter+3, 1) = h * interval + h_3;
%             mesh(counter+1, 2) = mesh(counter+1, 1) - h*interval;
%             mesh(counter+2, 2) = mesh(counter+2, 1) - mesh(counter+1, 1);
%             mesh(counter+3, 2) = mesh(counter+3, 1) - mesh(counter+2, 1);
            counter = counter + 3;
            interval = interval + 1;
        end
        index = counter+1; %save the index of first meshpoint in the left 'specialPoint' interval. 
        hTemp1 = specialPoint - h*interval;
        mesh(counter+1, 1) = h * interval + hTemp1 /4;
        mesh(counter+2, 1) = h * interval + hTemp1 /2;
        mesh(counter + 3, 1) = h * interval + 3*hTemp1 / 4;
        
%         mesh(counter+1, 2) = mesh(counter+1, 1) - h*interval;
%         mesh(counter+2, 2) = mesh(counter+2, 1) - mesh(counter+1, 1);
%         mesh(counter+3, 2) = mesh(counter+3, 1) - mesh(counter+2, 1);
         counter = counter + 3;
        %spInterval = interval; % save the mesh interval later to add to the mesh struct in the output
        interval = interval + 1;
        
        hTemp2 = h*interval - specialPoint;
        mesh(counter+1, 1) = specialPoint + hTemp2 /4;
        mesh(counter+2, 1) = specialPoint + hTemp2 /2;
        mesh(counter + 3, 1) = specialPoint + 3*hTemp2 / 4;
%         mesh(counter+1, 2) = mesh(counter+1, 1) - h*interval;
%         mesh(counter+2, 2) = mesh(counter+2, 1) - mesh(counter+1, 1);
%         mesh(counter+3, 2) = mesh(counter+3, 1) - mesh(counter+2, 1);
        counter = counter + 3;
        
        while interval < numPoints
             mesh(counter+1, 1) = h * interval + h_1;
            mesh(counter+2, 1) = h * interval + h_2;
            mesh(counter+3, 1) = h * interval + h_3;
%             mesh(counter+1, 2) = mesh(counter+1, 1) - h*interval;
%             mesh(counter+2, 2) = mesh(counter+2, 1) - mesh(counter+1, 1);
%             mesh(counter+3, 2) = mesh(counter+3, 1) - mesh(counter+2, 1);
            counter = counter + 3;
            interval = interval + 1;
        end
    end
    newMesh.index = index;
    newMesh.numPoints = numPoints;
    newMesh.mesh = mesh;
    newMesh.interval = spInterval;
    newMesh.h = h;
    newMesh.hTemp1 = hTemp1;
    newMesh.hTemp2 = hTemp2;
end