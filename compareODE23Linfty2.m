function [sol, Y, times] = compareODE23Linfty2(F, phi, tau, j, left, right)
    hist = strrep(char(phi), '@(t)', '');
    rhs = strrep(char(F), '@(x,d)', '');
    %phi is a function handle, multiply it against the gamma distribution
    %function
    a = j/tau;
    y_0 = zeros(j+1, 1);
    y_0(1) = phi(0);
    t = sym('t');
    Y=0;
    for i = 2 : j+1
    %the contents of this for loop are for finding the vector of initial
    %values for ode23
%        g = @(t) exp(-a.*t) .* t.^(i-2);
%        h = @(t) g(t).*phi(-t);
%        initialVal = integral(h, 0, Inf,'RelTol',1e-15);
       g(t) = exp(-a*t) .* t^(i-2) ;
       h(t) = g.*phi(-t);
       initialVal = eval(int(h, t, 0, Inf));
       y_0(i) = (a.^(i-1))./gamma(i-1) .*initialVal/a;
    end
    
    %these next two lines compute the ode23 and RK2Solver solutions
    options = odeset('RelTol',1e-15);
    sol = ode23(@(t,y) odefun2(t, y, j, tau), [left, right], y_0, options);
    %----------------------------------------------------------------------
    %This section is strictly for visualizing the graph of the solution
    %itself compared to the actual solution 
%     Y = RK2Solver_v2(F, phi, tau, j, left, right, steps);
%     
%     queryPoints = linspace(left, right, numQueryPoints);
%     
%     ode23Eval = deval(sol, queryPoints, 1)';
%     RKEval = evalSol(Y, queryPoints');
%     
%     ERR = zeros(10, 2); %premake the vector which will contain the error with step size
%---------------------------------------------------------------------------     
    times = zeros(10, 2);
    steps = ceil((right - left)/0.1);
    for i = 0 : 9
        if i ~=0
            steps = 2*steps;
        end
        disp(['steps: ', num2str(steps)])
        h = (right - left)/(steps);
        tic;
        Y = RK2Solver_v2(F, phi, tau, j, left, right, steps);
        t = toc;
        times(i+1, 2) = t; 
        times(i+1, 1) = steps;
        MAX = abs(evalSol(Y, left) - deval(sol, left, 1)); %set the first point to be the maximum
        time = left+h;
        while time <= right %loop through the rest of the mesh points
            err = abs(evalSol(Y, time) - deval(sol, time, 1));
            if  err > MAX 
                MAX = err;
            end
            time = time + h;
        end
        ERR(i + 1, 1) = log(h)./log(10); %fill first column with the step size
        ERR(i + 1, 2) = log(MAX)./log(10); %fill the second column with the error
    end

    FIG = figure();
    scatter(ERR(:, 1),ERR(:,2), '.', 'DisplayName', 'ode23');
    %scatter(queryPoints', ode23Eval(:,1), '.','DisplayName', 'ode23');
%     hold on 
%     plot(queryPoints', RKEval(:,1),'DisplayName', 'RK2Solver');
%     legend;
%     hold off 
    
    title(['RK2 Solver and ode23 error for history function ', hist, ', RHS: ', rhs, ', j= ', num2str(j), ', tau= ', num2str(tau)])
    xlabel('log_{10} of step size h');
    ylabel('log_{10} of max |y(t) - u_{h}(t)|');
    fileName = ['j=', num2str(j), 'tau=', num2str(tau), '.jpg'];
    path = 'C:\Users\pgill\Desktop\RK2 plots';
    saveas(FIG, fullfile(path, fileName), 'jpg');
end