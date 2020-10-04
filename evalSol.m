function solVals = evalSol(sol, pointsToEval)
%A script that receives a set of points to evaluate the solution of the IVP
%at, and returns the solution's function evaluations at that set of points

    if ~issorted(pointsToEval)
       ME = MException('InputException:DataFormatting', ...
            'Vector of points to evaluate not in ascending order.');
        throw(ME)
    end

    m = size(pointsToEval, 1);
    n = size(sol.x, 1);
    solVals = zeros(m, 1);
    order = sol.interpolantOrder;
    solValsCounter = 1;
   
    for i = 1 : m
       currPoint = pointsToEval(i) ;
       if currPoint <= sol.x(1)
           solVals(i) = sol.hist(currPoint);
       elseif currPoint> sol.x(1) && currPoint <= sol.x(n)
           while solValsCounter <n
              if currPoint > sol.x(solValsCounter) && currPoint <= sol.x(solValsCounter+1) 
                  switch order
                      case {3, 4}
                        solVals(i) = sol.y(solValsCounter, 2:5) * [ 1 ; currPoint ; currPoint.*currPoint  ; currPoint.*currPoint.*currPoint ];
                      case {1, 2}
                        solVals(i) = sol.y(solValsCounter, 2:3) * [ 1 ; currPoint ];
                  end
                break
              else
                  solValsCounter = solValsCounter+1;
              end
           end
       else
           oops = MException('evalSol:QueryPointOutOfBounds', [' Query point ', num2str(currPoint), ' is not less than', num2str(sol.x(n))]);
           throw(oops)
       end
    end
end