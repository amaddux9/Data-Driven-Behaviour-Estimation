function testcase = set_prices_shock_testcase()
% set_prices_dynamic_testcase: sets the initial variable for testing the demand
% estimation under better response dynamics

%% Filename
testcase.filename = 'observed_price_vectors';
 
% Variables for dynamic case
testcase.bound = 10;
testcase.grid = 7; % grid=2 is equivalent to using the function compute_error

% Ex 1 
% theta=(theta_i0, theta_ii, theta_i-i theta_i3)
testcase.theta1 = [1 -1.2 0.5 1];
testcase.theta2 = [1 -1 0.3 1];


testcase.seed = 1; % use zero if eps is fixed % with 3 some overshooting
rng(testcase.seed)

N = 5;
testcase.p_max = 8; % upper bound on price
testcase.p1 = zeros(N,1);
testcase.p2 = zeros(N,1);
testcase.p1_noisy = zeros(N,1);
testcase.p2_noisy = zeros(N,1);
testcase.xi1 = round(normrnd(5,1.5,[N,1]),2); % use this also in the idealized setting where xi is common knowledge
testcase.xi2 = round(normrnd(5,1.5,[N,1]),2);
testcase.xi_noisy = (round(normrnd(5,1.5,[N,1]),2)+testcase.xi1+testcase.xi2)/3;
for j=1:N
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % linear demand function (idealized setting)
    syms p1 p2
    eqn1 = -(testcase.theta1(1,1)+2*testcase.theta1(1,2)*p1+testcase.theta1(1,3)*p2+...
        testcase.theta1(1,4)*testcase.xi1(j,1)) == 0;
    eqn2 = -(testcase.theta2(1,1)+testcase.theta2(1,3)*p1+2*testcase.theta2(1,2)*p2+...
        testcase.theta2(1,4)*testcase.xi1(j,1)) == 0;
    [A,B] = equationsToMatrix([eqn1, eqn2], [p1,p2]);
    p = linsolve(A,B);
    p = round(double(p),4);
    if p(1) <= 0
        p(1) = 0;
    elseif p(1) >= testcase.p_max
            p(1) = testcase.p_max; 
    end
    if p(2) <= 0
        p(2) = 0;
    elseif p(2) >= testcase.p_max
        p(2) = testcase.p_max; 
    end
    testcase.p1(j,1) = p(1);
    testcase.p2(j,1) = p(2);
    
    % linear demand function (noisy setting)
    syms p1 p2
    eqn1 = -(testcase.theta1(1,1)+2*testcase.theta1(1,2)*p1+testcase.theta1(1,3)*p2+...
        testcase.theta1(1,4)*testcase.xi1(j,1)) == 0;
    eqn2 = -(testcase.theta2(1,1)+testcase.theta1(1,3)*p1+2*testcase.theta2(1,2)*p2+...
        testcase.theta2(1,4)*testcase.xi2(j,1)) == 0;
    [A,B] = equationsToMatrix([eqn1, eqn2], [p1,p2]);
    p = linsolve(A,B);
    p = round(double(p),4);
    if p(1) <= 0
        p(1) = 0;
    elseif p(1) >= testcase.p_max
            p(1) = testcase.p_max; 
    end
    if p(2) <= 0
        p(2) = 0;
    elseif p(2) >= testcase.p_max
        p(2) = testcase.p_max; 
    end
    testcase.p1_noisy(j,1) = p(1);
    testcase.p2_noisy(j,1) = p(2);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for demand function specified in Bertsimas Sec. 8.1 
    % idealized setting
%     step = 0.001;
%     max_niter = 500;
%     p = projection_algo(testcase.p_max,xi,step,max_niter,testcase.theta1(1),testcase.theta1(2),testcase.theta1(3),testcase.theta1(4),...
%         testcase.theta2(1),testcase.theta2(2),testcase.theta2(3),testcase.theta2(4));
%     testcase.p1(j,1) = p(1);
%     testcase.p2(j,1) = p(2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% scale theta by theta_{i,-i} and reduce to three dimensional vector
% testcase.theta1 = [testcase.theta1(1)/testcase.theta1(3) testcase.theta1(2)/testcase.theta1(3) testcase.theta1(4)/testcase.theta1(3)];
% testcase.theta2 = [testcase.theta2(1)/testcase.theta2(3) testcase.theta2(2)/testcase.theta2(3) testcase.theta2(4)/testcase.theta2(3)];

end

%%
function [p] = projection_algo(p_max,xi,step,max_niter,t10,t11,t12,t13,t20,t22,t21,t23)
    p = rand(2,1)*0.45;
    k = 0;
    isConverged = false;
    while (k < max_niter && isConverged == false)
        k = k+1;
        f = [-(log(p(1))+1+2*t11*p(1)+t12*p(2)+t13*xi+t10); -(log(p(2))+1+2*t22*p(2)+t21*p(1)+t23*xi+t20)];
        p = p-step*f;
        if p(1) <= 0
            p(1) = 0;
        elseif p(1) >= p_max
            p(1) = p_max; 
        end
        if p(2) <= 0
            p(2) = 0;
        elseif p(2) >= p_max
            p(2) = p_max; 
        end
        if (abs(f(1))<0.0001 && abs(f(2))<0.001)
            isConverged = true;
        end
    end      
    
end


