% Demand estimation under Bertrand-Nash competition via LP
testcase = set_prices_shock_testcase();

% save theta as:
% theta_i=(theta_i0/theta_i-i,theta_ii/theta_i-i,theta_i3/theta_i-i)

%% idealized setting M-NE
firm_i = 1;
% should define a grid, for now take p_i=0, p_i^*,\overline(p)
[X1, Y1, X1_constr, Y1_constr] = compute_error_grid(testcase.p1,testcase.p2,testcase.xi1,...
    firm_i, testcase.grid,testcase.bound,testcase.p_max);

% transform X and Y such that they are in standard form and solve LP
x = solve_LP(X1,Y1);

% find polyhedron of optimal solutions
t1_opt=x(4); % returns 4th entry
t_opt_vec=t1_opt*ones(size(X1,1),1);
add_zeros = size(Y1_constr,1)-size(Y1,1);
t_opt_vec = [t_opt_vec; zeros(add_zeros,1)];

% Construct Polyhedron from constraints
P1_theta = Polyhedron(X1_constr, (Y1_constr+t_opt_vec));
% P1_theta = P1_theta.computeVRep();
% P1_theta = P1_theta.minHRep();

firm_i = 2;
% should define a grid, for now take p_i=0, p_i^*,\overline(p)
[X2, Y2, X2_constr, Y2_constr] = compute_error_grid(testcase.p1,testcase.p2,testcase.xi1,...
    firm_i, testcase.grid,testcase.bound,testcase.p_max);

% transform X and Y such that they are in standard form and solve LP
x = solve_LP(X2,Y2);

% find polyhedron of optimal solutions
t2_opt=x(4); % returns 4th entry
t_opt_vec=t2_opt*ones(size(X2,1),1);
add_zeros = size(Y2_constr,1)-size(Y2,1);
t_opt_vec = [t_opt_vec; zeros(add_zeros,1)];
% t_opt = 0;
% t_opt_vec = t_opt*ones(size(X_constr,1),1);

% Construct Polyhedron from constraints
P2_theta = Polyhedron(X2_constr, (Y2_constr+t_opt_vec));
% P2_theta = P2_theta.minHRep();

%
P1_theta.contains([testcase.theta1(1,1)/testcase.theta1(1,3); testcase.theta1(1,2)/testcase.theta1(1,3);...
    testcase.theta1(1,4)/testcase.theta1(1,3)])

P2_theta.contains([testcase.theta2(1,1)/testcase.theta2(1,3); testcase.theta2(1,2)/testcase.theta2(1,3);...
    testcase.theta2(1,4)/testcase.theta2(1,3)])

% idealized setting inverse VI

[x,fval] = conic_opti(testcase.p1,testcase.p2,testcase.xi1,testcase.xi1,testcase.bound,testcase.p_max);

% use if we use classic conic optimization from Bertsimas
N = size(x,1)-2*4; 
N = N/3;
x(2*N+1:2*N+8)
P1 = [x(2*N+1,1) x(2*N+2,1) x(2*N+3,1) x(2*N+4,1)];
P2 = [x(2*N+5,1) x(2*N+6,1) x(2*N+7,1) x(2*N+8,1)];
P1_tilde = [P1(1,1)/P1(1,3) P1(1,2)/P1(1,3) P1(1,4)/P1(1,3)];
P2_tilde = [P2(1,1)/P2(1,3) P2(1,2)/P2(1,3) P2(1,4)/P2(1,3)];

% use if we use normalized conic optimization from Bertsimas
% N = size(x,1)-2*3;
% N = N/3;
% x(2*N+1:2*N+6)
% P1_tilde = [x(2*N+1,1) x(2*N+2,1) x(2*N+3,1)];
% P2_tilde = [x(2*N+4,1) x(2*N+5,1) x(2*N+6,1)];


% check whether the polyhedral solution set obtained from the M-NE-L_infty
% approach contains the conic solution
P1_theta.contains(transpose(P1_tilde))
P2_theta.contains(transpose(P2_tilde))

%%
fid = fopen('Figures/Stats_idealized.txt', 'at');
fprintf(fid, '\n M=%f \n true_theta1=%s true_theta2=%s \n conic_theta1=%s conic_theta2=%s \n lp_theta1 = %s \n lp_theta2=%s \n eps1=%f eps2=%f \n xi=%s \n p1=%s \n p2=%s \n',...
    size(testcase.p1,1), mat2str(testcase.theta1), mat2str(testcase.theta2), mat2str(x(2*N+1:2*N+4)), mat2str(x(2*N+5:2*N+8)),...
    mat2str(P1_theta.V), mat2str(P2_theta.V),...
    t1_opt, t2_opt,...
    mat2str(testcase.xi1), mat2str(testcase.p1), mat2str(testcase.p2));
fclose(fid);

%% Plot the polyhedral solution set, the true underlying preference vector
% (blue) and the conic solution (black)

fig = figure;

% Figure of Firm 1
subplot(1,2,1)
% plot solution set with xi
sol = P1_theta.plot('color',[0,0,0]+0.5);
% plot true underlying paramter
hold on
P_star = plot3(testcase.theta1(1,1)/testcase.theta1(1,3),testcase.theta1(1,2)/testcase.theta1(1,3),...
    testcase.theta1(1,4)/testcase.theta1(1,3), 'bo', 'MarkerSize', 9, 'MarkerFaceColor', 'b');
% plot true conic solution
hold on
P_conic = plot3(P1_tilde(1,1),P1_tilde(1,2),P1_tilde(1,3), 'ko', 'MarkerSize', 9, 'MarkerFaceColor', 'k');
xlabel('$\frac{\theta_{i,0}}{\theta_{i,-i}}$','interpreter','latex','FontSize',40);
ylabel('$\frac{\theta_{i,i}}{\theta_{i,-i}}$','interpreter','latex','FontSize',40);
zlabel('$\frac{\theta_{i,3}}{\theta_{i,-i}}$','interpreter','latex','FontSize',40);
ax = gca;
ax.FontSize = 20;
% legend([dynamic static theta_star],{'dynamic','static','$\theta_1^*$'},'Interpreter','latex','FontSize',11)
% legend([sol P_star P_conic],{'$\mathcal{P}_1^*$','$P_1^*$', '$P_{conic,1}$'},'Interpreter','latex','FontSize',14)
view(2)

% Figure Firm 2
subplot(1,2,2)
% fill in feasible region
sol = P2_theta.plot('color',[0,0,0]+0.5);
% plot ture underlying paramter
hold on
P_star = plot3(testcase.theta2(1,1)/testcase.theta2(1,3),testcase.theta2(1,2)/testcase.theta2(1,3),...
    testcase.theta2(1,4)/testcase.theta2(1,3), 'bo', 'MarkerSize', 9, 'MarkerFaceColor', 'b');
% plot true conic solution
hold on
P_conic = plot3(P2_tilde(1,1),P2_tilde(1,2),P2_tilde(1,3), 'ko', 'MarkerSize', 9, 'MarkerFaceColor', 'k');
xlabel('$\frac{\theta_{i,0}}{\theta_{i,-i}}$','interpreter','latex','FontSize',40);
ylabel('$\frac{\theta_{i,i}}{\theta_{i,-i}}$','interpreter','latex','FontSize',40);
zlabel('$\frac{\theta_{i,3}}{\theta_{i,-i}}$','interpreter','latex','FontSize',40);
ax = gca;
ax.FontSize = 20;
% legend([sol P_star P_conic],{'$\mathcal{P}_2^*$','$P_2^*$, $P_{conic,2}$'},'Interpreter','latex','FontSize',14)
% legend([sol P_star P_conic],{'$\mathcal{P}_1^*$','$P_1^*$', '$P_{conic,1}$'},'Interpreter','latex','FontSize',14)
view(2)

matlab2tikz('./Figures/polyhedron_M=5.tikz');

%% estimated marginal revenue
% fixed ξ to be its median value over the dataset, and fixed the other firm’s 
% price to be the price observed for this median value.
xi_med = median(testcase.xi1);
p1_med = median(testcase.p1);
p2_med = median(testcase.p2);
p = linspace(0,max(max(testcase.p1),max(testcase.p2)),50);

fig = figure();
subplot(1,2,1)
% true demand firm 1 (fixed to median)
true = plot(p,(testcase.theta1(1,1)+2*testcase.theta1(1,2)*p+testcase.theta1(1,3)*p2_med+...
    testcase.theta1(1,4)*xi_med),'-o','Color','b','MarkerSize',4,'MarkerFaceColor','b');
hold on 
conic = plot(p,(P1(1,1)+2*P1(1,2)*p+P1(1,3)*p2_med+P1(1,4)*xi_med),'-s','Color','k','MarkerSize',5,'MarkerFaceColor','k');
for j=1:size(P1_theta.V,1)-1
    y1 = P1_theta.V(j,1)+2*P1_theta.V(j,2)*p+p2_med+P1_theta.V(j,3)*xi_med;
    y2 = P1_theta.V(j+1,1)+2*P1_theta.V(j+1,2)*p+p2_med+P1_theta.V(j+1,3)*xi_med;
    lp = plot(p,y1,'-*','MarkerSize',4,'color',[0,0,0]+0.5);
    hold on
    lp = plot(p,y2,'color',[0,0,0]+0.5);
    hold on
    patch([p fliplr(p)], [y1 fliplr(y2)], 'k', 'FaceAlpha',0.2)
end
xlabel('Firm 1 Price','FontSize',18,'FontWeight','bold')
ylabel('Marginal Revenue','FontSize',18,'FontWeight','bold')
grid on
pbaspect([1 2 1])
legend([true conic lp],{'true MR','inverse VI fit','our'},'Interpreter','latex','FontSize',11)


subplot(1,2,2)
% true demand firm 1 (fixed to median)
true = plot(p,(testcase.theta2(1,1)+2*testcase.theta2(1,2)*p+testcase.theta2(1,3)*p1_med+...
    testcase.theta2(1,4)*xi_med),'-o','Color','b','MarkerSize',4,'MarkerFaceColor','b');
hold on 
conic = plot(p,(P2(1,1)+2*P2(1,2)*p+P2(1,3)*p1_med+P2(1,4)*xi_med),'-s','Color','k','MarkerSize',5,'MarkerFaceColor','k');
for j=1:size(P2_theta.V,1)-1
    y1 = P2_theta.V(j,1)+2*P2_theta.V(j,2)*p+p1_med+P2_theta.V(j,3)*xi_med;
    y2 = P2_theta.V(j+1,1)+2*P2_theta.V(j+1,2)*p+p1_med+P2_theta.V(j+1,3)*xi_med;
    lp = plot(p,y1,'-*','MarkerSize',4,'color',[0,0,0]+0.5);
    hold on
    lp = plot(p,y2,'color',[0,0,0]+0.5);
    hold on
    patch([p fliplr(p)], [y1 fliplr(y2)], 'k', 'FaceAlpha',0.1)
end
xlabel('Firm 2 Price', 'FontSize',18,'FontWeight','bold')
ylabel('Marginal Revenue','FontSize',18,'FontWeight','bold')
grid on
pbaspect([1 2 1])
%%
matlab2tikz('./Figures/marginal_revenue_M=10.tikz');
%% noisy setting
% to do: Gedanken machen über epsilon (rationality), Werte von Theta und
% warum conic estimation fast Null ist

firm_i = 1;
% should define a grid, for now take p_i=0, p_i^*,\overline(p)
[X1, Y1, X1_constr, Y1_constr] = compute_error_grid(testcase.p1_noisy,testcase.p2_noisy,testcase.xi_noisy,...
    firm_i, testcase.grid,testcase.bound,testcase.p_max);


% transform X and Y such that they are in standard form and solve LP
x = solve_LP(X1,Y1);

% find polyhedron of optimal solutions
t_opt=x(4); % returns 4th entry
t_opt_vec=t_opt*ones(size(X1,1),1);
add_zeros = size(Y1_constr,1)-size(Y1,1);
t_opt_vec = [t_opt_vec; zeros(add_zeros,1)];
% t_opt = 0;
% t_opt_vec = t_opt*ones(size(X_constr,1),1);

% Construct Polyhedron from constraints
P1_theta = Polyhedron(X1_constr, (Y1_constr+t_opt_vec));


firm_i = 2;
% should define a grid, for now take p_i=0, p_i^*,\overline(p)
[X2, Y2, X2_constr, Y2_constr] = compute_error_grid(testcase.p1_noisy,testcase.p2_noisy,testcase.xi_noisy,...
    firm_i, testcase.grid,testcase.bound,testcase.p_max);


% transform X and Y such that they are in standard form and solve LP
x = solve_LP(X2,Y2);

% find polyhedron of optimal solutions
t_opt=x(4); % returns 4th entry
t_opt_vec=t_opt*ones(size(X2,1),1);
add_zeros = size(Y2_constr,1)-size(Y2,1);
t_opt_vec = [t_opt_vec; zeros(add_zeros,1)];
% t_opt = 0;
% t_opt_vec = t_opt*ones(size(X_constr,1),1);

% Construct Polyhedron from constraints
P2_theta = Polyhedron(X2_constr, (Y2_constr+t_opt_vec));
%%
% check if true parameters lie in polyhedral solution
contains_P1 = P1_theta.contains([testcase.theta1(1,1)/testcase.theta1(1,3); testcase.theta1(1,2)/testcase.theta1(1,3);...
    testcase.theta1(1,4)/testcase.theta1(1,3)])

contains_P2 = P2_theta.contains([testcase.theta2(1,1)/testcase.theta2(1,3); testcase.theta2(1,2)/testcase.theta2(1,3);...
    testcase.theta2(1,4)/testcase.theta2(1,3)])

%% noisy setting inverse VI

[x,fval] = conic_opti(testcase.p1_noisy,testcase.p2_noisy,testcase.xi_noisy,testcase.xi_noisy,testcase.bound,testcase.p_max);
% % for the normalized theta version 
% N = size(x,1)-2*3;
% N = N/3;
% x(2*N+1:2*N+6)
% P1_tilde = [x(2*N+1,1) x(2*N+2,1) x(2*N+3,1)];
% P2_tilde = [x(2*N+4,1) x(2*N+5,1) x(2*N+6,1)];
% 
% % check whether the polyhedral solution set obtained from the M-NE-L_infty
% % approach contains the conic solution
% P1_theta.contains(transpose(P1_tilde))
% P2_theta.contains(transpose(P2_tilde))

% for the original theta version
N = size(x,1)-2*4;
N = N/3;
x(2*N+1:2*N+8)
P1 = [x(2*N+1,1) x(2*N+2,1) x(2*N+3,1) x(2*N+4,1)];
P2 = [x(2*N+5,1) x(2*N+6,1) x(2*N+7,1) x(2*N+8,1)];
P1_tilde = [P1(1,1)/P1(1,3) P1(1,2)/P1(1,3) P1(1,4)/P1(1,3)];
P2_tilde = [P2(1,1)/P2(1,3) P2(1,2)/P2(1,3) P2(1,4)/P2(1,3)];
fval
%%
% check whether the polyhedral solution set obtained from the M-NE-L_infty
% approach contains the conic solution
contains_conic_P1 = P1_theta.contains(transpose(P1_tilde))
contains_conic_P2 = P2_theta.contains(transpose(P2_tilde))

%%
fid = fopen('Figures/Stats_noisy.txt', 'at');
fprintf(fid, '\n M=%f \n true_theta1=%s true_theta2=%s \n conic_theta1=%s conic_theta2=%s \n lp_theta1 = %s \n lp_theta2=%s \n xi=%s \n xi1=%s \n xi2=%s \n p1=%s \n p2=%s \n theta_1=%f theta_2=%f theta_1_conic=%f theta_2_conic=%f \n',...
    size(testcase.p1,1), mat2str(testcase.theta1), mat2str(testcase.theta2), mat2str(x(2*N+1:2*N+4)), mat2str(x(2*N+5:2*N+8)),...
    mat2str(P1_theta.V), mat2str(P2_theta.V),...
    mat2str(testcase.xi_noisy), mat2str(testcase.xi1), mat2str(testcase.xi1),...
    mat2str(testcase.p1), mat2str(testcase.p2),...
    contains_P1, contains_P2, contains_conic_P1, contains_conic_P2);
fclose(fid);

%% Plot the polyhedral solution set, the true underlying preference vector
% (blue) and the conic solution (black)

fig = figure;

% Figure of Firm 1
subplot(1,2,1)
% plot solution set 
if size(P1_theta.V,1)==1
    sol = plot3(P1_theta.V(1,1),P1_theta.V(1,2),P1_theta.V(1,3),'o','MarkerSize', 5,... 
        'MarkerEdgeColor',[0,0,0]+0.5,'MarkerFaceColor',[0,0,0]+0.5);
else   
    sol = P1_theta.plot('color',[0,0,0]+0.5);
end
% plot true underlying paramter
hold on
P_star = plot3(testcase.theta1(1,1)/testcase.theta1(1,3),testcase.theta1(1,2)/testcase.theta1(1,3),...
    testcase.theta1(1,4)/testcase.theta1(1,3), 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
% plot true conic solution
hold on
P_conic = plot3(P1_tilde(1,1),P1_tilde(1,2),P1_tilde(1,3), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
xlabel('$\frac{\theta_{i,0}}{\theta_{i,-i}}$','interpreter','latex','FontSize',20);
ylabel('$\frac{\theta_{i,i}}{\theta_{i,-i}}$','interpreter','latex','FontSize',20);
zlabel('$\frac{\theta_{i,3}}{\theta_{i,-i}}$','interpreter','latex','FontSize',20);
% legend([dynamic static theta_star],{'dynamic','static','$\theta_1^*$'},'Interpreter','latex','FontSize',11)
% legend([sol P_star P_conic],{'$\mathcal{P}_1^*$','$P_1^*$', '$P_{conic,1}$'},'Interpreter','latex','FontSize',14)
grid on

% Figure Firm 2
subplot(1,2,2)
% plot solution
if size(P2_theta.V,1)==1
    sol = plot3(P2_theta.V(1,1),P2_theta.V(1,2),P2_theta.V(1,3),'o','MarkerSize', 5,... 
        'MarkerEdgeColor',[0,0,0]+0.5,'MarkerFaceColor',[0,0,0]+0.5);
else  
    sol = P2_theta.plot('color',[0,0,0]+0.5);
end
% plot ture underlying paramter
hold on
P_star = plot3(testcase.theta2(1,1)/testcase.theta2(1,3),testcase.theta2(1,2)/testcase.theta2(1,3),...
    testcase.theta2(1,4)/testcase.theta2(1,3), 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
% plot true conic solution
hold on
P_conic = plot3(P2_tilde(1,1),P2_tilde(1,2),P2_tilde(1,3), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
xlabel('$\frac{\theta_{i,0}}{\theta_{i,-i}}$','interpreter','latex','FontSize',20);
ylabel('$\frac{\theta_{i,i}}{\theta_{i,-i}}$','interpreter','latex','FontSize',20);
zlabel('$\frac{\theta_{i,3}}{\theta_{i,-i}}$','interpreter','latex','FontSize',20);
% legend([sol P_star P_conic],{'$\mathcal{P}_2^*$','$P_2^*$, $P_{conic,2}$'},'Interpreter','latex','FontSize',14)
% legend([sol P_star P_conic],{'$\mathcal{P}_1^*$','$P_1^*$', '$P_{conic,1}$'},'Interpreter','latex','FontSize',14)
grid on

%% estimated marginal revenue
% fixed ξ to be its median value over the dataset, and fixed the other firm’s 
% price to be the price observed for this median value.
xi_med = median(testcase.xi_noisy);
xi1_med = median(testcase.xi1);
xi2_med = median(testcase.xi2);
p1_med = median(testcase.p1_noisy);
p2_med = median(testcase.p2_noisy);
p = linspace(0,max(max(testcase.p1),max(testcase.p2)),50);

fig = figure();
subplot(1,2,1)
% true demand firm 1 (fixed to median)
true = plot(p,(testcase.theta1(1,1)+2*testcase.theta1(1,2)*p+testcase.theta1(1,3)*p2_med+...
    testcase.theta1(1,4)*xi1_med),'-o','Color','b','MarkerSize',4,'MarkerFaceColor','b');
hold on 
conic = plot(p,(P1(1,1)+2*P1(1,2)*p+P1(1,3)*p2_med+P1(1,4)*xi_med),'-s','Color','k','MarkerSize',5,'MarkerFaceColor','k');
if size(P1_theta.V,1)==1
    lp = plot(p,P1_theta.V(1,1)+2*P1_theta.V(1,2)*p+p2_med+P1_theta.V(1,3)*xi_med,'-*','MarkerSize',4,'color',[0,0,0]+0.5);
else
    for j=1:size(P1_theta.V,1)-1
        y1 = P1_theta.V(j,1)+2*P1_theta.V(j,2)*p+p2_med+P1_theta.V(j,3)*xi_med;
        y2 = P1_theta.V(j+1,1)+2*P1_theta.V(j+1,2)*p+p2_med+P1_theta.V(j+1,3)*xi_med;
        lp = plot(p,y1,'-*','MarkerSize',4,'color',[0,0,0]+0.5);
        hold on
        lp = plot(p,y2,'color',[0,0,0]+0.5);
        hold on
        patch([p fliplr(p)], [y1 fliplr(y2)], 'k', 'FaceAlpha',0.2)
    end
end
xlabel('Firm 1 Price','FontSize',18,'FontWeight','bold')
ylabel('Marginal Revenue','FontSize',18,'FontWeight','bold')
grid on

subplot(1,2,2)
% true demand firm 1 (fixed to median)
true = plot(p,(testcase.theta2(1,1)+2*testcase.theta2(1,2)*p+testcase.theta2(1,3)*p1_med+...
    testcase.theta2(1,4)*xi2_med),'-o','Color','b','MarkerSize',4,'MarkerFaceColor','b');
hold on 
conic = plot(p,(P2(1,1)+2*P2(1,2)*p+P2(1,3)*p1_med+P2(1,4)*xi_med),'-s','Color','k','MarkerSize',5,'MarkerFaceColor','k');
if size(P2_theta.V,1)==1
    lp = plot(p,P2_theta.V(1,1)+2*P2_theta.V(1,2)*p+p1_med+P2_theta.V(1,3)*xi_med,'-*','MarkerSize',4,'color',[0,0,0]+0.5);
else
    for j=1:size(P2_theta.V,1)-1
        y1 = P2_theta.V(j,1)+2*P2_theta.V(j,2)*p+p1_med+P2_theta.V(j,3)*xi_med;
        y2 = P2_theta.V(j+1,1)+2*P2_theta.V(j+1,2)*p+p1_med+P2_theta.V(j+1,3)*xi_med;
        lp = plot(p,y1,'-*','MarkerSize',4,'color',[0,0,0]+0.5);
        hold on
        lp = plot(p,y2,'color',[0,0,0]+0.5);
        hold on
        patch([p fliplr(p)], [y1 fliplr(y2)], 'k', 'FaceAlpha',0.2)
    end
end
xlabel('Firm 2 Price', 'FontSize',18,'FontWeight','bold')
ylabel('Marginal Revenue','FontSize',18,'FontWeight','bold')
grid on



%% 

% Bertsimas conic optimization (original version)
% look at note 10.01.21
function[x,fval] = conic_opti(p1,p2,xi1,xi2,bound,p_max)
    N = size(p1,1);
    % Initialize linear inequality constraints
    % x=(y1,y2,theta1,theta2,eps)
    % Note: we dont use xi

    % constraint y1^j>=0 for j=1,...,N
    Y1 = -eye(N);
    % constraint y2^j>=0 for j=1,...,N
    Y2 = -eye(N);
    % constraint  y1^j>=Demand1^j+theta11*p1^j for j=1,...,N
    Theta1 = [ones(N,1) 2*p1(:,1) p2(:,1) xi1(:,1)];
    % constraint  y2^j>=Demand2^j+theta22*p2^j for j=1,...,N
    Theta2 = [ones(N,1) 2*p2(:,1) p1(:,1) xi2(:,1)];
    % no constraints on epsilon
    Eps = zeros(N,N);

    Y1_add = p_max*eye(N);
    Y2_add = p_max*eye(N);
    Theta1_add = [];
    Theta2_add = [];
    Eps_add = -eye(N);
    for j=1:N
        Theta10 = -p1(j,1);
        Theta11 = -2*p1(j,1)^2;
        Theta12 = -p1(j,1)*p2(j,1);
        Theta13 = -p1(j,1)*xi1(j,1);
        Theta1_add = [Theta1_add;
                      Theta10 Theta11 Theta12 Theta13];
        Theta20 = -p2(j,1);
        Theta22 = -2*p2(j,1)^2;
        Theta21 = -p2(j,1)*p1(j,1);
        Theta23 = -p2(j,1)*xi2(j,1);
        Theta2_add = [Theta2_add;
                      Theta20 Theta22 Theta21 Theta23];

    end


    A = [Y1 zeros(N,N) zeros(N,4) zeros(N,4) Eps;
        zeros(N,N) Y2 zeros(N,4) zeros(N,4) Eps;
        Y1 zeros(N,N) Theta1 zeros(N,4) Eps;
        zeros(N,N) Y2 zeros(N,4) Theta2 Eps;
        Y1_add Y2_add Theta1_add Theta2_add Eps_add];

    b = zeros(size(A,1),1);

    % add constraints on the parameters -10<=theta_ik<=10 for i=1,2 k=0,1,2,3
    Theta_constr = [1 0 0 0;
                   -1 0 0 0;
                   0 1 0 0;
                   0 -1 0 0;
                   0 0 1 0;
                   0 0 -1 0;
                   0 0 0 1;
                   0 0 0 -1];
    other_params = zeros(8,N);

    A_constr = [A;
               other_params other_params Theta_constr zeros(8,4) other_params;
               other_params other_params zeros(8,4) Theta_constr other_params];

    % b_constr = [b; 10*ones(6*2,1)];
    b_constr = [b;
        bound;bound;0;bound;bound;0;bound;bound;bound;bound;0;bound;bound;0;bound;bound]; % assumption of concavity of U_i, 
    % hence theta_ii<=0

    % use quadprog as a solver
    % x = quadprog(H,f,A,b) minimizes 1/2*x'*H*x + f'*x subject to the restrictions A*x ≤ b
    % % x=(y1,y2,theta1,theta2,eps)

    % Initialize objective function
    f = zeros(1,size(A_constr,2));
    H = [zeros(N,N) zeros(N,N) zeros(N,4) zeros(N,4) zeros(N,N);
         zeros(N,N) zeros(N,N) zeros(N,4) zeros(N,4) zeros(N,N);
         zeros(4,N) zeros(4,N) zeros(4,4) zeros(4,4) zeros(4,N);
         zeros(4,N) zeros(4,N) zeros(4,4) zeros(4,4) zeros(4,N);
         zeros(N,N) zeros(N,N) zeros(N,4) zeros(N,4) 2*eye(N)];

    [x,fval] = quadprog(H,f,A_constr,b_constr); 
end

% Bertsimas conic optimization (normalized theta version)
% look at note 28.07.21
% function[x,fval] = conic_opti(p1,p2,xi1,xi2,bound,p_max)
%     N = size(p1,1);
%     % Initialize linear inequality constraints
%     % x=(y1,y2,theta1,theta2,eps)
%     % Note: we dont use xi
% 
%     % constraint y1^j>=0 for j=1,...,N
%     Y1 = -eye(N);
%     % constraint y2^j>=0 for j=1,...,N
%     Y2 = -eye(N);
%     % constraint  y1^j>=Demand1^j+theta11*p1^j for j=1,...,N
%     Theta1 = [ones(N,1) 2*p1(:,1) xi1(:,1)];
%     % constraint  y2^j>=Demand2^j+theta22*p2^j for j=1,...,N
%     Theta2 = [ones(N,1) 2*p2(:,1) xi2(:,1)];
%     % no constraints on epsilon
%     Eps = zeros(N,N);
% 
%     Y1_add = p_max*eye(N);
%     Y2_add = p_max*eye(N);
%     Theta1_add = [];
%     Theta2_add = [];
%     Eps_add = -eye(N);
%     for j=1:N
%         Theta10 = -p1(j,1);
%         Theta11 = -2*p1(j,1)^2;
%         Theta13 = -p1(j,1)*xi1(j,1);
%         Theta1_add = [Theta1_add;
%                       Theta10 Theta11 Theta13];
%         Theta20 = -p2(j,1);
%         Theta22 = -2*p2(j,1)^2;
%         Theta23 = -p2(j,1)*xi2(j,1);
%         Theta2_add = [Theta2_add;
%                       Theta20 Theta22 Theta23];
% 
%     end
% 
% 
%     A = [Y1 zeros(N,N) zeros(N,3) zeros(N,3) Eps;
%         zeros(N,N) Y2 zeros(N,3) zeros(N,3) Eps;
%         Y1 zeros(N,N) Theta1 zeros(N,3) Eps;
%         zeros(N,N) Y2 zeros(N,3) Theta2 Eps;
%         Y1_add Y2_add Theta1_add Theta2_add Eps_add];
% 
%     b = [zeros(N,1); zeros(N,1); -p2(:,1); -p1(:,1); 2*p1(:,1).*p2(:,1)];
% 
%     % add constraints on the parameters -10<=theta_ik<=10 for i=1,2 k=0,1,2,3
%     Theta_constr = [1 0 0;
%                    -1 0 0;
%                    0 1 0;
%                    0 -1 0;
%                    0 0 1;
%                    0 0 -1];
%     other_params = zeros(6,N);
% 
%     A_constr = [A;
%                other_params other_params Theta_constr zeros(6,3) other_params;
%                other_params other_params zeros(6,3) Theta_constr other_params];
% 
%     b_constr = [b; 10*ones(6*2,1)];
%     b_constr = [b;
%         bound;bound;0;bound;bound;bound;bound;bound;0;bound;bound;bound]; % assumption of concavity of U_i, hence theta_ii<=0
% 
%     % use quadprog as a solver
%     % x = quadprog(H,f,A,b) minimizes 1/2*x'*H*x + f'*x subject to the restrictions A*x ≤ b
%     % x=(y1,y2,theta1,theta2,eps)
% 
%     % Initialize objective function
%     f = zeros(1,size(A_constr,2));
%     H = [zeros(N,N) zeros(N,N) zeros(N,3) zeros(N,3) zeros(N,N);
%          zeros(N,N) zeros(N,N) zeros(N,3) zeros(N,3) zeros(N,N);
%          zeros(3,N) zeros(3,N) zeros(3,3) zeros(3,3) zeros(3,N);
%          zeros(3,N) zeros(3,N) zeros(3,3) zeros(3,3) zeros(3,N);
%          zeros(N,N) zeros(N,N) zeros(N,3) zeros(N,3) 2*eye(N)];
% 
%     [x,fval] = quadprog(H,f,A_constr,b_constr); 
% end


% function computes error e_i^j(p_i^j, theta_i)=U_i(p_i^j,p_-i^j,theta_i)-U_i(p_i^j*,p_-i^j*,theta_i)
% =(p_i^j-p_i^j*)theta_i0+((p_i^j)^2-(p_i^j*)^2)theta_i1 for all i=1,2 and
% all j=1,...,N for all p_i^j\in {0,overline(p)^j}
function [X,Y,X_constr,Y_constr] = compute_error_grid(p1,p2,xi, agent_i, grid,bound,p_max)
X = [];
Y = [];
for j=1:size(p1,1)
    grid_points = linspace(0,p_max,2^grid+1);
    if agent_i == 1
        for i=1:size(grid_points,2)
            s_0 = grid_points(i)-p1(j,1);
            s_i = grid_points(i)^2-p1(j,1)^2;
            s_3 = (grid_points(i)-p1(j,1))*xi(j,1);
            X = [X;
                 s_0 s_i s_3];
            s_minus_i = (grid_points(i)-p1(j,1))*p2(j,1);
            Y = [Y;
                -s_minus_i];
        end
    end
    if agent_i == 2
        for i=1:size(grid_points,2)
            s_0 = grid_points(i)-p2(j,1);
            s_i = grid_points(i)^2-p2(j,1)^2;
            s_3 = (grid_points(i)-p2(j,1))*xi(j,1);
            X = [X;
                 s_0 s_i s_3];
            s_minus_i = (grid_points(i)-p2(j,1))*p1(j,1);
            Y = [Y;
                -s_minus_i];
        end
    end
    % constrain parameter space Theta=\mathbb(R)^3 to [-10,10]^3
    X_constr = [];
    X_constr =[X;
              1 0 0;
              -1 0 0;
              0 1 0;
              0 -1 0;
              0 0 1;
              0 0 -1];
   Y_constr = [];
   Y_constr = [Y; bound; bound; bound; bound; bound; bound];

end
end

function [x] = solve_LP(X,Y)

% Augment X_constrained to rephrase as LP
% add column for variable t (we have inequalities of the form
% X*theta_i-t*1_{2^(N-1)}<=Y)
t_vec = [-ones((size(X,1)),1)];

A = [X, t_vec];
% add constraint -t<=0
A = [A; 0,0,0,-1];

% Augment Y_constrained to rephrase as LP: add 0 for constraint -t<=0
b = [Y; 0];

% Objective function 
% which is of the form c^T*\tilde(x)=0*theta_i,1+...+0*theta_i,3+1*t=t
f = [0, 0, 0, 1];

% Solve LP
% x=(theta_i,t)=(theta_i1,theta_i2,theta_i3,t)
x = linprog(f,A,b);
end