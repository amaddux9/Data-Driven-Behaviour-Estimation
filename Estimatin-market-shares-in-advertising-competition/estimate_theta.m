% estimate theta of both firma in a duopoly in a market share competition
testcase = dynamic_testcase_pepsicola();
lambda = 0;

%%
% yearly data gasmi

i = 1; % coke 
% compute error matrix on a yearly basis ( not quarterly)
[X1_gasmi,Y1_gasmi,X1_rows] = compute_error_lanchester(i,testcase.bound,testcase.adco_year,testcase.adpe_year,...
    testcase.mco_year,lambda);

% solve a LP to find minimal epsilon
x = solve_LP(X1_gasmi,Y1_gasmi,X1_rows);
pos = size(X1_gasmi,2)+1;
eps1 = x(pos);
scale = 1.05;
enlarged_eps1 = scale*eps1;

% Compute solution poylhedron Theta_i^* = {P_i\in Theta| e_i(t+1,P_i)<=eps_i, t=1,...,K-1}
eps1_vec= [enlarged_eps1*ones(X1_rows,1);
            zeros(size(X1_gasmi,1)-X1_rows,1)]; 
% Construct Polyhedron from constraints
P1_gasmi = Polyhedron(X1_gasmi,Y1_gasmi+eps1_vec);
% P1_theta = P1_theta.minHRep(); %remove
% P1_theta.computeVRep();
% disp(P1_theta.isEmptySet())
% disp(P1_theta.isFullDim())

i = 2; % pepsi 
% compute error matrix on a yearly basis ( not quarterly)
[X2_gasmi,Y2_gasmi,X2_rows] = compute_error_sorger(i,testcase.bound,testcase.adco_year,testcase.adpe_year,...
    testcase.mco_year,lambda);

% solve a LP to find minimal epsilon
x = solve_LP(X2_gasmi,Y2_gasmi,X2_rows);
pos = size(X2_gasmi,2)+1;
eps2 = x(pos);
scale = 1.05;
enlarged_eps2 = scale*eps2;

% Compute solution poylhedron Theta_i^* = {P_i\in Theta| e_i(t+1,P_i)<=eps_i, t=1,...,K-1}
eps2_vec= [enlarged_eps2*ones(X2_rows,1);
            zeros(size(X2_gasmi,1)-X2_rows,1)]; 
% Construct Polyhedron from constraints
P2_gasmi = Polyhedron(X2_gasmi,Y2_gasmi+eps2_vec);

disp(eps1)
disp(eps2)
disp(P1_gasmi.V)
disp(P2_gasmi.V)


%%
% use when error function with diminishing returns in advertising is used
% Check results, i.e. whether the market share is between 0 and 1
share_esti_coke = [testcase.mco_year(1,1)]; % market share of Coke from Coke's perspective
share_esti_pepsi = [testcase.mco_year(1,1)]; % market share of Coke from Pepsi's perspective
share_esti_ols = [testcase.mco_year(1,1)];
for t=2:length(testcase.adco_year)
    share_esti_coke = [share_esti_coke; share_esti_coke(t-1,1)+P1_gasmi.V(1)*sqrt(testcase.adco_year(t,1))*(1-share_esti_coke(t-1,1))-...
        P1_gasmi.V(2)*sqrt(testcase.adpe_year(t,1))*share_esti_coke(t-1,1)];
    share_esti_pepsi = [share_esti_pepsi; share_esti_pepsi(t-1,1)+P2_gasmi.V(2)*sqrt(testcase.adco_year(t,1))*(1-share_esti_pepsi(t-1,1))-...
        P2_gasmi.V(1)*sqrt(testcase.adpe_year(t,1))*share_esti_pepsi(t-1,1)];
    share_esti_ols = [share_esti_ols;share_esti_ols(t-1,1)+0.0106*sqrt(testcase.adco_year(t,1))*(1-share_esti_ols(t-1,1))-...
        0.0072*sqrt(testcase.adpe_year(t,1))*share_esti_ols(t-1,1)];
end
disp(share_esti_coke)
disp(share_esti_pepsi)
disp(share_esti_ols)

fig = figure();
subplot(1,2,1)
% true market share
true = plot(testcase.mco_year,'-*','MarkerSize',2);
hold on
% estimated market share
coke = plot(share_esti_coke,'-*','MarkerSize',2);
hold on
% OLS market share
ols = plot(share_esti_ols,'-*','MarkerSize',2);
ylim([0,1])
legend([true coke ols],{'True market share','Market share from Coke perspective','OLS'}, 'interpreter','latex','FontSize',9);
subplot(1,2,2)
% true market share
plot(testcase.mco_year,'-*','MarkerSize',2);
hold on
% estimated market share
pepsi = plot(share_esti_pepsi,'-*','MarkerSize',2);
hold on
% OLS market share
ols = plot(share_esti_ols,'-*','MarkerSize',2);
ylim([0,1])
% legend([true pepsi ols],{'True market share','Coke share from Pepsi perspective','OLS'}, 'interpreter','latex','FontSize',9);

% Check whether the transisition probabilities are between 0 and 1 [Horsky]
% if k_ij and A_i is larger or equal than zero for all i,j\in{1,2} it suffices to check whether k_ij*max A_j(t)<=1
disp(P1_gasmi.V(1)*sqrt(max(testcase.adco_year)))
disp(P1_gasmi.V(2)*sqrt(max(testcase.adpe_year)))
disp(P2_gasmi.V(1)*sqrt(max(testcase.adpe_year)))
disp(P2_gasmi.V(2)*sqrt(max(testcase.adco_year)))

%%
% use when error function with diminishing returns in advertising is used
% and word-of-mouth advertising
% Check results, i.e. whether the market share is between 0 and 1
share_esti_coke = [testcase.mco_year(1,1)]; % market share of Coke from Coke's perspective
share_esti_pepsi = [testcase.mco_year(1,1)]; % market share of Coke from Pepsi's perspective
share_esti_ols = [testcase.mco_year(1,1)];
for t=2:length(testcase.adco_year)
    share_esti_coke = [share_esti_coke; share_esti_coke(t-1,1)+P1_gasmi.V(1)*sqrt(testcase.adco_year(t,1))*sqrt((1-share_esti_coke(t-1,1)))-...
        P1_gasmi.V(2)*sqrt(testcase.adpe_year(t,1))*sqrt(share_esti_coke(t-1,1))];
    share_esti_pepsi = [share_esti_pepsi; share_esti_pepsi(t-1,1)+P2_gasmi.V(2)*sqrt(testcase.adco_year(t,1))*sqrt(1-share_esti_pepsi(t-1,1))-...
        P2_gasmi.V(1)*sqrt(testcase.adpe_year(t,1))*sqrt(share_esti_pepsi(t-1,1))];
    share_esti_ols = [share_esti_ols;share_esti_ols(t-1,1)+0.0124*sqrt(testcase.adco_year(t,1))*sqrt(1-share_esti_ols(t-1,1))-...
            0.0105*sqrt(testcase.adpe_year(t,1))*sqrt(share_esti_ols(t-1,1))];
end
disp(share_esti_coke)
disp(share_esti_pepsi)
disp(share_esti_ols)

fig = figure();
subplot(1,2,1)
% true market share
true = plot(testcase.mco_year,'-*','MarkerSize',2);
hold on
% estimated market share
coke = plot(share_esti_coke,'-*','MarkerSize',2);
hold on
% OLS market share
ols = plot(share_esti_ols,'-*','MarkerSize',2);
ylim([0,1])
legend([true coke ols],{'True market share','Market share from Coke perspective','OLS'}, 'interpreter','latex','FontSize',9);
subplot(1,2,2)
% true market share
plot(testcase.mco_year,'-*','MarkerSize',2);
hold on
% estimated market share
pepsi = plot(share_esti_pepsi,'-*','MarkerSize',2);
hold on
% OLS market share
ols = plot(share_esti_ols,'-*','MarkerSize',2);
ylim([0,1])
% legend([true pepsi ols],{'True market share','Coke share from Pepsi perspective','OLS'}, 'interpreter','latex','FontSize',9);



% Check whether the transisition probabilities are between 0 and 1 [Horsky]
% if k_ij and A_i is larger or equal than zero for all i,j\in{1,2} it suffices to check whether k_ij*max A_j(t)<=1
disp(P1_gasmi.V(1)*sqrt(max(testcase.adco_year)))
disp(P1_gasmi.V(2)*sqrt(max(testcase.adpe_year)))
disp(P2_gasmi.V(1)*sqrt(max(testcase.adpe_year)))
disp(P2_gasmi.V(2)*sqrt(max(testcase.adco_year)))


%% Plot beliefs on market share evolution when Coke applies Lanchester model and Pepsi applies Sorger modle
% plotting all vertices when error is increased by 5%
share_esti_coke = [testcase.mco_year(1,1);
    testcase.mco_year(1,1);
    testcase.mco_year(1,1)]; % market share of Coke from Coke's perspective
share_esti_pepsi = [testcase.mco_year(1,1);
    testcase.mco_year(1,1);
    testcase.mco_year(1,1)]; % market share of Coke from Pepsi's perspective
share_esti_coke_L2 = [testcase.mco_year(1,1)]; % market share of Coke from Coke's perspective
share_esti_pepsi_L2 = [testcase.mco_year(1,1)]; % market share of Coke from Pepsi's perspective
share_esti_ols_coke = [testcase.mco_year(1,1)];
share_esti_ols_pepsi = [testcase.mco_year(1,1)];
for t=2:length(testcase.adco_year)
    for i=1:3
        % Lanchester model
        share_esti_coke(i,t) = share_esti_coke(i,t-1)+P1_gasmi.V(i,1)*sqrt(testcase.adco_year(t,1))*(1-share_esti_coke(i,t-1))-...
            P1_gasmi.V(i,2)*sqrt(testcase.adpe_year(t,1))*share_esti_coke(i,t-1);
        % Sorger model
        share_esti_pepsi(i,t) = share_esti_pepsi(i,t-1)+P2_gasmi.V(i,2)*sqrt(testcase.adco_year(t,1))*sqrt(1-share_esti_pepsi(i,t-1))-...
            P2_gasmi.V(i,1)*sqrt(testcase.adpe_year(t,1))*sqrt(share_esti_pepsi(i,t-1));
    end
    share_esti_coke_L2 = [share_esti_coke_L2; share_esti_coke_L2(t-1,1)+0.0782*sqrt(testcase.adco_year(t,1))*(1-share_esti_coke_L2(t-1,1))-...
        0.0519*sqrt(testcase.adpe_year(t,1))*share_esti_coke_L2(t-1,1)];
    share_esti_pepsi_L2 = [share_esti_pepsi_L2; share_esti_pepsi_L2(t-1,1)+0.0208*sqrt(testcase.adco_year(t,1))*sqrt(1-share_esti_pepsi_L2(t-1,1))-...
        0.0158*sqrt(testcase.adpe_year(t,1))*sqrt(share_esti_pepsi_L2(t-1,1))];
    share_esti_ols_coke = [share_esti_ols_coke;share_esti_ols_coke(t-1,1)+0.0106*sqrt(testcase.adco_year(t,1))*(1-share_esti_ols_coke(t-1,1))-...
            0.0072*sqrt(testcase.adpe_year(t,1))*share_esti_ols_coke(t-1,1)];
    share_esti_ols_pepsi = [share_esti_ols_pepsi;share_esti_ols_pepsi(t-1,1)+0.0124*sqrt(testcase.adco_year(t,1))*sqrt(1-share_esti_ols_pepsi(t-1,1))-...
            0.0105*sqrt(testcase.adpe_year(t,1))*sqrt(share_esti_ols_pepsi(t-1,1))];
end

% use for patch function
t = linspace(1,size(share_esti_pepsi,2),size(share_esti_pepsi,2));

fig = figure();
subplot(1,2,1)
% true market share
true = plot(testcase.mco_year,'-o','Color','b','MarkerSize',4,'MarkerFaceColor','b');
hold on
% estimated market share
coke = plot(share_esti_coke(1,:),'-*','MarkerSize',4,'color',[0,0,0]+0.5);
hold on
plot(share_esti_coke(2,:),'-*','MarkerSize',4,'color',[0,0,0]+0.5);
hold on
plot(share_esti_coke(3,:),'-*','MarkerSize',4,'color',[0,0,0]+0.5);
hold on
% plot grey region
% hold on
% patch([t fliplr(t)], [share_esti_coke(1,:) fliplr(share_esti_coke(2,:))], 'k', 'FaceAlpha',0.2)
hold on
patch([t fliplr(t)], [share_esti_coke(2,:) fliplr(share_esti_coke(2,:))], 'k', 'FaceAlpha',0.2)
hold on
patch([t fliplr(t)], [share_esti_coke(3,:) fliplr(share_esti_coke(1,:))], 'k', 'FaceAlpha',0.2)
hold on
% L2 market share
L2 = plot(share_esti_coke_L2,'-*','MarkerSize',4,'color','r');
hold on
% OLS market share
ols = plot(share_esti_ols_coke,'-s','Color','k','MarkerSize',5,'MarkerFaceColor','k');
ylim([0,1])
xticks(1:4:18)
xticklabels({'1968','1972','1976','1980','1984'})
ax = gca;
ax.FontSize = 20;
xlabel('Years','FontSize',18,'FontWeight','bold')
ylabel('Cokes beliefes on market share evolution','FontSize',18,'FontWeight','bold')
%legend([true coke L2 ols],{'True market share','M-BRD: Cokes perspective','M-BR-$L_2$-D: Cokes perspective','OLS'}, 'interpreter','latex','FontSize',9);
% pbaspect([2 1 1])
grid on

subplot(1,2,2)
% true market share
plot(testcase.mco_year,'-o','Color','b','MarkerSize',4,'MarkerFaceColor','b');
hold on
% estimated market share
pepsi = plot(share_esti_pepsi(1,:),'-*','MarkerSize',4,'color',[0,0,0]+0.5);
hold on
plot(share_esti_pepsi(2,:),'-*','MarkerSize',4,'color',[0,0,0]+0.5);
hold on
plot(share_esti_pepsi(3,:),'-*','MarkerSize',4,'color',[0,0,0]+0.5);
hold on
% L2 market share
L2 = plot(share_esti_pepsi_L2,'-*','MarkerSize',4,'color','r');
% plot grey region
hold on
patch([t fliplr(t)], [share_esti_pepsi(1,:) fliplr(share_esti_pepsi(2,:))], 'k', 'FaceAlpha',0.2)
hold on
patch([t fliplr(t)], [share_esti_pepsi(2,:) fliplr(share_esti_pepsi(2,:))], 'k', 'FaceAlpha',0.2)
hold on
% OLS market share
ols = plot(share_esti_ols_pepsi,'-s','Color','k','MarkerSize',5,'MarkerFaceColor','k');
ylim([0,1])
xticks(1:4:18)
xticklabels({'1968','1972','1976','1980','1984'})
ax = gca;
ax.FontSize = 20;
xlabel('Years','FontSize',18,'FontWeight','bold')
ylabel('Pepsis beliefes on market share evolution','FontSize',18,'FontWeight','bold')
%legend([true pepsi L2 ols],{'True market share','M-BRD: Pepsis perspective','M-BR-$L_2$-D: Pepsis perspective','OLS'}, 'interpreter','latex','FontSize',9);
% pbaspect([2 1 1])
grid on

matlab2tikz('./Figures/market_share_evolution_error_5percent.tikz');
%%
% compute error assuming the new advertising expenditure of the competitior
% is known

% advertising expenditure with diminishing returns (square
% root)
% lambda = 0,1,2 --> not lagged, once lagged, twice lagged
function [X_constr, Y_constr,X_rows] = compute_error_lanchester(i,bound,adco,adpe,mco,lambda)
X = [];
Y = [];
if i == 1
    for t=2:length(adco)-1
        % build X matrix
        k_i = sqrt(adco(t,1))*(1-mco(t-1,1))-sqrt(adco(t+1,1))*(1-mco(t,1));
        k_minus_i = -(sqrt(adpe(t,1))*mco(t-1,1)-sqrt(adpe(t+1,1))*mco(t,1));
        k_lambda_i = sqrt(adco(t-1,1))*(1-mco(t-1,1))-sqrt(adco(t,1))*(1-mco(t,1));
        k_lambda_minus_i = -(sqrt(adpe(t-1,1))*mco(t-1,1)-sqrt(adpe(t,1))*mco(t,1));
        X = [X;
             k_i k_minus_i];
        % build Y_matrix
        const = mco(t-1,1)-mco(t,1);
        Y = [Y;
            -const];
    end
end

if i == 2
    for t=2:length(adco)-1
        % build X matrix
        k_i = sqrt(adpe(t,1))*mco(t-1,1)-sqrt(adpe(t+1,1))*mco(t,1);
        k_minus_i = -(sqrt(adco(t,1))*(1-mco(t-1,1))-sqrt(adco(t+1,1))*(1-mco(t,1)));
        k_lambda_i = sqrt(adpe(t-1,1))*mco(t-1,1)-sqrt(adpe(t,1))*mco(t,1);
        k_lambda_minus_i = -(sqrt(adco(t-1,1))*(1-mco(t-1,1))-sqrt(adco(t,1))*(1-mco(t,1)));
        X = [X;
            k_i k_minus_i];
        % build Y_matrix
        const = mco(t,1)-mco(t-1,1);
        Y = [Y;
             -const];
    end
end
    % constrain parameter space Theta=\mathbb(R)^2 to [-bound,bound]^2
    X_rows = size(X,1);
    X_constr = [];
    Y_constr = [];
        X_constr =[X;
              1 0;
              -1 0;
              0 1;
              0 -1];
        Y_constr = [Y;
               bound;
               bound;
               bound;
               bound];
end


% Chintagunta/Sorger: include word-of-mouth effects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X_constr, Y_constr,X_rows] = compute_error_sorger(i,bound,adco,adpe,mco,lambda)
X = [];
Y = [];
if i == 1
    for t=2:length(adco)-1
        % build X matrix
        k_i = sqrt(adco(t,1))*sqrt(1-mco(t-1,1))-sqrt(adco(t+1,1))*sqrt(1-mco(t,1));
        k_minus_i = -(sqrt(adpe(t,1))*sqrt(mco(t-1,1))-sqrt(adpe(t+1,1))*sqrt(mco(t,1)));
        if lambda == 0
            X = [X;
                k_i k_minus_i];
        end
         % build Y_matrix
         const = mco(t-1,1)-mco(t,1);
         Y = [Y;
             -const];
    end
end

if i == 2
    for t=2:length(adco)-1
        % build X matrix
        k_i = sqrt(adpe(t,1))*sqrt(mco(t-1,1))-sqrt(adpe(t+1,1))*sqrt(mco(t,1));
        k_minus_i = -(sqrt(adco(t,1))*sqrt(1-mco(t-1,1))-sqrt(adco(t+1,1))*sqrt(1-mco(t,1)));
        if lambda == 0
            X = [X;
                k_i k_minus_i];
        end
        % build Y_matrix
        const = mco(t,1)-mco(t-1,1);
        Y = [Y;
            -const];
    end
end
    % constrain parameter space Theta=\mathbb(R)^2 to [-bound,bound]^2
    X_rows = size(X,1);
    X_constr = [];
    Y_constr = [];
    if lambda == 0
        X_constr =[X;
              1 0;
              -1 0;
              0 1;
              0 -1];
        Y_constr = [Y;
               bound;
               bound;
               bound;
               bound];
    end
end

function [x] = solve_LP(X_constr,Y_constr,X_rows)

% Augment X to rephrase as LP by adding a column for eps
% (we have inequalities of the form e_i(t+1,P_i)<=eps)
eps_vec = [-ones(X_rows,1);
           zeros(size(X_constr,1)-X_rows,1)];

A = [X_constr, eps_vec];
% add constraint eps>=0
A = [A; 
    zeros(1,size(X_constr,2)), -1];
b = [Y_constr;
    0];

% Objective function 
% which is of the form
% f^T*x=0*alpha_ii+0*aplpha_i-i+0*gamma_ii+0*gamma_i-i+1+eps
f = [zeros(1,size(X_constr,2)) 1];

% Solve LP
x = linprog(f,A,b);
end


