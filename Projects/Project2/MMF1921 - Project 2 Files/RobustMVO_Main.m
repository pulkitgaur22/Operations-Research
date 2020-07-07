% MMF1921H - Robust MVO Matlab Example
% 
% Course Instructor: Professor Roy H. Kwon
%
% Prepared by:  Giorgio Costa 
% Email:        gcosta@mie.utoronto.ca
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear the memory and console
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
format short

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. DEFINE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          

% Number of assets
n = 3;

% Covariance matrix
Q  = [0.05 0.04 0.03; 0.04 0.06 0.03; 0.03 0.03 0.04];

% Estimates of the expected return
mu{1} = [0.0715; 0.0716; 0.06];            
mu{2} = [0.0716; 0.0715; 0.06];            

% Risk aversion coefficient
lambda = 0.02;

% Confidence level
alpha = 0.9;

% Number of return observations 
T = 20;

% Allocate space for our solutions
x  = NaN(n,4);

for i = 1:2
    %----------------------------------------------------------------------
    % 1. MVO
    %----------------------------------------------------------------------
    x(:,i) = MVO(mu{i}, Q, lambda);

    %----------------------------------------------------------------------
    % 2. Robust MVO
    %----------------------------------------------------------------------
    x(:,i+2) = robustMVO(mu{i}, Q, lambda, alpha, T);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Plot Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig1 = figure(1);
bar(x','stacked');
set(gca,'TickLabelInterpreter', 'latex','fontsize',22);
set(gca, 'XTickLabel', {'Nom. A','Nom. B','Rob. A', 'Rob. B'},'fontsize',22);
ylabel('Weights','interpreter', 'latex','FontSize',24);
title('Portfolio Weights','interpreter', 'latex','FontSize',24);


set(fig1,'Units','Inches', 'Position', [0 0 10, 10]);
    pos3 = get(fig1,'Position');
    set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pos3(3), pos3(4)])
    print(fig1,'Robust_Weights','-dpdf','-r0')


