%% MMF1921 (Summer 2020) - Project 1
% 
% The purpose of this program is to implement the following factor models
% a) Multi-factor OLS regression
% b) Fama-French 3-factor model
% c) LASSO
% d) Best Subset Selection
% 
% and to use these factor models to estimate the asset expected returns and 
% covariance matrix. These parameters will then be used to test the 
% out-of-sample performance using MVO to construct optimal portfolios.
% 
% Use can use this template to write your program.
%
% Student Name: Pulkit Gaur
% Student ID: 1005941947

clc
clear all
format short

% Program Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Read input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the stock weekly prices
adjClose = readtable('MMF1921_AssetPrices.csv');
adjClose.Properties.RowNames = cellstr(datetime(adjClose.Date));
adjClose.Properties.RowNames = cellstr(datetime(adjClose.Properties.RowNames));
adjClose.Date = [];

% Load the factors weekly returns
factorRet = readtable('MMF1921_FactorReturns.csv');
factorRet.Properties.RowNames = cellstr(datetime(factorRet.Date));
factorRet.Properties.RowNames = cellstr(datetime(factorRet.Properties.RowNames));
factorRet.Date = [];

riskFree = factorRet(:,9);
factorRet = factorRet(:,1:8);

% Identify the tickers and the dates 
tickers = adjClose.Properties.VariableNames';
dates   = datetime(factorRet.Properties.RowNames);

% Calculate the stocks' weekly EXCESS returns
prices  = table2array(adjClose);
returns = ( prices(2:end,:) - prices(1:end-1,:) ) ./ prices(1:end-1,:);
returns = returns - ( diag( table2array(riskFree) ) * ones( size(returns) ) );
returns = array2table(returns);
returns.Properties.VariableNames = tickers;
returns.Properties.RowNames = cellstr(datetime(factorRet.Properties.RowNames));

% Align the price table to the asset and factor returns tables by
% discarding the first observation.
adjClose = adjClose(2:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Define your initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial budget to invest ($100,000)
initialVal = 100000;

% Start of in-sample calibration period 
calStart = datetime('2008-01-01');
calEnd   = calStart + calyears(4) - days(1);

% Start of out-of-sample test period 
testStart = datetime('2012-01-01');
testEnd   = testStart + calyears(1) - days(1);

% Number of investment periods (each investment period is 1 year long)
NoPeriods = 5;

% Factor models
% Note: You must populate the functios OLS.m, FF.m, LASSO.m and BSS.m with your
% own code.
FMList = {'OLS' 'FF' 'LASSO' 'BSS' 'EqualWeight'};

FMList = cellfun(@str2func, FMList, 'UniformOutput', false);
NoModels = length(FMList);

% Tags for the portfolios under the different factor models
tags = {'OLS portfolio' 'FF portfolio' 'LASSO portfolio' 'BSS portfolio' 'EqualWeight portfolio'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Construct and rebalance your portfolios
%
% Here you will estimate your input parameters (exp. returns and cov. 
% matrix etc) from the Fama-French factor models. You will have to 
% re-estimate your parameters at the start of each rebalance period, and 
% then re-optimize and rebalance your portfolios accordingly. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiate counter for the number of observations per investment period
toDay = 0;

% Preallocate the space for the per period value of the portfolios 
currentVal = zeros(NoPeriods, NoModels);

%--------------------------------------------------------------------------
% Set the value of lambda and K for the LASSO and BSS models, respectively
%--------------------------------------------------------------------------
% lambda = we dont need to choose a lambda
% K      = 

for t = 1 : NoPeriods
  
    % Subset the returns and factor returns corresponding to the current
    % calibration period.
    periodReturns = table2array( returns( calStart <= dates & dates <= calEnd, :) );
    periodFactRet = table2array( factorRet( calStart <= dates & dates <= calEnd, :) );
    currentPrices = table2array( adjClose( ( calEnd - days(7) ) <= dates ... 
                                                    & dates <= calEnd, :) )';
    
    % Subset the prices corresponding to the current out-of-sample test 
    % period.
    periodPrices = table2array( adjClose( testStart <= dates & dates <= testEnd,:) );
    
    % Set the initial value of the portfolio or update the portfolio value
    if t == 1
        
        currentVal(t,:) = initialVal;
        
    else
        for i = 1 : NoModels
            
            currentVal(t,i) = currentPrices' * NoShares{i};
            
        end
    end
    
    % Update counter for the number of observations per investment period
    fromDay = toDay + 1;
    toDay   = toDay + size(periodPrices,1);
    
    % Calculate 'mu' and 'Q' using the 4 factor models.
    % Note: You need to write the code for the 4 factor model functions. 
    lambda = [];
    K= 4;
    for i = 1 : NoModels
        
        [mu{i}, Q{i}, R_squre_adj{i}] = FMList{i}(periodReturns, periodFactRet, lambda, K);
        
    end
            
    % Optimize your portfolios to get the weights 'x'
    % Note: You need to write the code for MVO with no short sales
    for i = 1 : NoModels
    % Define the target return as the geometric mean of the market 
    % factor for the current calibration period
    targetRet = geomean(periodFactRet(:,1) + 1) - 1;
    if (i==NoModels)
    x{i}(:,t) = ones(20,1)* (1/20);
    continue
    end
    x{i}(:,t) = MVO(mu{i}, Q{i}, targetRet); 
    end
    
    % Calculate the optimal number of shares of each stock you should hold
    for i = 1 : NoModels
        
        % Number of shares your portfolio holds per stock
        NoShares{i} = x{i}(:,t) .* currentVal(t,i) ./ currentPrices;
        
        % Weekly portfolio value during the out-of-sample window
        portfValue(fromDay:toDay,i) = periodPrices * NoShares{i};
        
    end

    % Update your calibration and out-of-sample test periods
    R_squre_list{t}=R_squre_adj;
    calStart = calStart + calyears(1);
    calEnd   = calStart + calyears(4) - days(1);
    
    testStart = testStart + calyears(1);
    testEnd   = testStart + calyears(1) - days(1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% 4.1 Evaluate any measures of fit of the regression models to assess their
% in-sample quality. You may want to modify Section 3 of this program to
% calculate the quality of fit each time the models are recalibrated.
%--------------------------------------------------------------------------
% Average adjusted R_squre
for t = 1:NoPeriods
    R_bar_different_model_mean{t} = [mean(R_squre_list{t}{1}) mean(R_squre_list{t}{2})...
                                mean(R_squre_list{t}{3}) mean(R_squre_list{t}{4})];
end

for t = 1:NoPeriods
    R_bar_different_model_mean{t} = [std(R_squre_list{t}{1}) std(R_squre_list{t}{2})...
                                std(R_squre_list{t}{3}) std(R_squre_list{t}{4})];
end

for t = 1:NoPeriods
    R_different_model{t} = [R_squre_list{t}{1}; R_squre_list{t}{2};...
                                R_squre_list{t}{3}; R_squre_list{t}{4}]';
end
%--------------------------------------------------------------------------
% 4.2 Calculate the portfolio average return, variance (or standard 
% deviation), and any other performance and/or risk metric you wish to 
% include in your report.
%--------------------------------------------------------------------------
portfValue = [100000*ones(1,size(portfValue,2));portfValue];
% Portfolio Returns
portfRet = (portfValue(2:end,:)-portfValue(1:end-1,:))./portfValue(1:end-1,:);

% Portfolio CAGR
portfCAGR = ((portfValue(end,:)./portfValue(1,:)).^(12/size(portfValue,1))-1);

% Portfolio Percentage Return
portfPercentageReturn = (portfValue(end,:)-portfValue(1,:))./portfValue(1,:);

% Portfolio average return
portfRet_bar = mean(portfRet,1);

% Portfolio standard deviation
portfRet_std = std(portfRet,1);

% Portfolio Sharp ratio
portf_sharpRatio = (portfRet_bar-mean(table2array(riskFree)))./portfRet_std;

% % Portfolio Maximum Drawdown
% portf_MaxDrawdown = maxdrawdown(portfValue);
% 
% % Portfolio skewness
% portf_skewness = skewness(portfRet);
% 
% % Portfolio Value at Risk (90% Threshold)
% RiskThreshold = 0.10;
% ValueAtRisk_90 = -portvrisk(portfRet_bar,portfRet_std, RiskThreshold);
% 
% % Portfolio Value at Risk (95% Threshold)
% RiskThreshold = 0.05;
% ValueAtRisk_95 = -portvrisk(portfRet_bar,portfRet_std, RiskThreshold);
%--------------------------------------------------------------------------
% 4.3 Plot the portfolio wealth evolution 
% 
% Note: The code below plots all portfolios onto a single plot. However,
% you may want to split this into multiple plots for clarity, or to
% compare a subset of the portfolios. 
%--------------------------------------------------------------------------
plotDates = dates(dates >= datetime('2012-01-01') );

fig1 = figure(1);

for i = 1 : NoModels
    
    plot( plotDates, portfValue(2:end,i) )
    hold on
    
end

legend(tags, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Portfolio value', 'FontSize', 14)
ylabel('Value','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig1,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
print(fig1,'fileName','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
% print(fig1,'fileName','-dpng','-r0');

%--------------------------------------------------------------------------
% 4.4 Plot the portfolio weights period-over-period
%--------------------------------------------------------------------------

% OLS portfolio weights

for i = 1 : NoModels
fig2 = figure(i+1);
nameofModel=string(tags{i});
area(x{i}','FaceColor','flat')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title(nameofModel + ' weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance period','interpreter','latex','FontSize',12);
% Define the plot size in inches
set(fig2,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig2,'Position');
set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches',...
'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
print(fig2,'weights '+ nameofModel ,'-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
% print(fig2,'weights '+ nameofModel ,'-dpng','-r0');
end
%% plot cumulative returns for each assets seperately
Returns = table2array( returns( datetime('2012-01-01') <= dates & dates <= datetime('2017-01-01'), :) );
Return_cum = cumprod(Returns+1)-1;
plot( plotDates, Return_cum )
asset_tags = {'F' 'CAT' 'DIS' 'MCD' 'KO' 'PEP' 'WMT' 'C' 'WFC' 'JPM'...
              'AAPL' 'IBM' 'PFE' 'JNJ' 'XOM' 'MRO' 'ED' 'T' 'VZ' 'NEM'};
legend(asset_tags, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Cumulative Returns for Assets', 'FontSize', 14)
ylabel('Cumulative Returns','interpreter','latex','FontSize',12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program End