%% MMF1921 (Summer 2020) - Project 2
% 
% The purpose of this program is to provide a template with which to
% develop Project 2. The project requires you to test different models 
% (and/or different model combinations) to create an asset management
% algorithm. 

% Use can use this template to write your program. You may modify this
% template as you see fit, EXCEPT for Section 1. Section 1 must remain
% intact so that different data sets may be tested during the competition.
%
% Student Name:
% Student ID:

clc
clear all
format short

% Program Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Read input files 
% Do not modify this section of the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input file names
assetData  = 'MMF1921_AssetPrices2.csv';
factorData = 'MMF1921_FactorReturns2.csv';

% Initial budget to invest ($100,000)
initialVal = 100000;

% Length of investment period (in months, either 6 or 12)
investPeriod = 6;

% Load the stock weekly prices
adjClose = readtable(assetData);
adjClose.Properties.RowNames = cellstr(datetime(adjClose.Date));
adjClose.Properties.RowNames = cellstr(datetime(adjClose.Properties.RowNames));
adjClose.Date = [];

% Load the factors weekly returns
factorRet = readtable(factorData);
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

% Start of out-of-sample test period 
testStart = datetime(returns.Properties.RowNames{1}) + calyears(5);

% End of the first investment period
testEnd = testStart + calmonths(investPeriod) - days(1);

% Total number of investment periods
NoPeriods = ceil( days(datetime(returns.Properties.RowNames{end}) - ...
            testStart) / (30.44*investPeriod) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Construct and rebalance your portfolios
%
% Here you will estimate your input parameters (exp. returns and cov. 
% matrix etc) from the Fama-French factor models. You will have to 
% re-estimate your parameters at the start of each rebalance period, and 
% then re-optimize and rebalance your portfolios accordingly. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the start and end of in-sample calibration period. The first 
% calibration window can be at most 5 years long (from Jan-XX until 
% Dec-XX), and must not overlap with the testStart date. 
% calStart = 
% calEnd   = 

% Define the models and the total number of models you wish to use/test
% NoModels = 

% Initiate counter for the number of observations per investment period
toDay = 0;

% Preallocate the space for the per period value of the portfolios 
currentVal = zeros(NoPeriods, NoModels);

% Preallocate the space for the per period turnover rate of the portfolios 
turnover = zeros(NoPeriods, NoModels);

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
            
            x0{i}(:,t) = (currentPrices .* NoShares{i}) ./ currentVal(t,i);
            
        end
    end
    
    % Update counter for the number of observations per investment period
    fromDay = toDay + 1;
    toDay   = toDay + size(periodPrices,1);
    
    % Calculate 'mu' and 'Q' using the 4 factor models.
    % Note: You need to write the code for the 4 factor model functions. 
    for i = 1 : NoModels
        
%         [mu{i}, Q{i}] = 
        
    end
            
    % Optimize your portfolios to get the weights 'x'
    % Note: You need to write the code for MVO with no short sales
    for i = 1 : NoModels
        
%         x{i}(:,t) = 

        if t > 1
            turnover(t,i) = sum( abs( x{i}(:,t) - x0{i}(:,t) ) );
        end
            
    end
    
    % Calculate the optimal number of shares of each stock you should hold
    for i = 1 : NoModels
        
        % Number of shares your portfolio holds per stock
        NoShares{i} = x{i}(:,t) .* currentVal(t,i) ./ currentPrices;
        
        % Weekly portfolio value during the out-of-sample window
        portfValue(fromDay:toDay,i) = periodPrices * NoShares{i};
        
    end

    % Update your calibration and out-of-sample test periods
%     calStart = 
%     calEnd   = 
    
    testStart = testStart + calmonths(investPeriod);
    testEnd   = testStart + calmonths(investPeriod) - days(1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% 3.1 Evaluate any measures of fit of the regression models to assess their
% in-sample quality. You may want to modify Section 3 of this program to
% calculate the quality of fit each time the models are recalibrated.
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% 3.2 Calculate the portfolio average return, variance (or standard 
% deviation), and any other performance and/or risk metric you wish to 
% include in your report.
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% 3.3 Plot the portfolio wealth evolution 
% 
% Note: The code below plots all portfolios onto a single plot. However,
% you may want to split this into multiple plots for clarity, or to
% compare a subset of the portfolios. 
%--------------------------------------------------------------------------

plotDates = dates(dates >= datetime(returns.Properties.RowNames{1}) + calyears(5) );

fig1 = figure(1);

for i = 1 : NoModels
    
    plot( plotDates, portfValue(:,i) )
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
% print(fig1,'fileName','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig1,'fileName','-dpng','-r0');

%--------------------------------------------------------------------------
% 3.4 Plot the portfolio weights period-over-period
%--------------------------------------------------------------------------

% OLS portfolio weights
fig2 = figure(2);
area(x{1}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('OLS portfolio weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig2,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig2,'Position');
set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig2,'fileName2','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig2,'fileName2','-dpng','-r0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program End