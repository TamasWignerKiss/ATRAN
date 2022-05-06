function [DFD, IFD, IFDS] = CalcFD(agents)
%This function calculates DFD and IFD values (as well as IFDS, which is intrapersonal diversity score, the agent-level measure of their functional
%diversity). Calculations are based on Bunderson and Sutcliffe, 2002, Academy of Management Journals paper.
%
%Usage: [DFD, IFD, IFDS] = CalcFD(agents)

%% Calculate DFD
maxDFD = 1-1/size(agents, 2); %This is the maxmial value DFD can take.
[~, I] = max(agents, [], 2);
T = tabulate(I);
DFD = 1 - sum((T(:,3)/100).^2); %This is the first formula in the Bunderson & Sutcliffe, 2002 paper
DFD = DFD/maxDFD; %This is the normalization they mention on page 885 below the formula

%% Calculate IFD
% First calculate IFDS
IFDS = 1 - sum((agents./sum(agents,2)).^2, 2);
% Second, mean IFDSs across agents to get IFD
IFD = mean(IFDS);

end

