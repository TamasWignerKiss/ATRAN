%% Parameters

%Simulation
par.nofunc = 9; %The number of functions and skills in the simulation
par.noagent = 100; %Number of agents to be generated
par.anorm = 10; %The sum of capabilities of agents
par.aDivm = 100; %This is where intra-agent diversity can be set: the higher, the more diverse
par.aDivd = 10; %This introduces some spread into inta-agent diversity (so that IDA is not the same for all agents)
par.gDiv = ones(par.nofunc,1); %This sets group diversity. [1 0 0 0 0 0 0 0 0] would result in DFD == 0
par.notask = 10; %Generate this many tasks when task generation is required
par.tnorm = 100; %The sum of skill requirement of tasks

%Etc
par.debug = true;

%% Generate agents
% Agents have par.nofunc functions. This is stored in an array and the agents' capability in any given function is the value in the array at a given
% position. Values can range from 0 to 10, but their sum is normalized to 10. Capabilities are generated by a Gaussian function. It's width specifies
% how universal the agent is: the wider the Gaussian the agent has high values (good capabilities) at multiple functions, if narrower, only at few
% functions. This will set the intra-personal functional diversity value.
%
% For the dominant functional diversity the maximum of the array is to be varied across agents: if multiple functions have the same probability, DFD
% will be high; if only a few functions have high probability, DFD will be low.

x = repmat(0:(par.nofunc-1), par.noagent, 1); %Create "functions"
domf = randsample(x(1,:), par.noagent, true, par.gDiv); %Find a dominant function for agents
agents = exp((-(x).^2)./(par.aDivm+par.aDivm*par.aDivd/100*rand(par.noagent,1)-par.aDivd/200)); %Initialize agents
agents = agents./sum(agents,2)*par.anorm; %Normalize agents to have a total capability of par.anorm
%Now some final touches
for aidx = 1:par.noagent
    %Randomly mix up functions so not all agents have the same strength for given functions
    agents(aidx,:) = agents(aidx, randperm(length(agents(aidx,:))));
    %And finally, put the strongest function to its place determined by domf
    swp = agents(aidx,domf(aidx)+1);
    [m, l] = max(agents(aidx,:));
    agents(aidx,domf(aidx)+1) = m;
    agents(aidx, l) = swp;
end
clear x aidx swp m l

%% Plot agents for checking
if par.debug
    numagents = size(agents, 1);
    DF = NaN(numagents,1);
    colors = jet(numagents);
    figure
    hold on
    for aidx = 1:numagents
        plot(agents(aidx, :), '*', 'Color', colors(aidx,:));
        DF(aidx) = find(agents(aidx,:) == max(agents(aidx, :))); %This finds the dominant function
        %     fprintf('%i\n', DF(aidx))
    end
    tabulate(DF)
    clear colors aidx DF
end

%% Calculate DFD
% For dominant functional diversity domf can be used.
maxDFD = 1-1/par.nofunc; %This is the maxmial value DFD can take.
T = tabulate(domf); %This is to see distribution of agents' DF across all functions
DFD = 1 - sum((T(:,3)/100).^2); %This is the first formula in the Bunderson & Sutcliffe, 2002 paper
DFD = DFD/maxDFD; %This is the normalization they mention on page 885 below the formula
fprintf('Dominant functional diversity of this group of agents is: %0.4f\n', DFD)
clear maxDFD T

%% Calculate IFD
% First calculate IFDS
IFDS = 1 - sum((agents./sum(agents,2)).^2, 2);
IFD = mean(IFDS);
fprintf('Group-average intrapersonal dominant functional diversity of this group of agents is: %0.4f\n', IFD)
clear IFDS

%% Generate tasks
tasks = rand(par.notask, par.nofunc);
tasks = tasks./sum(tasks, 2)*par.tnorm;

%% Assign tasks to agents
%The array below serves to connect tasks to agents. The index of the array element specifies the task's ID and the value at that element specifies the
%number of the agent
t2a = randperm(par.noagent, size(tasks,1));

%% This is a way to implement solving tasks
stepno = 1;
while any(any(tasks > 0))
    tasks = tasks - agents(t2a,:);
    if par.debug
        fprintf('Step #%i: %i tasks remaining\n', stepno, sum(any(tasks > 0,2)))
    end
    stepno = stepno + 1;
end
stepno = stepno -1;
fprintf('%i tasks were solved in %i steps.\n\n', size(tasks, 1), stepno)
