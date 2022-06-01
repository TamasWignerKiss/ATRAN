function [t2a, numPass] = passTasks_v2(t2a, agents, tasks, simTh, numPass, agdistmat)
%This function performs reorganization of which agent works on which task. This is version 2 because in the old version weird things could happen. For
%example if Agent A was working on Task 1 and Agent B was working on Task 2 and both of them were stuck, they could not exchange tasks, since they
%were both showing up as busy for the other agent.
%
%In the current version all agents that communicate first evaluate which is better in solving a given task. This way for each task a list of agents
%are assigned. Second, tasks are assigned to agents in order of their fitness.
%
%Usage: [Task2Agent, NumberOfPasses] = pasTasks(Task2Agent, Agents, Tasks, SimilarityThreshold, NumberOfPasses)
%

%% Perform pass calculation

activeTasks = find(not(isnan(t2a)));
agFit  = cell(numel(t2a), 1);
simAgs = cell(numel(t2a), 1);

% First, for each active task evaluate how good agents are in solving them, taking into account communication
for tidx = activeTasks

    %Check if task was not passed too many times and allow further pass-evaluation only if passes are still allowed
    if numPass(tidx) > 0

        %Who is similar enough to be worthy for communication? -- asked the agent working on task tidx?
        simAgs{tidx} = find(agdistmat(t2a(tidx), :) < simTh);

        %Calculate how good all these friends are in solving task tidx
        agFit{tidx} = CalcSolvFitness(agents(simAgs{tidx}, :), tasks(tidx, :));
    end
end

%Second, convert simAgs and agFit to matrices for faster computation
maxAg = max(cellfun(@length, simAgs));
tmpSA = NaN(numel(t2a), maxAg);
tmpAF = NaN(numel(t2a), maxAg);
for tidx = 1:numel(t2a)
    if ~isempty(simAgs{tidx})
        tmpSA(tidx, 1:length(simAgs{tidx})) = simAgs{tidx};
        tmpAF(tidx, 1:length(simAgs{tidx})) = agFit{tidx};
    end
end
simAgs = tmpSA;
agFit  = tmpAF;

% I had quite some trouble with how to implement this resulting in programming errors :) Hence, some debug variables
simAgs_orig = simAgs;
agFit_orig = agFit;
t2a_orig = t2a;
activeTasks_orig = activeTasks;

% Third, assign agents to tasks, based on how much real work they can do. First assign the one that can do the most, ... until all tasks have an agent
try
    while ~isempty(activeTasks)
        [I, J] = find(min(min(agFit)) == agFit); %Which agent(J) at which task(I) would do the most work
        if isempty(I) %This happens when all agents available for the remaining tasks have been used up for previous tasks
            %In this case the original agents assigned to remaining tasks and all friend agents have been used up for other tasks (the agFit array is
            %empty.
            remTasks = activeTasks; %These are the tasks that cannot be assigned to agents
            activeTasks = [];
        else
            I = I(1); J = J(1); %In case of a draw, pick the first one
            bestAg = simAgs(I, J); %The #number of the best agent
            t2a(I) = bestAg; %(Re-)assign the task
            numPass(I) = numPass(I) - 1;
            tmp = simAgs == bestAg; %Locate all tasks where this agent shows up
            simAgs(tmp) = NaN; %Remove agents from the list of available agents -- this line is unnecessary
            agFit(tmp) = NaN; %Remove its fitnes
            agFit(I, :) = NaN; %Remove task from list of tasks-to-be-assigned
        end
        activeTasks(activeTasks == I) = []; %Remove task from list of tasks-to-be-assigned
    end
    
    % If there are remaining tasks, assign them to available agents
    if exist('remTasks', 'var')
        avAgs = setdiff(1:size(agents, 1), t2a(setdiff(1:end, remTasks))); %These are the available agents (agents not in t2a)
        t2a(remTasks) = randsample(avAgs, length(remTasks)); %Randomly assign available agents to remaining tasks
    end
catch ME
    fprintf('\nDump:\nI=%i, J=%i\n', I, J)
    disp(I)
    disp(J)
    fprintf('********************\n')
    disp(activeTasks_orig)
    disp(t2a_orig)
    disp(simAgs_orig)
    disp(agFit_orig)
    fprintf('********************\n')
    disp(activeTasks)
    disp(t2a)
    disp(simAgs)
    disp(agFit)
    rethrow ME
end

end
