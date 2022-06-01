function [t2a, taskhist] = PassingSolveTasks(agents, tasks, etc, agdistmat)
%This function runs an instance of the simulation. In this version the simulation runs until the maximal number of steps is reached or all tasks are
%solved. The following processes are performed:
%1. Randomly assign tasks to agents (also, keep track of the assignment across stimulation steps)
%2. Figure out is another free agent is better in solving a task and if so, pass the task
%3. Agents work on the task
%4. Go to step 2.
%
% Input arguments:
% agents: is the agents matrix
% tasks: is the task matrix
% etc: other parameters, a struct containing:
%   emStop: a number that specifies the maximum number of steps allowed (to stop simulation in case the system gets stuck)
%   maxPass: maximum number of allowable passes
%   similThresh: a threshold of similarity (percent of max) above which agents don't pass

%% Some initialization
stepNo = 0; %Counter to track number of simulation steps needed to solve all tasks
numPass = etc.maxPass*ones(size(tasks,1), 1); %To assign the same numpass to all tasks at the start of the simulation
maxDist = max(pdist(agents)); %The maximum distance among agents. Used to calculate absolute passing threshold from relative

%% Solving the tasks

%Assign tasks to agents
t2a = NaN(etc.emStop+1, size(tasks, 1));
t2a(1, :) = randperm(size(agents,1), size(tasks,1));

%This variable will store task history
taskhist = NaN(size(tasks, 1), size(tasks, 2), etc.emStop+1);
taskhist(:, :, 1) = tasks;

% Iterate until there are tasks or emergeny stop reached
while any(not(isnan(t2a(stepNo+1, :)))) && stepNo < etc.emStop
    
    %Tasked agents check other agents they are willing to talk to, and pass their task if other is better and free
    [t2a(stepNo+2, :), numPass] = passTasks_v2(t2a(stepNo + 1, :), agents, tasks, maxDist*etc.similThresh/100, numPass, agdistmat);
    %[t2a(stepNo+2, :), numPass] = passTasks(t2a(stepNo + 1, :), agents, tasks, maxDist*etc.similThresh/100, numPass);
    
    %Agents work on tasks, mark by NaN if a component is completed
    %Need to find busy agents and active tasks because NaN t2a cannot be used for index
    ActiveTaskIdx = not(all(isnan(tasks), 2));
    BusyAgentIdx = not(isnan(t2a(stepNo+2, :)));
    tasks(ActiveTaskIdx, :) = tasks(ActiveTaskIdx, :) - agents(t2a(stepNo+2, BusyAgentIdx), :);

    tasks(tasks <= 0) = NaN;
    
    %Find finished tasks
    isFinished = all(isnan(tasks), 2);

    %Free agent
    t2a(stepNo+2, isFinished) = NaN;

    %Store current tasks
    taskhist(:, :, stepNo + 2) = tasks;
    
    %Increase step counter
    stepNo = stepNo + 1;
end

% No work was done using the initial task assignment, only passing happened
t2a = t2a(2:end, :);
