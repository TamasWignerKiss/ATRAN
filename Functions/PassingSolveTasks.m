function [stepNo, t2a, taskhist] = PassingSolveTasks(agents, tasks, etc)
% Input arguments:
% agents: is the agents matrix
% tasks: is the task matrix
% etc: other parameters
%  emStop: a number that specifies the maximum number of steps allowed (to stop simulation in case the system gets stuck)
%  maxPass: maximum number of allowable passes
%  similThresh: a threshold of similarity (percent of max) above which agents don't pass

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
while all(not(isnan(t2a(stepNo+1, :)))) && stepNo < etc.emStop
    
    %Tasked agents check other agents they are willing to talk to, and pass their task if other is better and free
    [t2a(stepNo+2, :), numPass] = passTasks(t2a(stepNo + 1, :), agents, tasks, maxDist*etc.similThresh/100, numPass);
    
    %Agents work on tasks, mark by NaN if a component is completed
    tasks = tasks - agents(t2a(stepNo+2, :),:);
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