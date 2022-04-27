function stepNo = PassingSolveTasks(agents, tasks, etc)
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
t2a = randperm(size(agents,1), size(tasks,1));

% Iterate until there are tasks or emergeny stop reached
while size(tasks,1) > 0 && stepNo < etc.emStop
        
    %Tasked agents check other agents they are willing to talk to, and pass their task if other is better and free
    [t2a, numPass] = passTasks(t2a, agents, tasks, maxDist*etc.similThresh/100, numPass);
    
    %Agents work on tasks
    tasks = tasks - agents(t2a,:);
    
    %Find finished tasks
    isFinished = all(tasks <= 0, 2);
    %Erase task
    tasks(isFinished, :) = [];
    %Free agent
    t2a(isFinished) = [];

    
    %Increase step counter
    stepNo = stepNo + 1;
end
