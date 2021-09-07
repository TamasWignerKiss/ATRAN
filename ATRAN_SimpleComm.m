%% Description
% The second spimplest model
%
% Assumptions:
% - Tasks are non-specific. (They could suite the dominant function of the group!, which might be important, who knows...)
% - Agents all know each other:
%   - They know skills of all other agents
%   - Consequently, they know how similar (Euclidean distance in function space) they are to each other
%   - They know if other agents are free or busy
% - All agents use the same distance measure and threshold (why complicate things if not necessary)
%
% Simulation steps
% - Generate some agents with given DFD and IFD
% - Give some of them a task to work on
% - Keep doing the following until tasks are all solved or EmergenyStop
% - Tasked agents check if
%   - Other agents they are willing to talk to are better in solving their task than themselves
%   - a better agent is free, pass the task to it | if not, work on the task
%   
% I think, if multiple passing is allowed the question will be how quickly a task gets to the most optimal agent. It might be more interesting to
% allow for only 1 pass. I'll try that too.

%% Parameters
% Basic setup parameters
par.numfuncs = 9; %Number of domain functions
par.numagents = 100; %Number of agents in the system
par.numtasks = 10; %Number of tasks injected
par.agspread = 10; %A parameter specifying a meta-spread in agent similarity (could possibly be eliminated)
par.anorm = 10; %Total skill of agents
par.tnorm = 10; %Total labour requirement of tasks

% Internal simulation control
par.maxPass = 8; %Each task can be passed this many times
par.similThresh = 80; %The threshold of similarity (as percentage of maximal distance) above which agents are willing to communicate w/ each other

% External simulation control
par.numrepeats = 10;
par.EmergencyStop = 5e2;

%% Some initialization
adivvals = logspace(-1, 3, 20);
gdivvals = logspace(-1, 3, 10);
DFD = NaN(20,10);
IFD = NaN(20,10);
minsn = NaN(20, 10);
maxsn = NaN(20, 10);
meansn = NaN(20, 10);
etc.maxPass = par.maxPass;
etc.similThresh = par.similThresh;
etc.emStop = par.EmergencyStop;

%% Loop through DFD and IFD values
tic
T = 0;
ai = 1;
for adiv = adivvals
    gi = 1;
    fprintf('%i ', ai)
    for gdiv = gdivvals
        %Generate agents
        agents = GenAgent(par.numfuncs, par.numagents, [adiv, par.agspread], exp((-(0:par.numfuncs-1).^2)/gdiv), par.anorm);
        %Calculate and store diversity values
        [DFD(ai, gi), IFD(ai, gi)] = CalcFD(agents);
        
        %Run some repetitions of task solving with the same parameters
        sn = NaN(par.numrepeats, 1);
        for ridx = 1:par.numrepeats
            %Generate tasks
            tasks = GenTask(par.numfuncs, par.numtasks, par.tnorm);
            sn(ridx) = PassingSolveTasks(agents, tasks, etc); %This returns number of steps required to solve all tasks
        end
        minsn(ai, gi) = min(sn);
        maxsn(ai, gi) = max(sn);
        meansn(ai, gi) = mean(sn);
        gi = gi + 1;
        fprintf('.')
    end
    T_old = T;
    T = toc;
    fprintf(' (%0.2f sec)\n', T-T_old)
    ai = ai + 1;
end
toc

%% Plot what we got :-)
figure
surf(DFD, IFD, meansn)
title(['Performance Plot | Simple Task Passing w/ passing limit=' num2str(par.maxPass) ' | SimThresh=' num2str(par.similThresh) '%'])
xlabel('Dominant Functional Diversity')
ylabel('Average Individual Functional Diversity')
zlabel('Average Time of Completion [steps]')