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
% IO
par.functiondir = '/home/umat/bognor/ATRAN/Functions';

% Basic setup parameters
par.numfuncs = 9; %Number of domain functions
par.numagents = 10; %Number of agents in the system
par.numtasks = 5; %Number of tasks injected
par.agspread = 10; %A parameter specifying a meta-spread in agent similarity (could possibly be eliminated)
par.anorm = 10; %Total skill of agents
par.tnorm = 10; %Total labour requirement of tasks
par.TaskType = 2; %Sets what type of task to generate. See help GenTask
par.noavals = 20; %Number of agent diversity steps
par.nogvals = 10; %Number of group diversity steps

% Internal simulation control
par.maxPass = Inf; %Each task can be passed this many times
par.similThresh = 80; %The threshold of similarity (as percentage of maximal distance) above which agents are willing to communicate w/ each other

% External simulation control
par.numrepeats = 100;
par.EmergencyStop = 2.5e2;

% Etc
par.debug = 0;

%% Set the path
warning('off', 'MATLAB:rmpath:DirNotFound')
restoredefaultpath;
rmpath('/home/umat/Documents/MATLAB')
clear RESTOREDEFAULTPATH_EXECUTED
addpath(par.functiondir)
rehash
warning('on', 'MATLAB:rmpath:DirNotFound')

%% Some initialization
adivvals = 10.^(linspace(0.1, 0.2, par.noavals).*(linspace(-6, 7, par.noavals)));
gdivvals = 10.^(0.2*linspace(-3, 5, par.nogvals));
DFD = NaN(par.noavals, par.nogvals);
IFD = NaN(par.noavals, par.nogvals);
minsn = NaN(par.noavals, par.nogvals);
maxsn = NaN(par.noavals, par.nogvals);
meansn = NaN(par.noavals, par.nogvals);
stdsn = NaN(par.noavals, par.nogvals);
minstn = NaN(par.noavals, par.nogvals);
maxstn = NaN(par.noavals, par.nogvals);
meanstn = NaN(par.noavals, par.nogvals);
stdstn = NaN(par.noavals, par.nogvals);
etc.maxPass = par.maxPass;
etc.similThresh = par.similThresh;
etc.emStop = par.EmergencyStop;

%% Loop through DFD and IFD values
tic
T = 0;
ai = 1;
for adiv = adivvals
    gi = 1;
    if par.debug == 0
        fprintf('%i ', ai)
    end
    for gdiv = gdivvals
        
        %Generate agents
        agents = GenAgent(par.numfuncs, par.numagents, [adiv, par.agspread], exp((-(0:par.numfuncs-1).^2)/gdiv), par.anorm);
        if par.debug > 1
            PlotAgents(agents); %This is how you can visualize agents
        end

        %Calculate and store diversity values
        [DFD(ai, gi), IFD(ai, gi)] = CalcFD(agents);
        if par.debug > 0
            fprintf('DFD: %0.4f, IFD: %0.4f\n', DFD(ai, gi), IFD(ai, gi))
        end
        
        %The following variables will only live for an [ai, gi] point
        sn = NaN(par.numrepeats, 1); %Initialize number of steps nedded to solve all tasks
        t2a = NaN(par.EmergencyStop, par.numtasks, par.numrepeats); %Initialize variable keeping track of which agent works on which task
        taskhists = NaN(par.numtasks, par.numfuncs, par.EmergencyStop+1, par.numrepeats); %Initialize variable keeping track of how tasks progress
        stn = zeros(par.numrepeats, 1); %Initialize variable storing the amount of work done in repeats

        %The following piece of code runs in parallel on multiple cores. Remove "par" to allow debugging of the inside.
        parfor ridx = 1:par.numrepeats
            
            %Generate tasks
            tasks = GenTask(par.numfuncs, par.numtasks, par.tnorm, par.EmergencyStop, par.TaskType, agents); %#ok<PFBNS> 

            %PassingSolveTasks is the function that runs the while loop running the whole game of passing-working
            [t2a(:, :, ridx), tmpth] = PassingSolveTasks(agents, tasks, etc); %This returns number of steps required to solve all tasks
            taskhists(:, :, :, ridx) = tmpth; %This is technically needed for the parallel machine

            % Calculate steps taken to solve all tasks
            tmp = find(all(isnan(t2a(:, :, ridx)), 2), 1, 'first')';
            if isempty(tmp)
                sn(ridx) = par.EmergencyStop;
            else
                sn(ridx) = tmp;
            end

            % Calculate the amount of solved sub-tasks
            tmp1 = tmpth(:, :, 1);
            tmp2 = tmpth(:, :, 251);
            tmp1(isnan(tmp1)) = 0; %There can be NaNs initially due to sub-tasks that the system cannot solve!
            tmp2(isnan(tmp2)) = 0;
            stn(ridx) = sum(sum(tmp1 - tmp2));
        end

        % An example to plot how task assignment evolved in repeat repeatNum, and Task # taskNum
        if par.debug > 5
            repeatNum = 1;
            taskNum = 2;
            PlotRepeatEvents(repeatNum, taskNum, t2a, taskhists, par);
            clear repeatNum taskNum
        end
        
        % Average steps
        minsn(ai, gi) = min(sn);
        maxsn(ai, gi) = max(sn);
        meansn(ai, gi) = mean(sn);
        stdsn(ai, gi) = std(sn);

        % Average solved subtasks
        minstn(ai, gi) = min(stn);
        maxstn(ai, gi) = max(stn);
        meanstn(ai, gi) = mean(stn);
        stdstn(ai, gi) = std(stn);

        % Step counter
        gi = gi + 1;
        if par.debug == 0
            fprintf('.')
        end
    end

    % Step counter
    T_old = T;
    T = toc;
    fprintf(' (%0.2f sec)\n', T-T_old)
    ai = ai + 1;
end
toc

%% Plot what we've got :-)

%Due to the way agents are generated, DFD and IFD values are not sorted. Sorting them here.
[tmpDFD, sI] = sort(DFD, 2);
tmpIFD = NaN(size(IFD));
tmpMsn = NaN(size(meansn));
tmpSsn = NaN(size(stdsn));
for idx = 1:size(DFD, 1)
    tmpIFD(idx,:) = IFD(idx, sI(idx, :));
    tmpMsn(idx,:) = meansn(idx, sI(idx, :));
    tmpSsn(idx,:) = stdsn(idx, sI(idx, :));
end

%Time required to complete all tasks
%Mean
figure
surf(tmpDFD, tmpIFD, meansn)
title(['Simple Task Passing w/ passing limit=' num2str(par.maxPass) ' | SimThresh=' num2str(par.similThresh) '% | NoAgs= ' ...
    num2str(par.numagents) ' | NoTasks= ' num2str(par.numtasks) ' | NoRepeats= ' num2str(par.numrepeats)])
xlabel('Dominant Functional Diversity')
ylabel('Individual Functional Diversity')
zlabel('Average Time of Completion [steps]')
%Std
figure
surf(tmpDFD, tmpIFD, stdsn)
title(['Simple Task Passing w/ passing limit=' num2str(par.maxPass) ' | SimThresh=' num2str(par.similThresh) '% | NoAgs= ' ...
    num2str(par.numagents) ' | NoTasks= ' num2str(par.numtasks) ' | NoRepeats= ' num2str(par.numrepeats)])
xlabel('Dominant Functional Diversity')
ylabel('Individual Functional Diversity')
zlabel('Std of Time of Completion [steps]')

%Amount of work done
%Mean
figure
surf(tmpDFD, tmpIFD, meanstn)
title(['Simple Task Passing w/ passing limit=' num2str(par.maxPass) ' | SimThresh=' num2str(par.similThresh) '% | NoAgs= ' ...
    num2str(par.numagents) ' | NoTasks= ' num2str(par.numtasks) ' | NoRepeats= ' num2str(par.numrepeats)])
xlabel('Dominant Functional Diversity')
ylabel('Individual Functional Diversity')
zlabel('Average Solved Subtasks')
%Std
figure
surf(tmpDFD, tmpIFD, stdstn)
title(['Simple Task Passing w/ passing limit=' num2str(par.maxPass) ' | SimThresh=' num2str(par.similThresh) '% | NoAgs= ' ...
    num2str(par.numagents) ' | NoTasks= ' num2str(par.numtasks) ' | NoRepeats= ' num2str(par.numrepeats)])
xlabel('Dominant Functional Diversity')
ylabel('Individual Functional Diversity')
zlabel('Std of Solved Subtasks')
