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

% Simulation parameters
par.numfuncs = 9; %Number of domain functions
par.numagents = 10; %Number of agents in the system
par.numtasks = 7; %Number of tasks injected
par.agspread = 10; %A parameter specifying a meta-spread in agent similarity (could possibly be eliminated)
par.anorm = 10; %Total skill of agents (agent normalizing factor)
par.tnorm = 10; %Total labour requirement of tasks (task normalizing factor)
par.TaskType = 2; %Sets what type of task to generate. See help GenTask
par.noavals = 40; %Number of agent diversity steps
par.nogvals = 40; %Number of group diversity steps

% Internal simulation control parameters
par.maxPass = Inf; %Each task can be passed this many times
par.similThresh = 80; %The threshold of similarity (as percentage of maximal distance) below which agents are willing to communicate w/ each other
par.passCost = 100; % Steps, the cost of communication: the receiving agent cannot work until this amount of time stamps

% External simulation control parameters
par.numrepeats = 10;
par.EmergencyStop = 2.5e2;

% Plotting
par.IFDbins = linspace(0, 1, 11); %For coarse-graining the values calculated
par.DFDbins = linspace(0, 1, 11);

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
etc.passCost = par.passCost; %Steps

%% Loop through DFD and IFD values
tic
T = 0;
for ai = 1:length(adivvals)
    adiv = adivvals(ai);
    fprintf('%i ', ai)
    
    parfor gi = 1:length(gdivvals) %This is for parallel computing
        gdiv = gdivvals(gi);

        %Generate agents
        agents = GenAgent(par.numfuncs, par.numagents, [adiv, par.agspread], exp((-(0:par.numfuncs-1).^2)/gdiv), par.anorm); %#ok<PFBNS> 
        agdistmat = squareform(pdist(agents));

        %Calculate and store diversity values
        [DFD(ai, gi), IFD(ai, gi)] = CalcFD(agents);
        
        %The following variables will only live for an [ai, gi] point
        t2a = NaN(par.EmergencyStop+1, par.numtasks, par.numrepeats); %Initialize variable keeping track of which agent works on which task
        taskhists = NaN(par.numtasks, par.numfuncs, par.EmergencyStop+1, par.numrepeats); %Initialize variable keeping track of how tasks progress
        sn = NaN(par.numrepeats, 1); %Initialize number of steps nedded to solve all tasks
        stn = zeros(par.numrepeats, 1); %Initialize variable storing the amount of work done in repeats
        ca = NaN(par.numrepeats, 1); %Init variable for communication amount

        %The following piece of code runs in parallel on multiple cores. Remove "par" to allow debugging of the inside.
        for ridx = 1:par.numrepeats
            
            %Generate tasks
            tasks = GenTask(par.numfuncs, par.numtasks, par.tnorm, par.EmergencyStop, par.TaskType, agents);

            %PassingSolveTasks is the function that runs the while loop running the whole game of passing-working
            [t2a(:, :, ridx), taskhists(:, :, :, ridx)] = PassingSolveTasks(agents, tasks, etc, agdistmat); %This returns number of steps required to solve all tasks

            % Calculate steps taken to solve all tasks
            tmp = find(all(isnan(t2a(:, :, ridx)), 2), 1, 'first')';
            if isempty(tmp)
                sn(ridx) = par.EmergencyStop;
            else
                sn(ridx) = tmp;
            end

            % Calculate the amount of solved sub-tasks
            tmp1 = taskhists(:, :, 1, ridx);
            tmp2 = taskhists(:, :, par.EmergencyStop+1, ridx);
            tmp1(isnan(tmp1)) = 0; %There can be NaNs initially due to sub-tasks that the system cannot solve!
            tmp2(isnan(tmp2)) = 0;
            stn(ridx) = sum(sum(tmp1 - tmp2));

            % Calculate amount of communication
            tmp = diff(t2a(:, :, ridx));
            ca(ridx) = sum(sum(tmp ~= 0 & not(isnan(tmp))));

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

        % Average communication amount
        minca(ai, gi) = min(ca);
        maxca(ai, gi) = max(ca);
        meanca(ai, gi) = mean(ca);
        stdca(ai, gi) = std(ca);

        fprintf('.')
    end

    % Step counter
    T_old = T;
    T = toc;
    fprintf(' (%0.2f sec)\n', T-T_old)
end
toc

% Interim cleanup
clear adiv agents ai gdiv gi sn stn T t2a T_old taskhists

%% Coarse-grain (average in bins) results and normalize max of amount of solved subtasks to 1
cgmstn = NaN(length(par.IFDbins)-1, length(par.DFDbins)-1);
cgmsn = NaN(size(cgmstn));
cgca = NaN(size(cgmstn));
cgIFD = NaN(size(cgmstn));
cgDFD = NaN(size(cgmstn));
for iidx = 1:length(par.IFDbins)-1
    for didx = 1:length(par.DFDbins)-1
        cgmstn(iidx, didx) = mean(meanstn(par.IFDbins(iidx)<=IFD & IFD<par.IFDbins(iidx+1) & par.DFDbins(didx)<=DFD & DFD<par.DFDbins(didx+1)));
        cgmsn(iidx, didx)  = mean(meansn(par.IFDbins(iidx)<=IFD & IFD<par.IFDbins(iidx+1) & par.DFDbins(didx)<=DFD & DFD<par.DFDbins(didx+1)));
        cgca(iidx, didx)   = mean(meanca(par.IFDbins(iidx)<=IFD & IFD<par.IFDbins(iidx+1) & par.DFDbins(didx)<=DFD & DFD<par.DFDbins(didx+1)));
        cgIFD(iidx, didx)  = par.IFDbins(iidx);
        cgDFD(iidx, didx)  = par.DFDbins(didx);
    end
end

cgmstn = cgmstn/(par.numtasks*par.tnorm);

% Interim cleanup
clear iidx didx

%% Plot amount of work solved vs. IFD, DFD (coarse-grained), also amount of communication
figure
p = panel();
p.pack(1,2)

% The work plot
p(1,1).select()
ax1 = gca;
surf(cgDFD, cgIFD, cgmstn)
grid on
title(['Passing cost=' num2str(par.passCost) ' | SimThresh=' num2str(par.similThresh) '% | NoAgs= ' ...
    num2str(par.numagents) ' | NoTasks= ' num2str(par.numtasks) ' | NoRepeats= ' num2str(par.numrepeats)])
xlabel('Dominant Functional Diversity')
ylabel('Individual Functional Diversity')
zlabel('Average Solved Subtasks')
% cbh = colorbar;
% ylabel(cbh, 'Average Solved Subtasks')
xlim([0 1])
ylim([0 1])
zlim([0 1])
clim([0 1])

% The communication plot
p(1,2).select()
ax2 = gca;
surf(cgDFD, cgIFD, cgca)
grid on
xlabel('Dominant Functional Diversity')
ylabel('Individual Functional Diversity')
zlabel('Amount of Communication')
% cbh = colorbar;
% ylabel(cbh, 'Amount of Communication')
xlim([0 1])
ylim([0 1])

%Linking the two plots
Link = linkprop([ax1, ax2],{'View'});
setappdata(gcf, 'StoreTheLink', Link);

%% Cuts

%Due to the way agents are generated, DFD and IFD values are not sorted. Sorting them here.
% [tmpDFD, sI] = sort(DFD, 2);
% tmpIFD = NaN(size(IFD));
% tmpMsn = NaN(size(meansn));
% tmpSsn = NaN(size(stdsn));
% for idx = 1:size(DFD, 1)
%     tmpIFD(idx,:) = IFD(idx, sI(idx, :));
%     tmpMsn(idx,:) = meansn(idx, sI(idx, :));
%     tmpSsn(idx,:) = stdsn(idx, sI(idx, :));
% end

%Time required to complete all tasks
%Mean
% figure
% %surf(tmpDFD, tmpIFD, meansn)
% surf(cgDFD, cgIFD, cgmsn)
% title(['Simple Task Passing w/ passing limit=' num2str(par.maxPass) ' | SimThresh=' num2str(par.similThresh) '% | NoAgs= ' ...
%     num2str(par.numagents) ' | NoTasks= ' num2str(par.numtasks) ' | NoRepeats= ' num2str(par.numrepeats)])
% xlabel('Dominant Functional Diversity')
% ylabel('Individual Functional Diversity')
% zlabel('Average Time of Completion [steps]')
%Std
% figure
% surf(tmpDFD, tmpIFD, stdsn)
% title(['Simple Task Passing w/ passing limit=' num2str(par.maxPass) ' | SimThresh=' num2str(par.similThresh) '% | NoAgs= ' ...
%     num2str(par.numagents) ' | NoTasks= ' num2str(par.numtasks) ' | NoRepeats= ' num2str(par.numrepeats)])
% xlabel('Dominant Functional Diversity')
% ylabel('Individual Functional Diversity')
% zlabel('Std of Time of Completion [steps]')

%Amount of work done
%Mean
% figure
% %surf(tmpDFD, tmpIFD, meanstn)
% surf(cgDFD, cgIFD, cgmstn)
% title(['Simple Task Passing w/ passing limit=' num2str(par.maxPass) ' | SimThresh=' num2str(par.similThresh) '% | NoAgs= ' ...
%     num2str(par.numagents) ' | NoTasks= ' num2str(par.numtasks) ' | NoRepeats= ' num2str(par.numrepeats)])
% xlabel('Dominant Functional Diversity')
% ylabel('Individual Functional Diversity')
% zlabel('Average Solved Subtasks')
%Std
% figure
% surf(tmpDFD, tmpIFD, stdstn)
% title(['Simple Task Passing w/ passing limit=' num2str(par.maxPass) ' | SimThresh=' num2str(par.similThresh) '% | NoAgs= ' ...
%     num2str(par.numagents) ' | NoTasks= ' num2str(par.numtasks) ' | NoRepeats= ' num2str(par.numrepeats)])
% xlabel('Dominant Functional Diversity')
% ylabel('Individual Functional Diversity')
% zlabel('Std of Solved Subtasks')

% Interim cleanup
% clear idx tmp*