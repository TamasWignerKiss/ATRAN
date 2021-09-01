%% Description
% The spimplest model:
% - Generate some agents with given DFD and IFD
% - Assign them some tasks
% - And see how much time they need to finish

%% Parameters
par.numfuncs = 9;
par.numagents = 100;
par.numtasks = 10;
par.agspread = 10;
par.anorm = 10;
par.tnorm = 10;
par.numrepeats = 10;
par.EmergencyStop = 5e2;

%% Check DFD and IFD values
% adivvals = logspace(-1, 3, 20);
% gdivvals = logspace(-1, 3, 20);
% DFD = NaN(20,20);
% IFD = NaN(20,20);
% ai = 1;
% for adiv = adivvals
%     gi = 1;    
%     for gdiv = gdivvals
%         agents = GenAgent(par.numfuncs, par.numagents, [adiv, par.agspread], exp((-(0:par.numfuncs-1).^2)/gdiv), par.anorm);
%         [DFD(ai, gi), IFD(ai, gi)] = CalcFD(agents);
%         gi = gi + 1;
%     end
%     ai = ai + 1;
% end
% 
% figure
% surf(gdivvals, adivvals, DFD)
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% axis tight
% title('Dominant Functional Diversity')
% xlabel('Group Diversity Parameter')
% ylabel('Individual Diversity Parameter')
% zlabel('DFD')
% 
% figure
% surf(gdivvals, adivvals, IFD)
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% axis tight
% title('Intrapersonal Functional Diversity')
% xlabel('Group Diversity Parameter')
% ylabel('Individual Diversity Parameter')
% zlabel('Average IFD')

%% Solve tasks
adivvals = logspace(-1, 3, 20);
gdivvals = logspace(-1, 3, 10);
DFD = NaN(20,10);
IFD = NaN(20,10);
minsn = NaN(20, 10);
maxsn = NaN(20, 10);
meansn = NaN(20, 10);
ai = 1;
for adiv = adivvals
    gi = 1;    
    for gdiv = gdivvals
        %Generate agents
        agents = GenAgent(par.numfuncs, par.numagents, [adiv, par.agspread], exp((-(0:par.numfuncs-1).^2)/gdiv), par.anorm);
        %Calculate and store diversity values
        [DFD(ai, gi), IFD(ai, gi)] = CalcFD(agents);
        
        %I will not regenerate agents for each run because div values then
        %would be different and I'd have to bin them (an extra few lines).
        %Instead with fixed agents tasks and task allocation is repeated.
        sn = NaN(par.numrepeats, 1);
        for ridx = 1:par.numrepeats
            %Generate tasks
            tasks = GenTask(par.numfuncs, par.numtasks, par.tnorm);
            %Assign tasks to agents
            t2a = randperm(par.numagents, size(tasks,1));
            sn(ridx) = SimpleSolveTasks(agents, tasks, t2a, par.EmergencyStop); %This returns number of steps required to solve all tasks, w/ upper limit
        end
        minsn(ai, gi) = min(sn);
        maxsn(ai, gi) = max(sn);
        meansn(ai, gi) = mean(sn);
        
        gi = gi + 1;
    end
    ai = ai + 1;
end

figure
surf(DFD, IFD, meansn)
title('Performance Plot')
xlabel('Dominant Functioanl Diversity')
ylabel('Average Individual Functional Diversity')
zlabel('Average Time of Completion [steps]')
