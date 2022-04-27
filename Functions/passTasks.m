function [t2a, numPass] = passTasks(t2a, agents, tasks, simTh, numPass)

%Sanity check
if length(t2a) ~= size(tasks, 1)
    error('Something went wrong and more agents are assigned to tasks than tasks available.')
end

% For each task evaluate if the agent working on it should pass it or keep it
noTasks = size(tasks, 1);
for tidx = 1:noTasks
    
    %Check if task was not passed too many times and allow further pass-evaluation only if passes are still allowed
    if numPass(tidx) > 0
        
        %Who is similar enough to be worthy for communication? -- asked the agent working to task tidx?
        %simAgs = find(pdist2(agents(t2a(tidx), :), agents) < simTh);
        simAgs = find(sqrt(sum((agents - agents(t2a(tidx), :)).^2, 2)) < simTh); %This is more efficient
        
        %Which are the free similar agents?
        freeAgs = simAgs(not(ismember(simAgs, t2a)));
        
        %How good am I to solve my task?
        %myFit = pdist2(agents(t2a(tidx), :), tasks(tidx,:));
        myFit = sqrt(sum((agents(t2a(tidx), :) - tasks(tidx,:)).^2)); %Again, this is more than 10x faster
        
        %How good are free agents I talk to?
        %theirFit = pdist2(agents(freeAgs, :), tasks(tidx,:));
        theirFit = sqrt(sum((agents(freeAgs, :) - tasks(tidx, :)).^2, 2));
        
        %Who's best?
        [~, bestidx] = min([myFit; theirFit]);
        
        %If I'm not the best, re-assign my task to the best agent
        if bestidx ~= 1
            t2a(tidx) = freeAgs(bestidx-1);
            numPass(tidx) = numPass(tidx) - 1;
        end
    end
end
