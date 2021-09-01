function stepno = SimpleSolveTasks(agents, tasks, t2a, emstop)
%% This is a way to implement solving tasks
stepno = 1;
while any(any(tasks > 0)) && stepno <= emstop
    tasks = tasks - agents(t2a,:);
    stepno = stepno + 1;
end
stepno = stepno -1;

end
