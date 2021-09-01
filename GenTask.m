function tasks = GenTask(nf, nt, no)

%% Generate tasks
tasks = rand(nt, nf);
tasks = tasks./sum(tasks, 2)*no;

end
