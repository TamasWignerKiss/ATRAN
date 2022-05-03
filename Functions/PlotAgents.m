function fh = PlotAgents(agents)

fh = figure('WindowStyle','docked', 'Name','Agents');
bar(agents, 'stacked')
xlabel('Agent''s index')
ylabel('Skills')