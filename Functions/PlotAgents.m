function fh = PlotAgents(agents)
%This is a helper function, to quickly visualize what agents of a simulation instance look like.
%
%Usage: PlotAgents(agents)

%% Plot agents
fh = figure('WindowStyle','docked', 'Name','Agents');
bh = bar(agents, 'stacked');

set(bh, 'FaceColor', 'Flat')
colors =  mat2cell(jet(numel(bh)),ones(numel(bh),1), 3); 
set(bh, {'CData'}, colors) % using jet colormap

xlabel('Agent''s index')
ylabel('Skills')

legend(strcat("Skill #", num2str((1:numel(bh))')))
grid on
