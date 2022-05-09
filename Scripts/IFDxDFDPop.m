%% Description
% This function shows how the IFD x DFD plane is populated in the simulations
%
%% Parameters
par.numfuncs = 9; %Number of domain functions
par.numagents = 10; %Number of agents in the system
par.agspread = 10; %A parameter specifying a meta-spread in agent similarity (could possibly be eliminated)
par.anorm = 10; %Total skill of agents
par.noavals = 40; %Number of agent diversity steps
par.nogvals = 40; %Number of group diversity steps

%% Some initialization
adivvals = 10.^(linspace(0.1, 0.2, par.noavals).*(linspace(-6, 7, par.noavals)));
gdivvals = 10.^(0.2*linspace(-3, 5, par.nogvals));
DFD = NaN(par.noavals,par.nogvals);
IFD = NaN(par.noavals,par.nogvals);

%% Loop through DFD and IFD values
ai = 1;
for adiv = adivvals
    gi = 1;
    fprintf('\n%i ', ai)
    for gdiv = gdivvals
        %Generate agents
        agents = GenAgent(par.numfuncs, par.numagents, [adiv, par.agspread], exp((-(0:par.numfuncs-1).^2)/gdiv), par.anorm);
        %Calculate and store diversity values
        [DFD(ai, gi), IFD(ai, gi)] = CalcFD(agents);
         gi = gi + 1;
        fprintf('.')
    end
    ai = ai + 1;
end
fprintf('\n')
clear ai adiv gi gdiv agents

%% Plot IFD x DFD space
figure
subplot(2,2,1)
plot(adivvals, '-*')
hold on
plot(gdivvals, '-*')
grid on
legend('Agent diversity parameter', 'Group diversity parameter')
xlabel('Index')
ylabel('Value')

subplot(2,2,2)
surf(IFD, 'Marker', '.')
ylabel('Agent diversity parameter index')
xlabel('Group diversity parameter index')
zlabel('IFD')
zlim([0, 1])

subplot(2,2,4)
surf(DFD, 'Marker', '.')
ylabel('Agent diversity parameter index')
xlabel('Group diversity parameter index')
zlabel('DFD')
zlim([0, 1])

subplot(2,2,3)
scatter(IFD(:), DFD(:), '.')
xlabel('IFD')
ylabel('DFD')
grid on
xlim([0, 1])
ylim([0, 1])
