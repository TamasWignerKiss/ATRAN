%% Description
% This function shows how the IFD x DFD plane is populated in the simulations
%
%% Parameters
% Basic setup parameters
par.numfuncs = 9; %Number of domain functions
par.numagents = 100; %Number of agents in the system
par.agspread = 10; %A parameter specifying a meta-spread in agent similarity (could possibly be eliminated)
par.anorm = 10; %Total skill of agents

%% Some initialization
adivvals = 10.^(linspace(0.1, 0.2, 20).*(linspace(-6, 7, 20)));
gdivvals = 10.^(0.2*linspace(-3, 5, 10));
DFD = NaN(20,10);
IFD = NaN(20,10);

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
subplot(2,2,4)
surf(DFD, 'Marker', '.')
ylabel('Agent diversity parameter index')
xlabel('Group diversity parameter index')
zlabel('DFD')
subplot(2,2,3)
scatter(IFD(:), DFD(:), '.')
xlabel('IFD')
ylabel('DFD')
