%% Radiative electron-hole recombination across the band gap
% Top: band diagram with electrons (blue, conduction band) and holes
% (red rings, valence band) doing random walks. When an electron and hole
% meet, they recombine and emit a photon with energy ~ Eg.
% Bottom: number of electron-hole pairs vs time, compared with the
% mean-field rate equation dN/dt = -B*N^2.
clear; clc; close all
rng(1)

%% Parameters (play with these)
Eg     = 1.12;     % band gap (eV), silicon-like
N0     = 60;       % initial electron-hole pairs (photo-excited)
nSteps = 400;      % number of time steps
sigma  = 0.03;     % random-walk step in position (device length = 1)
relax  = 0.12;     % fraction of excess energy lost per step (phonon emission)
rCap   = 0.05;     % capture distance for recombination
pRec   = 0.003;    % recombination probability per step for a close pair

%% Initial carriers: created by a light pulse, so they start "hot"
xe = rand(N0,1);   Ee = 0.5*rand(N0,1);   % electron position, energy above Ec
xh = rand(N0,1);   Eh = 0.4*rand(N0,1);   % hole position, energy below Ev
ev = zeros(0,2);   % recombination events: [x, age]
nPh = 0;           % photons emitted so far

%% Figure
figure('Position',[100 100 900 680],'Color','w');
ax1 = subplot(3,1,[1 2]); hold on
fill([0 1 1 0],[Eg Eg Eg+0.75 Eg+0.75],[0.80 0.88 1.00],'EdgeColor','none')
fill([0 1 1 0],[-0.55 -0.55 0 0],[1.00 0.85 0.85],'EdgeColor','none')
plot([0 1],[Eg Eg],'b','LineWidth',2)
plot([0 1],[0 0],'r','LineWidth',2)
text(1.01,Eg,'E_c','FontSize',12); text(1.01,0,'E_v','FontSize',12)
text(0.5,Eg/2,sprintf('E_g = %.2f eV',Eg),'HorizontalAlignment','center', ...
    'Color',[.4 .4 .4],'FontSize',12)
flash = plot(NaN,NaN,'-','Color',[1 0.6 0],'LineWidth',2);
hE = scatter(xe, Eg+Ee, 60, 'b', 'filled');
hH = scatter(xh, -Eh, 60, 'r', 'LineWidth', 1.5);
hP = scatter(NaN, NaN, 160, [1 0.7 0], 'p', 'filled');
xlim([0 1.08]); ylim([-0.7 Eg+0.9]); ylabel('Energy (eV)'); xlabel('position')
ttl = title('');
legend([hE hH hP],{'electron','hole','photon'},'Location','northeast')

ax2 = subplot(3,1,3); hold on; grid on
k    = pRec*(2*rCap - rCap^2);              % mean-field rate constant per step
Nth  = N0 ./ (1 + k*N0*(0:nSteps));
plot(0:nSteps, Nth, 'k--', 'LineWidth', 1.5)
hN = plot(0, N0, 'b', 'LineWidth', 2);
xlim([0 nSteps]); ylim([0 N0*1.05])
xlabel('time step'); ylabel('pairs left')
legend('mean-field  N_0/(1+BN_0t)','simulation')
Nhist = N0;

%% Time loop
for step = 1:nSteps
    % 1) thermal random walk in position, reflecting walls
    xe = abs(xe + sigma*randn(size(xe)));  xe(xe>1) = 2 - xe(xe>1);
    xh = abs(xh + sigma*randn(size(xh)));  xh(xh>1) = 2 - xh(xh>1);

    % 2) hot carriers lose energy to the lattice (electrons sink to Ec,
    %    holes float up to Ev)
    Ee = Ee*(1-relax) + 0.01*rand(size(Ee));
    Eh = Eh*(1-relax) + 0.01*rand(size(Eh));

    % 3) recombination: close electron-hole pairs may annihilate
    [ie, ih] = find(abs(xe - xh.') < rCap);
    usedE = false(size(xe));  usedH = false(size(xh));
    for q = randperm(numel(ie))
        if ~usedE(ie(q)) && ~usedH(ih(q)) && rand < pRec
            usedE(ie(q)) = true;  usedH(ih(q)) = true;
            ev(end+1,:) = [(xe(ie(q)) + xh(ih(q)))/2, 0]; %#ok<SAGROW>
            nPh = nPh + 1;
        end
    end
    xe(usedE) = [];  Ee(usedE) = [];
    xh(usedH) = [];  Eh(usedH) = [];

    % 4) age the events (photon flies upward and out of the crystal)
    ev(:,2) = ev(:,2) + 1;
    ev(ev(:,2) > 25, :) = [];

    % 5) update graphics
    hE.XData = xe;  hE.YData = Eg + Ee;
    hH.XData = xh;  hH.YData = -Eh;
    hP.XData = ev(:,1);  hP.YData = Eg/2 + 0.05*ev(:,2);
    f  = ev(ev(:,2) <= 4, 1);                 % flash the transition briefly
    fx = [f'; f'; nan(1,numel(f))];
    fy = [Eg*ones(1,numel(f)); zeros(1,numel(f)); nan(1,numel(f))];
    flash.XData = fx(:);  flash.YData = fy(:);

    Nhist(end+1) = numel(xe); %#ok<SAGROW>
    hN.XData = 0:step;  hN.YData = Nhist;
    ttl.String = sprintf('step %d  |  pairs left: %d  |  photons emitted: %d', ...
        step, numel(xe), nPh);

    drawnow
    pause(0.03)
    if isempty(xe) || isempty(xh), break, end
end