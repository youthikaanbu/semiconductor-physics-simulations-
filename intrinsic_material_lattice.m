%% Intrinsic semiconductor, part 1: crystal -> thermal energy -> electron-hole pairs
% A 2D picture of the silicon crystal (as in Streetman & Banerjee): every atom
% shares its 4 valence electrons in 4 covalent bonds.
%   * At 0 K every bond is intact: nothing can move.
%   * Thermal energy occasionally breaks a bond. That frees ONE electron and
%     leaves ONE hole (the incomplete bond), so carriers are made in pairs.
%   * Free electrons wander; holes move when a neighboring bond electron
%     hops in. If an electron meets a hole the bond is restored (recombination).
%   * At each temperature, bonds breaking (G) and bonds restoring (R) balance,
%     which fixes how many pairs exist. That equilibrium number is n_i.
%
% NOTE: "temperature" here is the dimensionless kT/Eg. At real 300 K in Si,
% kT/Eg = 0.023 and only ~1 bond in 10^13 is broken, which you could never see,
% so this toy model uses hotter, scaled-up values. Part 2 uses real Si numbers.
clear; clc; close all
rng(3)

%% Settings (edit these)
L          = 10;                       % atoms per side
Eg         = 1.12;                     % band gap in eV (only used for labels)
tauLevels  = [0 0.07 0.09 0.11 0.13];  % kT/Eg at each stage (0 = absolute zero)
stepsPer   = 300;                      % animation steps at each temperature
A          = 15;                       % bond-breaking attempt prefactor (toy model)
pR         = 0.05;                     % recombination probability per step when close
rRec       = 0.6;                      % "close" = within this many lattice spacings
pHop       = 0.6;                      % chance per step that a hole hops to a neighbor bond
sE         = 0.35;                     % random-walk step of a free electron
pauseTime  = 0.05;                     % seconds per frame (bigger = slower)

%% Build the lattice: atoms and bonds
[AX,AY] = meshgrid(1:L,1:L);
atomX = AX(:);  atomY = AY(:);
atomId = reshape(1:L*L, L, L);          % atomId(row i, column j)
b1 = zeros(0,1);  b2 = zeros(0,1);
for i = 1:L
    for j = 1:L
        if j < L, b1(end+1,1) = atomId(i,j); b2(end+1,1) = atomId(i,j+1); end %#ok<SAGROW>
        if i < L, b1(end+1,1) = atomId(i,j); b2(end+1,1) = atomId(i+1,j); end %#ok<SAGROW>
    end
end
nB = numel(b1);
mx = (atomX(b1) + atomX(b2))/2;   my = (atomY(b1) + atomY(b2))/2;   % bond midpoints
dirx = atomX(b2) - atomX(b1);     diry = atomY(b2) - atomY(b1);
posA = [mx - 0.14*dirx, my - 0.14*diry];   % where the two bonding electrons sit
posB = [mx + 0.14*dirx, my + 0.14*diry];

% neighbor bonds (bonds that share an atom): a hole hops between these
Inc = sparse([b1; b2], [1:nB 1:nB]', 1, L*L, nB);
Adj = full(Inc'*Inc);  Adj(logical(eye(nB))) = 0;
nbr = cell(nB,1);
for b = 1:nB, nbr{b} = find(Adj(b,:) > 0); end

reflect = @(z) 0.5 + L*(1 - abs(1 - abs((z-0.5)/L)));   % keep electrons inside the crystal
totSteps = numel(tauLevels)*stepsPer;

%% Figure
fig = figure('Position',[40 40 1250 760],'Color','w');
hBan = annotation('textbox',[0.03 0.945 0.94 0.05],'String','','EdgeColor','none', ...
    'HorizontalAlignment','center','FontSize',16,'FontWeight','bold');
hExp = annotation('textbox',[0.03 0.885 0.94 0.055],'String','','EdgeColor','none', ...
    'HorizontalAlignment','center','FontSize',12,'Color',[0.25 0.25 0.25]);

% ---- crystal ----
ax1 = axes('Position',[0.03 0.09 0.47 0.78]); hold(ax1,'on')
plot(ax1,[atomX(b1) atomX(b2)]',[atomY(b1) atomY(b2)]','-','Color',[0.80 0.80 0.80],'LineWidth',1)
hBrk = plot(ax1,NaN,NaN,'-','Color',[0.9 0.2 0.2],'LineWidth',2.5);
scatter(ax1,atomX,atomY,650,[0.84 0.86 0.90],'filled','MarkerEdgeColor',[0.55 0.55 0.6]);
for k = 1:L*L
    text(ax1,atomX(k),atomY(k),'+4','HorizontalAlignment','center','FontSize',8,'Color',[0.3 0.3 0.4]);
end
hBE   = scatter(ax1,NaN,NaN,18,[0.15 0.15 0.15],'filled');           % bonding electrons
hHole = scatter(ax1,NaN,NaN,110,'r','LineWidth',2);                  % hole = missing bond electron
hGenF = scatter(ax1,NaN,NaN,420,[1 0.65 0],'p','filled');            % flash: bond just broke
hRecF = scatter(ax1,NaN,NaN,420,[0.2 0.7 0.3],'o','LineWidth',2.5);  % flash: bond restored
hEl1  = scatter(ax1,NaN,NaN,90,'b','filled','MarkerEdgeColor','k');  % free electrons
axis(ax1,'equal'); xlim(ax1,[0.4 L+0.6]); ylim(ax1,[0.4 L+0.6]); box(ax1,'on')
set(ax1,'XTick',[],'YTick',[])
title(ax1,'Silicon crystal, 2D bond picture   (black dots = bond electrons, +4 = ion cores)', ...
    'FontSize',10,'FontWeight','normal')
hCnt = text(ax1,(L+1)/2,-0.15,'','HorizontalAlignment','center','FontSize',12,'VerticalAlignment','top');

% ---- band diagram ----
ax2 = axes('Position',[0.55 0.60 0.42 0.27]); hold(ax2,'on')
fill(ax2,[0.5 L+0.5 L+0.5 0.5],[Eg Eg Eg+0.5 Eg+0.5],[0.82 0.90 1.00],'EdgeColor','none')
fill(ax2,[0.5 L+0.5 L+0.5 0.5],[-0.5 -0.5 0 0],[1.00 0.88 0.88],'EdgeColor','none')
plot(ax2,[0.5 L+0.5],[Eg Eg],'b','LineWidth',2)
plot(ax2,[0.5 L+0.5],[0 0],'r','LineWidth',2)
text(ax2,L+0.6,Eg,'E_c','FontSize',11); text(ax2,L+0.6,0,'E_v','FontSize',11)
text(ax2,L/2+0.5,Eg/2,sprintf('E_g = %.2f eV',Eg),'HorizontalAlignment','center','Color',[0.45 0.45 0.45])
text(ax2,0.6,Eg+0.38,'conduction band','FontSize',9,'Color',[0.3 0.4 0.6])
text(ax2,0.6,-0.4,'valence band','FontSize',9,'Color',[0.6 0.3 0.3])
hArrG  = plot(ax2,NaN,NaN,'-','Color',[1 0.55 0],'LineWidth',2);
hHeadG = plot(ax2,NaN,NaN,'^','Color',[1 0.55 0],'MarkerFaceColor',[1 0.55 0],'MarkerSize',6);
hArrR  = plot(ax2,NaN,NaN,'-','Color',[0.2 0.7 0.3],'LineWidth',2);
hHeadR = plot(ax2,NaN,NaN,'v','Color',[0.2 0.7 0.3],'MarkerFaceColor',[0.2 0.7 0.3],'MarkerSize',6);
hE2 = scatter(ax2,NaN,NaN,50,'b','filled');
hH2 = scatter(ax2,NaN,NaN,50,'r','LineWidth',1.6);
xlim(ax2,[0.5 L+1.1]); ylim(ax2,[-0.5 Eg+0.5]); box(ax2,'on'); set(ax2,'XTick',[])
ylabel(ax2,'Energy (eV)')
title(ax2,'Same carriers on the band diagram (orange = pair created, green = pair recombines)', ...
    'FontSize',10,'FontWeight','normal')

% ---- n and p vs time ----
ax3 = axes('Position',[0.55 0.34 0.42 0.17]); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on')
for k = 1:numel(tauLevels)
    plot(ax3,[(k-1)*stepsPer (k-1)*stepsPer],[0 80],':k')
    text(ax3,(k-0.5)*stepsPer,73,sprintf('kT/E_g=%.2f',tauLevels(k)),'HorizontalAlignment','center','FontSize',8)
end
hn = plot(ax3,NaN,NaN,'b-','LineWidth',3.5);
hp = plot(ax3,NaN,NaN,'r--','LineWidth',1.8);
xlim(ax3,[0 totSteps]); ylim(ax3,[0 80]); ylabel(ax3,'carriers in crystal')
legend(ax3,[hn hp],{'n (free electrons)','p (holes)'},'Location','west','FontSize',9)
title(ax3,'n and p always equal: pairs are the only source','FontSize',10,'FontWeight','normal')

% ---- generation vs recombination ----
ax4 = axes('Position',[0.55 0.09 0.42 0.17]); hold(ax4,'on'); box(ax4,'on'); grid(ax4,'on')
hG = plot(ax4,NaN,NaN,'-','Color',[1 0.55 0],'LineWidth',2.5);
hR = plot(ax4,NaN,NaN,'-','Color',[0.2 0.7 0.3],'LineWidth',2.5);
xlim(ax4,[0 totSteps]); ylim(ax4,[0 2]); xlabel(ax4,'time step'); ylabel(ax4,'events per step')
legend(ax4,[hG hR],{'bonds breaking (G)','bonds restoring (R)'},'Location','northwest','FontSize',9)
title(ax4,'At each temperature G = R: that balance sets n_i','FontSize',10,'FontWeight','normal')

%% State
isHole = false(nB,1);  holeB = zeros(0,1);  holeU = zeros(0,1);
E   = zeros(0,3);                 % free electrons: [x y u]
evG = zeros(0,3);  evR = zeros(0,3);   % flashes: [x y age]
nHist = zeros(1,totSteps);  pHist = nHist;  gHist = nHist;  rHist = nHist;
step = 0;

%% Animation
for lev = 1:numel(tauLevels)
    tau = tauLevels(lev);
    pG  = A*exp(-1/tau);                       % bond-breaking probability per step (0 when tau = 0)
    for s = 1:stepsPer
        step = step + 1;

        % 1) thermal generation: intact bonds may break, making an electron AND a hole
        cand = find(~isHole);
        newB = cand(rand(numel(cand),1) < pG);
        for q = 1:numel(newB)
            b = newB(q);
            isHole(b) = true;  holeB(end+1,1) = b;  holeU(end+1,1) = rand; %#ok<SAGROW>
            E(end+1,:) = [mx(b) my(b) rand];                              %#ok<SAGROW>
            evG(end+1,:) = [mx(b) my(b) 0];                               %#ok<SAGROW>
        end
        gHist(step) = numel(newB);

        % 2) motion: electrons random-walk, holes hop between neighboring bonds
        if ~isempty(E)
            E(:,1) = reflect(E(:,1) + sE*randn(size(E,1),1));
            E(:,2) = reflect(E(:,2) + sE*randn(size(E,1),1));
        end
        for k = 1:numel(holeB)
            if rand < pHop
                nb = nbr{holeB(k)};
                c  = nb(randi(numel(nb)));
                if ~isHole(c)
                    isHole(holeB(k)) = false;  isHole(c) = true;  holeB(k) = c;
                end
            end
        end

        % 3) recombination: an electron next to a hole can fill it (bond restored)
        nRec = 0;
        if ~isempty(E) && ~isempty(holeB)
            D = hypot(E(:,1) - mx(holeB).', E(:,2) - my(holeB).');
            [ie, ih] = find(D < rRec);
            usedE = false(size(E,1),1);  usedH = false(numel(holeB),1);
            for q = randperm(numel(ie))
                if ~usedE(ie(q)) && ~usedH(ih(q)) && rand < pR
                    usedE(ie(q)) = true;  usedH(ih(q)) = true;
                    evR(end+1,:) = [mx(holeB(ih(q))) my(holeB(ih(q))) 0]; %#ok<SAGROW>
                    nRec = nRec + 1;
                end
            end
            isHole(holeB(usedH)) = false;
            holeB(usedH) = [];  holeU(usedH) = [];  E(usedE,:) = [];
        end
        rHist(step) = nRec;
        nHist(step) = size(E,1);  pHist(step) = numel(holeB);

        % 4) age the flashes
        evG(:,3) = evG(:,3) + 1;  evG(evG(:,3) > 10,:) = [];
        evR(:,3) = evR(:,3) + 1;  evR(evR(:,3) > 10,:) = [];

        % ---- draw crystal ----
        keep = ~isHole;
        set(hBE,'XData',[posA(:,1); posB(keep,1)],'YData',[posA(:,2); posB(keep,2)])
        set(hHole,'XData',posB(holeB,1),'YData',posB(holeB,2))
        sx = [atomX(b1(holeB))'; atomX(b2(holeB))'; nan(1,numel(holeB))];
        sy = [atomY(b1(holeB))'; atomY(b2(holeB))'; nan(1,numel(holeB))];
        set(hBrk,'XData',sx(:),'YData',sy(:))
        set(hEl1,'XData',E(:,1),'YData',E(:,2))
        set(hGenF,'XData',evG(:,1),'YData',evG(:,2))
        set(hRecF,'XData',evR(:,1),'YData',evR(:,2))
        hCnt.String = sprintf('free electrons n = %d      holes p = %d      broken bonds: %d of %d', ...
            size(E,1), numel(holeB), numel(holeB), nB);

        % ---- draw band diagram ----
        set(hE2,'XData',E(:,1),'YData',Eg + 0.05 + 0.3*E(:,3))
        set(hH2,'XData',mx(holeB),'YData',-0.05 - 0.3*holeU)
        gx = evG(:,1)';  rx = evR(:,1)';
        gX = [gx; gx; nan(size(gx))];  gY = [zeros(size(gx)); Eg*ones(size(gx)); nan(size(gx))];
        rX = [rx; rx; nan(size(rx))];  rY = [Eg*ones(size(rx)); zeros(size(rx)); nan(size(rx))];
        set(hArrG,'XData',gX(:),'YData',gY(:));   set(hHeadG,'XData',gx,'YData',Eg*ones(size(gx)))
        set(hArrR,'XData',rX(:),'YData',rY(:));   set(hHeadR,'XData',rx,'YData',zeros(size(rx)))

        % ---- time series ----
        set(hn,'XData',1:step,'YData',nHist(1:step))
        set(hp,'XData',1:step,'YData',pHist(1:step))
        w = 40;
        set(hG,'XData',1:step,'YData',movmean(gHist(1:step),[w-1 0]))
        set(hR,'XData',1:step,'YData',movmean(rHist(1:step),[w-1 0]))

        % ---- text ----
        if tau == 0
            hBan.String = sprintf('Stage %d of %d:   T = 0 K', lev, numel(tauLevels));
            hBan.Color  = [0.10 0.35 0.65];
            hExp.String = ['Every bond is intact: valence band full, conduction band empty. ' ...
                'No electron can move, so there are no free carriers and the crystal is an insulator.'];
        else
            hBan.String = sprintf('Stage %d of %d:   kT/E_g = %.2f   (hotter)', lev, numel(tauLevels), tau);
            hBan.Color  = [0.75 0.30 0.10];
            hExp.String = ['Thermal energy breaks bonds at random (it takes about E_g). Each break frees ONE electron and ' ...
                'leaves ONE hole, so n = p. Electrons and holes that meet recombine; the two rates balance (G = R).'];
        end

        drawnow
        pause(pauseTime)
    end
end