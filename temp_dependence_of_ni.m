%% Temperature dependence of extrinisic carrier concentration in n-type silicon
% Im applying cold to hot to the semiconductor silicon to see how the carrier concentration is affected
%   top left  : band diagram (donor level, electrons, holes, thermal generation)
%   top right : the material itself (what is actually free to move)
%   bottom    : n0 vs 1000/T curve with a marker at the current temperature
% The curve is computed by solving charge neutrality
%   n = p + Nd+      with   Nd+ = Nd / (1 + 2 (n/Nc) exp(dE/kT))
clear; clc; close all
rng(2)

%% Physical parameters (silicon doped with (P) phosphorus-like donors)
kB    = 8.617e-5;     % Boltzmann constant (eV/K)
Eg    = 1.12;         % band gap of silicon (eV)
dE    = 0.05;         % donor binding energy Ec - Ed (eV)
Nd    = 1e15;         % donor concentration (cm^-3)
Nc300 = 2.8e19;       % conduction-band effective DOS at 300 K (cm^-3)
Nv300 = 1.04e19;      % valence-band effective DOS at 300 K (cm^-3)

%% Solve charge neutrality at every temperature
nF   = 400;                                  % number of animation frames
Tvec = logspace(log10(35), log10(1000), nF); % cold -> hot
nvec = zeros(1,nF); pvec = nvec; Ndp = nvec; nivec = nvec;
for i = 1:nF
    T  = Tvec(i);  kT = kB*T;
    Nc = Nc300*(T/300)^1.5;   Nv = Nv300*(T/300)^1.5;
    ni2 = Nc*Nv*exp(-Eg/kT);
    f = @(x) exp(x) - ni2*exp(-x) - Nd/(1 + 2*exp(x)/Nc*exp(dE/kT));
    x = fzero(f, [-50 50]);                  % x = ln(n)
    nvec(i)  = exp(x);
    pvec(i)  = ni2/nvec(i);
    Ndp(i)   = Nd/(1 + 2*nvec(i)/Nc*exp(dE/kT));
    nivec(i) = sqrt(ni2);
end

% regime boundaries: n/Nd < 0.9 = freeze-out, 0.9..1.1 = extrinsic, >1.1 = intrinsic
r  = nvec/Nd;
i1 = find(r >= 0.9, 1);   i2 = find(r > 1.1, 1);
xb1 = 1000/Tvec(i1);      xb2 = 1000/Tvec(i2);
xmin = 0.9;  xmax = 1000/Tvec(1) + 1;

%% Display settings
M     = 12;           % donor atoms drawn
maxE  = 80;           % max electrons / holes drawn
drawEd = Eg - 0.12;   % donor level as DRAWN (real value is 0.05 eV, exaggerated for clarity)
xdb   = linspace(0.07,0.93,M)';                       % donor x positions, band diagram
[DX,DY] = meshgrid(linspace(0.15,0.85,4), linspace(0.2,0.8,3));
dxb = DX(:);  dyb = DY(:);                            % donor positions, material box
ord = randperm(M);                                    % order in which donors ionize
xE = rand(maxE,1); yE = rand(maxE,1); uE = rand(maxE,1);   % electron pool
xH = rand(maxE,1); yH = rand(maxE,1); uH = rand(maxE,1);   % hole pool
xth = [0.2 0.8 0.33 0.67 0.08 0.92];                  % x of thermal-generation arrows
refl = @(z) 1 - abs(1 - abs(z));                      % reflect into [0,1]

%% Figure
fig = figure('Position',[60 40 1150 760],'Color','w');
hBan = annotation('textbox',[0.03 0.945 0.94 0.05],'String','','EdgeColor','none', ...
    'HorizontalAlignment','center','FontSize',16,'FontWeight','bold');
hExp = annotation('textbox',[0.03 0.885 0.94 0.055],'String','','EdgeColor','none', ...
    'HorizontalAlignment','center','FontSize',12,'Color',[0.25 0.25 0.25]);

% ---- band diagram ----
ax1 = axes('Position',[0.07 0.43 0.45 0.42]); hold(ax1,'on')
fill(ax1,[0 1 1 0],[Eg Eg 1.7 1.7],[0.82 0.90 1.00],'EdgeColor','none')
fill(ax1,[0 1 1 0],[-0.35 -0.35 0 0],[1.00 0.88 0.88],'EdgeColor','none')
plot(ax1,[0 1],[Eg Eg],'b','LineWidth',2)
plot(ax1,[0 1],[0 0],'r','LineWidth',2)
plot(ax1,[0 1],[drawEd drawEd],'--','Color',[0 0.6 0.6],'LineWidth',1.2)
text(ax1,1.02,Eg,'E_c','FontSize',12); text(ax1,1.02,0,'E_v','FontSize',12)
text(ax1,1.02,drawEd,'E_D','FontSize',12,'Color',[0 0.5 0.5])
text(ax1,0.5,Eg/2+0.05,sprintf('E_g = %.2f eV',Eg),'HorizontalAlignment','center', ...
    'Color',[0.4 0.4 0.4],'FontSize',12)
text(ax1,0.01,drawEd-0.1,'donor level drawn larger than real 0.05 eV','FontSize',9,'Color',[0.45 0.45 0.45])
hArrI  = plot(ax1,NaN,NaN,'-','Color',[0.25 0.45 0.75],'LineWidth',1.5);
hHeadI = plot(ax1,NaN,NaN,'^','Color',[0.25 0.45 0.75],'MarkerFaceColor',[0.25 0.45 0.75],'MarkerSize',6);
hArrT  = plot(ax1,NaN,NaN,'-','Color',[0.85 0.3 0.25],'LineWidth',1.5);
hHeadT = plot(ax1,NaN,NaN,'^','Color',[0.85 0.3 0.25],'MarkerFaceColor',[0.85 0.3 0.25],'MarkerSize',6);
hDonB  = scatter(ax1,NaN,NaN,45,'b','filled');
hDonSq = scatter(ax1,NaN,NaN,110,'k','s','LineWidth',1.2);
hDonPl = scatter(ax1,NaN,NaN,45,'k','+','LineWidth',1.2);
hE1    = scatter(ax1,NaN,NaN,45,'b','filled');
hH1    = scatter(ax1,NaN,NaN,45,'r','LineWidth',1.4);
xlim(ax1,[0 1.1]); ylim(ax1,[-0.35 1.7]); box(ax1,'on')
ylabel(ax1,'Energy (eV)'); xlabel(ax1,'position in crystal')
text(ax1,0.60,-0.17,'valence band filled with electrons','FontSize',9,'Color',[0.45 0.45 0.45])
text(ax1,0.89,0.33,'thermal generation','FontSize',9,'Color',[0.85 0.3 0.25],'HorizontalAlignment','right')
legend(ax1,[hE1 hH1 hDonB hDonSq],{'mobile electron','mobile hole','bound electron','ionized donor'}, ...
    'Location','northoutside','Orientation','horizontal','FontSize',9)

% ---- material box: what is free to move ----
ax2 = axes('Position',[0.58 0.43 0.38 0.42]); hold(ax2,'on')
rectangle(ax2,'Position',[0 0 1 1],'EdgeColor',[0.2 0.35 0.55],'LineWidth',1.5)
hDonN2  = scatter(ax2,NaN,NaN,110,[0.7 0.7 0.7],'s');
hDonSq2 = scatter(ax2,NaN,NaN,110,'k','s','LineWidth',1.2);
hDonPl2 = scatter(ax2,NaN,NaN,45,'k','+','LineWidth',1.2);
hE2     = scatter(ax2,NaN,NaN,60,'b','filled');
hH2     = scatter(ax2,NaN,NaN,60,'r','LineWidth',1.4);
hNone   = text(ax2,0.5,0.5,'','HorizontalAlignment','center','FontSize',13, ...
    'Color',[0.4 0.4 0.4],'FontAngle','italic');
xlim(ax2,[-0.02 1.02]); ylim(ax2,[-0.02 1.02]); set(ax2,'XTick',[],'YTick',[])
title(ax2,'What is actually free to move (gray squares = neutral donors, boxed + = ionized, fixed in lattice)', ...
    'FontSize',9,'FontWeight','normal')
hNum = text(ax2,0.5,-0.13,'','HorizontalAlignment','center','FontSize',11,'VerticalAlignment','top');

% ---- n0 vs 1000/T ----
ax3 = axes('Position',[0.08 0.07 0.86 0.29]); hold(ax3,'on'); box(ax3,'on')
set(ax3,'YScale','log'); ylo = 1e11; yhi = 1e18;
patch(ax3,[xb1 xmax xmax xb1],[ylo ylo yhi yhi],[0.3 0.5 0.9],'FaceAlpha',0.12,'EdgeColor','none')
patch(ax3,[xb2 xb1 xb1 xb2],[ylo ylo yhi yhi],[0.3 0.8 0.4],'FaceAlpha',0.12,'EdgeColor','none')
patch(ax3,[xmin xb2 xb2 xmin],[ylo ylo yhi yhi],[0.9 0.4 0.3],'FaceAlpha',0.12,'EdgeColor','none')
text(ax3,mean([xb1 xmax]),4e17,'Freeze-out / ionization','HorizontalAlignment','center','FontSize',11)
text(ax3,mean([xb2 xb1]),4e17,'Extrinsic (plateau)','HorizontalAlignment','center','FontSize',11)
text(ax3,mean([xmin xb2]),1e12,'Intrinsic','Rotation',90,'HorizontalAlignment','left','FontSize',11)
plot(ax3,[xmin xmax],[Nd Nd],':','Color',[0.4 0.4 0.4])
text(ax3,xmax-0.3,Nd*1.5,'N_d','HorizontalAlignment','right','Color',[0.4 0.4 0.4])
plot(ax3,1000./Tvec,nivec,'k--','LineWidth',1.3)
text(ax3,3.1,3e12,'n_i','FontSize',12)
plot(ax3,1000./Tvec,nvec,'-','Color',[0.75 0.85 0.95],'LineWidth',2)
hRev = plot(ax3,NaN,NaN,'-','Color',[0.1 0.5 0.85],'LineWidth',2.5);
hMk  = plot(ax3,NaN,NaN,'o','MarkerFaceColor',[0.9 0.2 0.1],'MarkerEdgeColor','k','MarkerSize',11);
xlim(ax3,[xmin xmax]); ylim(ax3,[ylo yhi])
xlabel(ax3,'1000 / T  (K^{-1})        hot  \leftarrow                    \rightarrow  cold')
ylabel(ax3,'n_0  (cm^{-3})')

%% Optional: record a video (uncomment these and the writeVideo line below)
% vw = VideoWriter('carrier_vs_temperature.mp4','MPEG-4'); vw.FrameRate = 30; open(vw)

%% Animation loop: temperature rises from cold to hot
for i = 1:nF
    T = Tvec(i);  kT = kB*T;
    n = nvec(i);  p = pvec(i);

    % how many dots to draw (1 dot = a fixed concentration, enlarged only
    % when carriers become too numerous to draw)
    scale = max(Nd/M, n/maxE);
    nE = min(round(n/scale), maxE);
    nH = min(round(p/scale), maxE);
    kIon = round(M*Ndp(i)/Nd);
    ion = ord(1:kIon);  neu = ord(kIon+1:end);

    % thermal jitter of mobile carriers (faster when hotter)
    s  = 0.01 + 0.03*sqrt(T/300);
    xE = refl(xE + s*randn(maxE,1));      yE = refl(yE + s*randn(maxE,1));
    xH = refl(xH + 0.5*s*randn(maxE,1));  yH = refl(yH + 0.5*s*randn(maxE,1));

    % band diagram: donors
    set(hDonB, 'XData',xdb(neu),'YData',drawEd*ones(numel(neu),1))
    set(hDonSq,'XData',xdb(ion),'YData',drawEd*ones(numel(ion),1))
    set(hDonPl,'XData',xdb(ion),'YData',drawEd*ones(numel(ion),1))
    xa = xdb(ion)';
    ax_ = [xa; xa; nan(size(xa))];
    ay_ = [drawEd*ones(size(xa)); Eg*ones(size(xa)); nan(size(xa))];
    set(hArrI,'XData',ax_(:),'YData',ay_(:))
    set(hHeadI,'XData',xa,'YData',Eg*ones(size(xa)))

    % band diagram: thermal generation arrows (valence -> conduction band)
    nTh = min(nH,6);
    xt  = xth(1:nTh);
    tx_ = [xt; xt; nan(size(xt))];
    ty_ = [zeros(size(xt)); Eg*ones(size(xt)); nan(size(xt))];
    set(hArrT,'XData',tx_(:),'YData',ty_(:))
    set(hHeadT,'XData',xt,'YData',Eg*ones(size(xt)))

    % band diagram: mobile carriers
    set(hE1,'XData',xE(1:nE),'YData',Eg + 0.02 + 1.5*kT*uE(1:nE))
    set(hH1,'XData',xH(1:nH),'YData',-0.02 - 1.5*kT*uH(1:nH))

    % material box
    set(hDonN2, 'XData',dxb(neu),'YData',dyb(neu))
    set(hDonSq2,'XData',dxb(ion),'YData',dyb(ion))
    set(hDonPl2,'XData',dxb(ion),'YData',dyb(ion))
    set(hE2,'XData',xE(1:nE),'YData',yE(1:nE))
    set(hH2,'XData',xH(1:nH),'YData',yH(1:nH))
    if nE == 0, hNone.String = 'nothing mobile'; else, hNone.String = ''; end

    % curve and marker
    set(hRev,'XData',1000./Tvec(1:i),'YData',nvec(1:i))
    set(hMk, 'XData',1000/T,'YData',n)

    % labels
    if r(i) < 0.9
        name = 'FREEZE-OUT / IONIZATION';  col = [0.10 0.35 0.65];
        msg  = ['Not enough thermal energy to free the electrons from the donors. Most stay bound ' ...
                'at E_D, so few carriers are mobile.'];
    elseif r(i) <= 1.1
        name = 'EXTRINSIC (all donors ionized)';  col = [0.10 0.50 0.30];
        msg  = ['Every donor has released its electron, so n \approx N_d and barely changes with T. ' ...
                'This is the regime devices are designed for.'];
    else
        name = 'INTRINSIC (thermal generation takes over)';  col = [0.75 0.25 0.20];
        msg  = ['kT is now large enough to lift electrons across the whole gap. Electron-hole pairs ' ...
                'outnumber the donor electrons, so n \approx p \approx n_i.'];
    end
    hBan.String = sprintf('T = %.0f K    (1000/T = %.1f)    %s', T, 1000/T, name);
    hBan.Color  = col;
    hExp.String = msg;
    hNum.String = sprintf('n = %.1e cm^{-3}    p = %.1e cm^{-3}    ionized donors: %.0f%%\n1 dot \\approx %.1e cm^{-3}', ...
        n, p, 100*Ndp(i)/Nd, scale);

    drawnow
    pause(0.02)
    % writeVideo(vw, getframe(fig))
end
% close(vw)