%% Electron wave packet on a 1D atomic chain, with electric field
N  = 200;  t = 1;  F = 0.05;
n  = (1:N)';
H  = -t*(diag(ones(N-1,1),1) + diag(ones(N-1,1),-1)) + diag(F*n);
% H(100,100) = H(100,100) - 1.5;      % uncomment: add an impurity atom

x0 = 140; sigma = 8; k0 = 0;
psi = exp(-(n-x0).^2/(4*sigma^2)) .* exp(1i*k0*n);
psi = psi/norm(psi);

dt = 0.5;  U = expm(-1i*H*dt);

figure('Position',[100 100 1000 450]);
subplot(2,1,1); hold on
plot(n, zeros(N,1), 'o', 'Color',[.6 .6 .6], 'MarkerFaceColor',[.8 .8 .8])
prob = area(n, abs(psi).^2, 'FaceAlpha',0.5, 'FaceColor','r');
xlim([1 N]); ylim([0 0.06]); ylabel('|\psi|^2'); xlabel('lattice site')
ttl = title('');

subplot(2,1,2); hold on
plot(n, F*n - 2*t, 'b'); plot(n, F*n + 2*t, 'b')
cm = plot(x0, F*x0, 'ro', 'MarkerFaceColor','r');
xlim([1 N]); ylabel('E'); xlabel('site (band bottom/top vs. position)')

for step = 1:600
    psi = U*psi;
    prob.YData = abs(psi).^2;
    xm = sum(n.*abs(psi).^2);
    cm.XData = xm;  cm.YData = F*xm - 2*t;
    ttl.String = sprintf('t = %.0f', step*dt);
    drawnow
    pause(0.03)
end