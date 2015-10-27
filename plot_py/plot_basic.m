

clf
chdir('~/models/PS_1D/results/exp')
p = load('profiles06.txt');
t = load('turbulence06.txt');
m = load('microphys06.txt');
%%%%  CONSTANTS %%%%%%%%%%%%%%%%%%%
hlat = 2.501e6;
cp = 1005.;

%%%% MEANS %%%%%%%%%%%%%%%%%%%%%%%%
Zm = p(:,1);
U = p(:,4);
TH = p(:,5);
Qc = m(:,6)*1000.;
Qv = p(:,7)*1000.;
THl = p(:,6);

%%% TURBULENCE %%%%%%%%%%%%%%%%%%%%
Zf = t(:,1);
Km = p(:,8);
Kt = p(:,9);
thw = t(:,2);
H = t(:,8);
H2 = t(:,9);
Bf = t(:,7);
Qtotf = t(:,10);
wq = t(:,6);
wqc = t(:,11);
wthv = t(:,12);
ri = p(:,10);
wqt = t(:,14);
wthl= t(:,13);


set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [21 29.7]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 21 29.7]);
set(gcf, 'renderer', 'painters');


subplot('Position',[0.15 0.77 0.36 0.18]);
plot(TH,Zm);
hold on;
plot(TH,Zm,'*');
axis([min(TH)-1 max(TH) 0 max(Zm)]);
xlabel('\theta (K)')

subplot('Position',[0.6 0.77 0.36 0.18]);
plot(THl,Zm);
hold on
plot(THl,Zm,'*');
axis([min(THl)-1 max(THl) 0 max(Zm)]);
xlabel('\theta_l (K)')

subplot('Position',[0.15 0.5 0.36 0.18]);
plot(Qv,Zm);
hold on
plot(Qv,Zm,'*');
axis([0 max(Qv) 0 max(Zm)]);
xlabel ('q_v (g/kg)')

subplot('Position',[0.6 0.5 0.36 0.18]);
plot(Qv+Qc,Zm);
hold on
plot(Qv+Qc,Zm,'*');
axis([0 max(Qv+Qc) 0 max(Zm)]);
xlabel('q_{total} (g/kg)')

subplot('Position',[0.15 0.23 0.36 0.18]);
plot(Qc,Zm);
hold on
plot(Qc,Zm,'*');
%axis([0 max(Qc)+0.1*max(Qc) 0 max(Zm)]);
xlabel('q_c (g/kg)')

subplot('Position',[0.6 0.23 0.36 0.18]);
plot(U,Zm);
hold on
plot(U,Zm,'*');
axis([min(U) max(U)+0.1*max(U) 0 max(Zm)]);
xlabel('|V| (m/s)')

print -f1 -depsc2 -tiff 'means';
clf


subplot('Position',[0.15 0.77 0.36 0.18]);

plot(Km,Zm);
hold on;
plot(Kt,Zm,'r');
axis([min(Km)-1 max(Kt)+5 0 max(Zm)]);
xlabel('Eddy diffusivities')
legend('K_m','K_t')

subplot('Position',[0.6 0.77 0.36 0.18]);
plot(thw,Zf);
axis([min(thw)-0.1*abs(min(thw)) max(thw)+0.1*max(thw) 0 max(Zf)]);
xlabel('\overline{w\theta}')


subplot('Position',[0.15 0.5 0.36 0.18]);
%plot(H,Zf);
%hold on
%plot(H2,Zf,'r');
axis([min(H)-0.1*abs(min(H)) max(H)+0.1*max(H) 0 max(Zf)]);
xlabel('H')

subplot('Position',[0.6 0.5 0.36 0.18]);
%plot(-(0.61*wq-wqc)*290+thw,Zf);
%plot(Bf,Zf)
%plot(wq*0.61*290,Zf)
%plot(wthv,Zf,'r')
plot(wqc,Zf,'*');
hold on
plot(wqt,Zf,'r');
%axis([min(wq-wqc)*290*0.61-0.1*abs(min((wq-wqc)*290*0.61)) max((wq-wqc)*290*0.61)+0.1*abs(max((wq-wqc)*290*0.61)) 0 max(Zf)]);
xlabel('wqc')

subplot('Position',[0.15 0.23 0.36 0.18]);
%plot(Qtotf,Zf);
plot(wq,Zf,'*');
hold on
plot(wqt,Zf,'r');
axis([min(wq)-0.1*abs(min(wq)) max(wq)+0.1*max(wq) 0 max(Zf)]);
xlabel('wq')

subplot('Position',[0.6 0.23 0.36 0.18]);
%plot(Bf,Zf);
%axis([min(Bf)-0.1*abs(min(Bf)) max(Bf)+0.1*max(Bf) 0 max(Zf)]);
%xlabel('Buoyancy flux')
%plot(wqc*290,Zf);
%plot(ri,Zm)
%hold on
%plot(ri,Zm,'*')
%xlabel ('Richardson number')
%plot(wthv,Zf)
print -f1 -depsc2 -tiff 'turb';