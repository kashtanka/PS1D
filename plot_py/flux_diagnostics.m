clf
chdir('~/models/PS_1D/results/exp')
p = load('profiles06.txt');
t = load('turbulence06.txt');
m = load('microphys06.txt');

%%%% MEANS %%%%%%%%%%%%%%%%%%%%%%%%
Zm = p(:,1);
U = p(:,4);
TH = p(:,5);
Qc = m(:,6);
Qv = p(:,7);
THl = p(:,6);

%%% TURBULENCE %%%%%%%%%%%%%%%%%%%%
Zf = t(:,1);
Km = p(:,8);
Kt = p(:,9);
thw = t(:,2);
wth_l = t(:,9);
Bf = t(:,7);
wq_t = t(:,10);
wq = t(:,6);
wqc = t(:,11);
wthv = t(:,12);
ri = p(:,10);

%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hlat=2.501e6;
r=287.05;
cp=1005.;
g=9.8066;
akapa=r/cp;
p00=1.e5;
rv=461.51;

%%%%%%%%%% THERMODYNAMICS AND MOISTURE %%%%%%%%%%%%%%%%
Qs = m(:,3);
temp = p(:,17);
pres = p(:,18);
wq2(1:179)=0;
wq3(1:179)=0;
dqdz = p(:,14);

for i = 2:size(TH,1)-2
    dthdz = (TH(i) - TH(i-1))/(Zm(i)- Zm(i-1));
    Kt_hl = 0.5*(Kt(i) + Kt(i-1));
    t_hl = 0.5*(temp(i) + temp(i-1));
    th_hl =  0.5*(TH(i) + TH(i-1));
    qs_hl = 0.5*(Qs(i) + Qs(i-1));
    alpha = qs_hl*hlat/rv/t_hl/t_hl;
    pot = (0.5*(pres(i)+pres(i-1))/p00)^akapa;
%%%%%%% SPECIFIC HUMIDITY FLUX AT SATURATION %%%%%%%%%%
    wq2(i) = -Kt_hl*(alpha*dthdz*pot - alpha*g/cp+g/r/t_hl*qs_hl);
    dqsdz = (Qs(i) - Qs(i-1))/(Zm(i) - Zm(i-1));
%    dqdz2 = (Qv(i) - Qv(i-1))/(Zm(i) - Zm(i-1));
    wq3(i) = -Kt_hl*dqsdz;
    wq4(i) = -Kt_hl*(alpha*dthdz*pot - alpha*g/cp);
    wq5 (i) = -Kt_hl*(alpha*dthdz*pot);
    wq6(i) = -Kt_hl*(dqsdz+qs_hl*g/r/t_hl);
    wq7(i) = alpha*(pot*wth_l(i)-hlat/cp*wq_t(i))/(1+alpha*hlat/cp);

%%%%%%%%%% THETA FLUX FROM FLUXES OF CONSERVED VARIABLES %%%%%%%%

 %   wth6(i) = -Kt_hl*(dthdz - g/cp/pot +g/r/temp(i)*Qs(i)/pot/alpha);
 %   wth(i) = (H2(i) + pot^(-1)*hlat/cp*wq_t(i))/(1 + hlat/cp*alpha);
    if (0.5*(Qc(i)+Qc(i-1)) > 0)
        wth(i) = (wth_l(i) - pot^(-1)*hlat/cp*wq_t(i))/(1 + hlat/cp*alpha);
    wthv2(i) =  wth6(i) + 0.5*(TH(i)+TH(i-1))*(0.61*wq(i));
    else
        wth(i)= wth_l(i);
    wthv2(i) = thw(i) -0.5*(TH(i)+TH(i-1))*(0.61*wq(i));
    end
end
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [21 29.7]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 21 29.7]);
set(gcf, 'renderer', 'painters');

subplot('Position',[0.15 0.77 0.36 0.18]);

plot(-wq,Zf,'b')
hold on
plot(wq2,Zf,'r')
plot(wq3,Zf,'m')
plot(wq4,Zf,'g')
plot(wq5,Zf,'m')
plot(wq6,Zf,'k')
plot(wq7,Zf,'c')
l = legend('$$\overline{w \prime q_v \prime}$$', ...
    '$$\overline{w \prime q_s \prime}$$ as func of $$\overline{w \prime \theta \prime}$$ full' ...
    , '$$\overline{w \prime q_s \prime}$$ as func of $$\overline{w \prime \theta \prime}$$' ...
    , '$$\overline{w \prime q_s \prime}$$ as func of $$\overline{w \prime \theta \prime}$$ reduced' )
set(l,'interpreter','latex')

subplot('Position',[0.6 0.77 0.36 0.18]);
plot(thw,Zf)
hold on
plot (wth, Zf,'r')
%plot(wth6,Zf,'g')
subplot('Position',[0.15 0.5 0.36 0.18]);

plot (wq_t, Zf)
hold on
%plot (wthv2,Zf,'r')

