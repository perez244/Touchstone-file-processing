
% 
% % % plot S21
% figure('Renderer', 'painters', 'Position', [200 300 1200 500])
% plot(F,S21_50,F,S21_70,F,S21_100,F,S21_250,F,S21_365,F,S21_500,F,S21_750,F,S21_1000,'--r','LineWidth',1.5)
% set(gca, 'FontName', 'Consolas')
% xlabel('Frequency(GHz)','FontSize',11,'FontWeight','bold')
% ylabel('Transmission(dB)','FontSize',11,'FontWeight','bold')
% title('S21 for varying Cu_2O thickness ')
% grid on
% lgd2 = legend('50nm','70nm','100nm','250nm','365nm','500nm','750nm','1000nm');
% lgd2.FontSize = 14;

% figure('Renderer', 'painters', 'Position', [300 300 1200 500])
% plot(F,S11_500,F,newS11_500,'--r','LineWidth',1.5)
% set(gca, 'FontName', 'Consolas')
% xlabel('Frequency(GHz)','FontSize',11,'FontWeight','bold')
% ylabel('Reflection(dB)','FontSize',11,'FontWeight','bold')
% title('S11 for varying Cu_2O thickness ')
% grid on
% lgd = legend('500nm','new 500nm');
% lgd.FontSize = 14;

% % % plot S21
% %figure('Renderer', 'painters', 'Position', [200 300 1200 500])
% plot(freq,total_att_constant_dB,freq,quarterUM.*Ksr,freq,OneUM.*Ksr,freq,halfUM.*Ksr,freq,ThreeQuarterUM.*Ksr,'LineWidth',1.5)
% set(gca, 'FontName', 'Consolas')
% grid on
% xlabel('Frequency(GHz)','FontSize',11,'FontWeight','bold')
% ylabel('Attenuation (dB)','FontSize',11,'FontWeight','bold')
% title('Simulation of a 350nm Cu_2O/CuO oxide layer')
% grid on
% lgd2 =legend('"Atomicly smooth"','.25 \mum','.5\mum','.75\mum','1 \mum ');
% lgd2.FontSize = 14;



% plot(F,CorrectedMetalContribution,':',F,CorrectedDielecricContributiondB,':',F,TotalS21for350nmCuO_Cu2OmixCuSiO2Si,':','LineWidth',2.5)
% set(gca, 'FontName', 'Consolas')
% xlabel('Frequency(GHz)','FontSize',11,'FontWeight','bold')
% ylabel('Reflection(dB)','FontSize',11,'FontWeight','bold')
% title(' New S21 signal ')
% grid on
% lgd = legend('Metal','Dielectric','Total');
% lgd.FontSize = 14;

hold all
%plot figure S21
figure(1)
plot(F,CuOmaxS21,'--b',F,CuOminS21,'--r','LineWidth',1.5)
set(gca, 'FontName', 'Consolas')
xlabel('Frequency(GHz)','FontSize',11,'FontWeight','bold')
ylabel('Insertion loss(dB)','FontSize',11,'FontWeight','bold')
title('CuO layer on TSV with Leti \sigma values')
grid on
lgd = legend('Si \sigma = 10^4 S/m','Si \sigma = 6.667 S/m');
lgd.FontSize = 14;

figure(2)
plot(F,Cu2OmaxS21,'--b',F,Cu2OminS21,'--r','LineWidth',1.5)
set(gca, 'FontName', 'Consolas')
xlabel('Frequency(GHz)','FontSize',11,'FontWeight','bold')
ylabel('Insertion loss(dB)','FontSize',11,'FontWeight','bold')
title('Cu_2O layer on TSV with Leti \sigma values')
grid on
lgd = legend('Si \sigma = 10^4 S/m','Si \sigma = 6.667 S/m');
lgd.FontSize = 14;

figure(3)
plot(F,Cu2OmaxS21,'-.b*',F,CuOmaxS21,'--c',F,Cu2OminS21,'--m',F,CuOminS21,':rs',F,zeroCu2OS21,':g','LineWidth',1.5)
set(gca, 'FontName', 'Consolas')
xlabel('Frequency(GHz)','FontSize',11,'FontWeight','bold')
ylabel('Insertion loss(dB)','FontSize',11,'FontWeight','bold')
title('Coppe Oxide layers on TSV with Leti \sigma values')
grid on
lgd = legend('Cu_2O with Si \sigma = 10^4 S/m','CuO with Si \sigma = 10^4 S/m','Cu_2O with Si \sigma = 6.667 S/m','CuO withSi \sigma = 6.667 S/m','Cu_2O with Si \sigma = 0 S/m');
lgd.FontSize = 14;

figure(4)
plot(F,zeroCu2OS21,'--b','LineWidth',1.5)
set(gca, 'FontName', 'Consolas')
xlabel('Frequency(GHz)','FontSize',11,'FontWeight','bold')
ylabel('Insertion loss(dB)','FontSize',11,'FontWeight','bold')
title('Cu_2O layer on TSV with\sigma = 0')
grid on
lgd = legend('Cu_2O');
lgd.FontSize = 14;