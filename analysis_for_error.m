%Script for calculating and plotting the SEM DIC Challenge data versus the
%output from FIDIC

close all

%Find mean and sd across x2 for each x1 point (as Blaber 2015)
u1_11m = mean(u{1}{1});
u1_11sd = std(u{1}{1});
px_dist1 = 64:8:2049-64;
px_dist2 = 1:8:2049;

u1_11m_n = mean(u_norm{1}{1});
u1_11sd_n = std(u_norm{1}{1});

%Desample to plot errorbars less densely
u1_desamp_m = u1_11m(1:8:end);
u1_desamp_sd =  u1_11sd(1:8:end);
desamp_px1 = px_dist1(1:8:end);
desamp_px2 = px_dist2(1:8:end);
u1_desamp_m_n = u1_11m_n(1:8:end);
u1_desamp_sd_n =  u1_11sd_n(1:8:end);

%Build the plots
figure
plot(L5,'-r','linewidth',2)
hold on
plot(px_dist1,u1_11m_n,':b','linewidth',2)
errorbar(desamp_px1,u1_desamp_m_n,u1_desamp_sd_n,'o');
hold on
plot(px_dist2,-u1_11m,'--m','linewidth',2)
errorbar(desamp_px2,-u1_desamp_m,u1_desamp_sd,'o');
xlabel('x, pixels')
ylabel('u, pixels')
legend('Commanded Displacement','FIDIC-Normalized Result',...
    'Errorbars-Normalizated','FIDIC Result','Errorbars')
title('SEM Challenge Sample14 L3')
set(gca,'FontSize',18,'FontWeight','bold');



% u2_11m = -mean(u{2}{1});
% u2_11sd = std(u{2}{1});
% u2_desamp_m = u2_11m(1:8:end);
% u2_desamp_sd =  u2_11sd(1:8:end);

% figure
% plot(L3,'-r','linewidth',2)
% hold on
% plot(px_dist,u2_11m,'b','linewidth',1)
% errorbar(desamp_px,u2_desamp_m,u2_desamp_sd,'o');
% xlabel('x, pixels')
% ylabel('u, pixels')
% legend('Commanded Displacement','FIDIC Result','Example Errorbars')
% title('SEM Challenge Sample14 L3')
% set(gca,'FontSize',18,'FontWeight','bold');
