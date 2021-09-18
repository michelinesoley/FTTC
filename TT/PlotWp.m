function [t] = PlotWp(G_tt,Gc_tt,nx,d,r1,r2,k,t)

figure

str = {strcat('Time = ' , num2str(t(k)*2.418E-2),' fs')};

% Visualize wavepacket
if (d == 2)
    hold on;
    v1=reshape(abs(full(G_tt(:,:))).^2,[nx nx]);
    contour(r1,r2,v1','linewidth', 2,'linestyle',':','color','k');
    v3=reshape(abs(full(Gc_tt(:,:))).^2,[nx nx]);
    contour(r1,r2,v3','linewidth', 2,'linestyle','--','color','b');
    hold off;
    grid on;
    xlabel('x_1 [au]','fontsize',18);
    ylabel('x_2 [au]','fontsize',18);
elseif(d == 3)
    hold on;
    v1=reshape(abs(full(G_tt(:,:,nx/2))).^2,[nx nx]);
    contour(r1,r2,v1','linewidth', 2,'linestyle',':','color','k');
    v3=reshape(abs(full(Gc_tt(:,:,nx/2))).^2,[nx nx]);
    contour(r1,r2,v3','linewidth', 2,'linestyle','--','color','b');
    hold off;
    grid on;
    xlabel('x_1 [au]','fontsize',18);
    ylabel('x_2 [au]','fontsize',18);
elseif(d == 4)
    hold on;
    v1=reshape(abs(full(G_tt(:,:,nx/2,nx/2))).^2,[nx nx]);
    contour(r1,r2,v1','linewidth', 2,'linestyle',':','color','k');
    v3=reshape(abs(full(Gc_tt(:,:,nx/2,nx/2))).^2,[nx nx]);
    contour(r1,r2,v3','linewidth', 2,'linestyle','--','color','b');
    hold off;
    grid on;
    xlabel('x_1 [au]','fontsize',18);
    ylabel('x_2 [au]','fontsize',18);
elseif(d == 6)
    hold on;
    v1=reshape(abs(full(G_tt(:,:,nx/2,nx/2,nx/2,nx/2))).^2,[nx nx]);
    contour(r1,r2,v1','linewidth', 2,'linestyle',':','color','k');
    v3=reshape(abs(full(Gc_tt(:,:,nx/2,nx/2,nx/2,nx/2))).^2,[nx nx]);
    contour(r1,r2,v3','linewidth', 2,'linestyle','--','color','b');
    hold off;
    grid on;
    xlabel('x_1 [au]','fontsize',18);
    ylabel('x_2 [au]','fontsize',18);
elseif(d == 12)
    hold on;
    v1=reshape(abs(full(G_tt(:,:,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2))).^2,[nx nx]);
    contour(r1,r2,v1','linewidth', 2,'linestyle',':','color','k');
    v3=reshape(abs(full(Gc_tt(:,:,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2))).^2,[nx nx]);
    contour(r1,r2,v3','linewidth', 2,'linestyle','--','color','b');
    hold off;
    grid on;
    xlabel('x_1 [au]','fontsize',18);
    ylabel('x_2 [au]','fontsize',18);
elseif(d == 36)
    hold on;
    v1=reshape(abs(full(G_tt(:,:,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2))).^2,[nx nx]);
    contour(r1,r2,v1','linewidth', 2,'linestyle',':','color','k');
    v3=reshape(abs(full(Gc_tt(:,:,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2))).^2,[nx nx]);
    contour(r1,r2,v3','linewidth', 2,'linestyle','--','color','b');
    hold off;
    grid on;
    xlabel('x_1 [au]','fontsize',18);
    ylabel('x_2 [au]','fontsize',18);
elseif(d == 50)
    hold on;
    v1=reshape(abs(full(G_tt(:,:,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2))).^2,[nx nx]);
    contour(r1,r2,v1','linewidth', 2,'linestyle',':','color','k');
    v3=reshape(abs(full(Gc_tt(:,:,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,...
        nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2,nx/2))).^2,[nx nx]);
    contour(r1,r2,v3','linewidth', 2,'linestyle','--','color','b');
    hold off;
    grid on;
    xlabel('x_1 [au]','fontsize',18);
    ylabel('x_2 [au]','fontsize',18);
end
title(str,'fontsize',18);
legend('SOFT','Chebyshev')
set(gca,'FontSize',15);
end

