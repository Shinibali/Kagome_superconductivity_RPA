%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shinibali, NOV 4,2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clearvars

fres = './bareresults/'; addpath('./functions'); fin = './input_files/';
img_dir = './images/';

fntsz =15; mksz=15; hex='hexagon'; zone = 'first'; result_type = 'bareU_';

fout = 'a1a2_b1b2.mat'; load([fin,fout],'b1b2','hexagon_vert','HSP','-mat');
area = norm(cross([b1b2(1,:),0],[b1b2(2,:),0]));

fn = 'scatvec.mat'; load([fres,fn],'rearrklist','rearrbdin','-mat');
totk = max(size(rearrklist));

fn = [result_type,'gapstructure.mat']; load([fres,fn],'U','J','Vnn','GammaU','-mat');
GammaU = -GammaU./area ;

noU = numel(U); SCgk = zeros(totk,4,4,noU); lambdaSC = zeros(4,4,noU);

for uu = 1:noU
    for l = 1:4
        [vector,value] = eig( GammaU(:,:,l,l,uu));
        [eigval,ind] = sort( real(diag(value)) );
        lambdaSC(:,l,uu) = eigval(end:-1:end-3);
        SCgk(:,1:4,l,uu) = vector(:,ind(end:-1:end-3)) ;
    end
end

lind = {'00','xx','yy','zz'};
figure('units','normalized','outerposition',[0 0 1 1]); box on; hold on;
end_id = noU-1;
for l = 1:4
    y(:)=lambdaSC(1,l,:);
    plot(U,y,'-o','Displayname',['$\lambda_{',lind{l},'}$'],'Linewidth',2);
end
lgd = legend('orientation','vertical','location','northwest');
set(lgd,'interpreter','Latex'); legend boxoff;
xlabel('$U$','Interpreter','latex') ;
ylabel('$\lambda$','Interpreter','latex');
axis tight; axis square;
% set(gca,'xlim',[0 U(end)],'ylim',[0 2.5]);
set(gca,'Fontsize',fntsz*1.5); hold off
print([img_dir,result_type,'lambda_vs_U'],'-dpng','-r300');

for uu = 1:noU
    figure('units','normalized','outerposition',[0 0 1 1]);
    for plt = 1:8
        subplot(2,4,plt) ; box on; hold on;
        if plt <= 4
            lltext=lind{1}; qty = real(SCgk(:,plt,1,uu)) ;
            lambda = ['$\lambda_',num2str(plt-1),'=',num2str(round(lambdaSC(plt,1,uu),6)),'$'] ;
        else
            lltext=lind{2}; qty = real(SCgk(:,plt-4,2,uu)) ;
            lambda = ['$\lambda_',num2str(plt-1-4),'=',num2str(round(lambdaSC(plt-4,2,uu),6)),'$'] ;
        end
        title( ['$ll=',lltext,', U = ',num2str(round(U(uu),2)),'\, , V = ',...
            num2str(round(Vnn(uu),2)),'$'],'interpreter','Latex')
        text( -1.3*pi,1.3*pi,lambda,'interpreter','Latex','Fontsize',fntsz)
        [xlo, full_q, full_qty, hex_vert] = ...
            reproduce_hexagon(rearrklist,qty,b1b2, hexagon_vert);

        scatter(full_q(:,1),full_q(:,2),mksz,full_qty,'Filled');
        if plt == 1 || plt == 5
            ylabel('$k_y$','Interpreter','latex');
            set(gca,'yticklabel',{'-\pi','0','\pi'},'ytick',[-pi,0,pi])
        else
            set(gca,'ytick',[])
        end
        if plt > 4
            xlabel('$k_x$','Interpreter','latex') ;
            set(gca,'xticklabel',{'-\pi','0','\pi'},'xtick',[-pi,0,pi]);
        else
            set(gca,'xtick',[])
        end

        for iv = 2:size(hex_vert,1)
            start = hex_vert(iv-1,:); stop = hex_vert(iv,:);
            line = [start;stop];
            plot(line(:,1),line(:,2),'-','Linewidth',1.5,'color',[17,17,17,255/3]/255);
        end

        set(gca,'xlim',[xlo -xlo],'ylim',[xlo -xlo]);
        caxis([-max(abs(qty)) max(abs(qty)) ]./2)
        colormap(redblue); colorbar; shading interp; axis square;
        set(gca,'Fontsize',fntsz); hold off
    end
    print([img_dir,[result_type,'gap_U',num2str(uu)]],'-dpng','-r300');
end
