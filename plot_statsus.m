clearvars; close all

fres = './bareresults/'; fin = './input_files/';
img_dir = './images/'; fntsz = 15; addpath('./functions/')
U = 2.45; Vnn = 0.; J = U/10; mksz=10; hex='hexagon'; zone = 'first'; 

fwkq = 'NsiteNorbNkxNqxNr.bin'; fid = fopen([fin,fwkq],'r','n');
ins = fread(fid,'double'); fclose(fid);
nsite = ins(1); norb = ins(2); nqgrid = ins(4); nR = ins(7); totq = nqgrid^2;
fqp = 'qpoints.bin'; fid = fopen([fin,fqp],'r','n');
regq = fread(fid,[totq 2],'double'); fclose(fid);

fout = 'a1a2_b1b2.mat'; load([fin,fout],'b1b2','hexagon_vert','HSP','-mat');

[fullchio,lindhardt_chio] = copy_chi(nsite,norb);
[Intmat,Identity,R0ind] = nn_interaction(nsite,norb,U,J,Vnn);
lindhardt_chio = reshape( lindhardt_chio, (norb*nsite),(norb*nsite),(norb*nsite),(norb*nsite),[]);

totbare = zeros(totq,4); % 3 for charge, xx, yy and zz components
[~,~,Int_ind] = basis_formation(nsite,norb); dimen = numel(Int_ind);
[rows,cols,val_Vnn_fullq] = nn_interact_R0block(regq(:,1),regq(:,2),...
                                        nsite,norb,Vnn,Int_ind,R0ind);
for i = 1:totq
    int_nn_R0 = sparse(rows,cols,val_Vnn_fullq(:,i),dimen,dimen);
    Interaction = Intmat + int_nn_R0 ;
    Achi(:,:) = fullchio(:,:,i);
%     totbare(i,1) = lindhardt_chio(1,1,1,1,i) + lindhardt_chio(1,1,2,2,i) + ...
%         lindhardt_chio(1,1,3,3,i) + lindhardt_chio(2,2,1,1,i) + lindhardt_chio(2,2,2,2,i) + ...
%         lindhardt_chio(2,2,3,3,i) + lindhardt_chio(3,3,1,1,i) + lindhardt_chio(3,3,2,2,i) + ...
%         lindhardt_chio(3,3,3,3,i) ;
% [~,val] = eig( [ lindhardt_chio(1,1,1,1,i) , lindhardt_chio(1,1,2,2,i) , ...
%         lindhardt_chio(1,1,3,3,i) ; lindhardt_chio(2,2,1,1,i),  lindhardt_chio(2,2,2,2,i) , ...
%         lindhardt_chio(2,2,3,3,i) ; lindhardt_chio(3,3,1,1,i), lindhardt_chio(3,3,2,2,i) , ...
%         lindhardt_chio(3,3,3,3,i)] );
%     totbare(i,2)=max(real(diag(val)));
    chirpa = Achi/( Identity - Interaction * Achi );
    [chi_xyz] = physical_chi(nsite,norb,R0ind,Int_ind,chirpa) ;
    totbare(i,1:4) = chi_xyz([1,2,3,4]);
end
totbare = real(totbare);

figure('units','normalized','outerposition',[0 0 1 1]);
for f = 1:4
    if f==1
        tit = '$\chi^{CDO}$';
    elseif f==2
        tit = '$\chi^{CBO}$';
    elseif f==3
        tit = '$\chi^{xx}$';
    else
        tit = '$\chi^{zz}$';
    end
    tb = totbare(:,f) ;
    [xlo, full_q, full_qty, hex_vert] = ...
        reproduce_hexagon(regq,tb,b1b2, hexagon_vert);
    
    subplot (2,2,f); box on; hold on
    scatter(full_q(:,1),full_q(:,2),mksz,full_qty,'Filled');
    xlabel('$q_x$','Interpreter','latex') ;
    
    for iv = 2:size(hex_vert,1)
        start = hex_vert(iv-1,:); stop = hex_vert(iv,:);
        line = [start;stop]; plot(line(:,1),line(:,2),'-k','Linewidth',2);
    end
    
    if f==1 || f==3
        ylabel('$q_y$','Interpreter','latex');
        set(gca,'yticklabel',{'-\pi','0','\pi'},'ytick',[-pi,0,pi])
    else
        set(gca,'yticklabel',[],'ytick',[])
    end
    set(gca,'xticklabel',{'-\pi','0','\pi'},'Xtick',[-pi,0,pi])
    set(gca,'xlim',[xlo -xlo],'ylim',[xlo -xlo]);
    title(tit,'interpreter','Latex')
    colormap parula; colorbar; shading interp; axis square;
    set(gca,'Fontsize',fntsz); hold off
end
print([img_dir,'susceptibility'],'-dpng','-r300');