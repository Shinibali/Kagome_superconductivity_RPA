clearvars; close all;
fin = './input_files/'; addpath('./functions/'); img_dir = './images/';

modify_hr();

fntsz =25; nx = 99; totkf_in_ring = 6*25; % multiple of 6 per hexagonal ring
orbitals = 1; Ef =0.08; % top van Hove
% orbitals = 1; Ef =-2.02; % bottom van Hove

% orbitals = 2; Ef =2.19; % top most van Hove 
% orbitals = 2; Ef =0.16; % top second van Hove 
% orbitals = 2; Ef =-0.03; % top third van Hove 

str_orb= [num2str(orbitals),'orb'];
kT = 8.613e-5 * 11.605; bzedge = 1; hex='hexagon'; 

fid = fopen([fin,'unit_cell_info.txt']);
lines = textscan(fid,'%s','delimiter','\n'); fclose(fid);
lines = lines{1}; nol = max(size(lines));

for i = 1:nol
    temp = strsplit(strtrim(lines{i}),' ');
    if cellfun(@isempty,temp) ~= 1
        numbers = strsplit( erase( erase( temp{3},')' ),'(' ),',' ) ;
        if strcmp('a1',temp{1}) == 1
            a1 = [str2double(numbers{1}), str2double(numbers{2})];
        elseif strcmp('a2',temp{1}) == 1
            a2 = [str2double(numbers{1}), str2double(numbers{2})];
        elseif strcmp('TA',temp{1}) == 1
            TA = [str2double(numbers{1}), str2double(numbers{2})];
        elseif strcmp('TB',temp{1}) == 1
            TB = [str2double(numbers{1}), str2double(numbers{2})];
        elseif strcmp('TC',temp{1}) == 1
            TC = [str2double(numbers{1}), str2double(numbers{2})];
        end
    end
end

a1a2 = [a1;a2] ; b1b2 = 2. * pi * eye(2,2)/ transpose(a1a2) ;
TABC = [TA;TB;TC];

hexagon_vert = [ 1/3 2/3;
    2/3 1/3;
    1/3 -1/3;
    -1/3 -2/3;
    -2/3 -1/3;
    -1/3 1/3;
    1/3 2/3] ;
HSP = [0 0;
    1/3 2/3;
    1/2 1/2];
fout = 'a1a2_b1b2.mat'; save([fin,fout],'a1a2','b1b2','TABC','hexagon_vert','HSP','-mat');

kx = (0 : bzedge/nx : bzedge); kx=kx(1:end-1) ; ky=kx;
[KX,KY]=meshgrid(kx,ky) ; [nwan,Ek,~,HK] = eigenThomale(KX(:),KY(:),orbitals);
Ek = reshape(Ek,nwan,nx,nx);

Ek=Ek-Ef; nfill = 1./( 1 + exp((Ek)./kT) ) ;
nfill = sum( sum( sum( nfill,3 ), 2), 1) ./nx^2 ;

[b1b2klist,bandind] = kpcontour(orbitals,Ef,totkf_in_ring,b1b2,hexagon_vert) ;

totk = size(b1b2klist,1);
[~,~,ak] = eigenThomale(b1b2klist(:,1),b1b2klist(:,2),orbitals); weight = zeros(nwan,totk);
for j = 1:totk
    weight(:,j) = ak(:,bandind(j),j) .* conj( ak(:,bandind(j),j) ) ;
end
for i = 1:nwan
    [xlo, klist, weight(i,:), hex_vert] = ...
    reproduce_hexagon(b1b2klist, weight(i,:) ,b1b2, hexagon_vert) ;
end
color = zeros(size(weight,2),3);
if orbitals == 2
    for j = 1:size(weight,2)
        color(j,:) = (sum(weight([1,2],j))).* [0 1 0] + ...
            (sum(weight([3,4],j))).* [1 0 0] + ...
            (sum(weight([5,6],j))).* [0 0 1] ;
    end
else
    for j = 1:size(weight,2)
        color(j,:) = (sum(weight(1,j))).* [0 1 0] + ...
            (sum(weight(2,j))).* [1 0 0] + ...
            (sum(weight(3,j))).* [0 0 1] ;
    end
end

plotbands(50,fntsz,Ef,orbitals);

figure('units','normalized','outerposition',[0 0 1 1]); hold on; box on;
lbl = xlabel('$k_x$','Interpreter','Latex');
lbl.Position(1) = 0; % change horizontal position of xlabel
lbl.Position(2) = -1.5*pi; % change vertical position of xlabel
lbl =  ylabel('$k_y$','Interpreter','Latex');
lbl.Position(1) = -1.5*pi; % change horizontal position of ylabel
lbl.Position(2) = 0; % change vertical position of ylabel

mksiz = 6;
for j = 1:size(weight,2)
    plot( klist(j,1),klist(j,2),'o','markerfacecolor',color(j,:),...
        'markeredgecolor','none', 'MarkerSize',mksiz);
end

% Creating my own colormap
norows = 64; [myMap] = custom_colormap(norows);
colormap(myMap); cb= colorbar;
set(cb, 'ticks', [0,0.5,1], 'ticklabels', {'s_A','s_B','s_C'},'Fontsize',fntsz );

for iv = 2:size(hex_vert,1)
    start = hex_vert(iv-1,:); stop = hex_vert(iv,:);
    line = [start;stop]; plot(line(:,1),line(:,2),'-k','Linewidth',2);
end

HSP = HSP * b1b2;
HStext = {'$\Gamma$','$K$','$M$'};
text(HSP(1,1)-0.4,HSP(1,2),HStext{1},'Fontsize',fntsz,'Interpreter','Latex')
text(HSP(2,1)-0.9,HSP(2,2)+0.2,HStext{2},'Fontsize',fntsz,'Interpreter','Latex')
text(HSP(3,1)+0.2,HSP(3,2)+0.2,HStext{3},'Fontsize',fntsz,'Interpreter','Latex')

txt_pos = [-0.1,-0.25]*b1b2;
text(txt_pos(1),txt_pos(2),['$\mu_0=',num2str(Ef),'$'],'Fontsize',fntsz,'Interpreter','Latex')

scatter(HSP(:,1),HSP(:,2),50,'mo','Filled');

set(gca,'xticklabel',{'-\pi','\pi'},'Xtick',[-1*pi,1*pi])
set(gca,'yticklabel',{'-\pi','\pi'},'ytick',[-1*pi,1*pi])
set(gca,'xlim',[xlo -xlo],'ylim',[xlo -xlo],'Fontsize',fntsz);
axis square; hold off;
print([img_dir,'FS_',str_orb],'-dpng','-r300');
