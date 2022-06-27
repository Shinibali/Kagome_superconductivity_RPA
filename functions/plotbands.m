%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name: RPA function along high symmtery cut
% Shinibali, Dec 28,2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [color] = plotbands(nkp_seg,fntsz,Ef,orbitals)

[klist,~,~,index] = generate_kpath(nkp_seg);
str_xt={'$\Gamma$',"$M$","$K$",'$\Gamma$'};

nks = size(klist,1);
[nwan,Ek,ak] = eigenThomale(klist(:,1),klist(:,2),orbitals);
color = zeros(nwan,3,nks);
if orbitals == 2
    for kl = 1:nks
        for io = 1:nwan
            weight = ak(:,io,kl) .* conj( ak(:,io,kl) ) ;
            color(io,:,kl) = (sum(weight([1,2]))).* [0 1 0] + ...
                (sum(weight([3,4]))).* [1 0 0] + ...
                (sum(weight([5,6]))).* [0 0 1] ;
        end
    end
else
    for kl = 1:nks
        for io = 1:nwan
            weight = ak(:,io,kl) .* conj( ak(:,io,kl) ) ;
            color(io,:,kl) = (sum(weight(1))).* [0 1 0] + ...
                (sum(weight(2))).* [1 0 0] + ...
                (sum(weight(3))).* [0 0 1] ;
        end
    end
end
Ek = Ek-Ef;

x_lim = [1 nks-1]; y_lim = [min(min(Ek)) max(max(Ek))].*1.1;

figure; box on; hold on
set(gca,'YLim',y_lim,'XLim',x_lim,'XTick',[]);%,'YTick',[-1,0,1]);
plot([index';index'],y_lim','k-'); plot(x_lim,[0 0],'k-');
text(index,y_lim(ones(size(index)))-0.05*(y_lim(2)-y_lim(1)),...
    str_xt,'Fontsize',fntsz,'Interpreter','Latex');
ylabel('$E/|t|$','Interpreter','Latex')

% for kl = 1:nks
    for io = 1:nwan
        plot(1:kl ,Ek(io,:),'-b','markerfacecolor','b',...%color(io,:,kl),...
            'Linewidth',2,'markeredgecolor','none', 'MarkerSize',4);
    end
% end
set(gca,'Fontsize',fntsz ); hold off;

end