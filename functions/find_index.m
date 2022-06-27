function [nx,ny] = find_index(qx,qy,kx,ky,b1b2)

BZedge = kx(end) + kx(2)-kx(1); xup = kx(end); etol = 1e-2;
[X,Y] = meshgrid(kx,ky) ; klist = [X(:), Y(:)] ;

qxsus = max(size(qx)) ; kxsus = max(size(kx));

nx = zeros(kxsus,qxsus); ny = nx; kq_k_ind = zeros(kxsus,2);
for iq = 1:qxsus
    kq_k_ind(:,1) = mod( kx(:)+qx(iq) , BZedge);
    kq_k_ind(:,2) = mod( ky(:)+qy(iq) , BZedge);
    for ik = 1:kxsus
        [~,nx(ik,iq)] = min( abs( kx - kq_k_ind(ik,1) )) ;
        [~,ny(ik,iq)] = min( abs( ky - kq_k_ind(ik,2) )) ;
    end
end

% nx = zeros(kxsus^2,qxsus); ny = nx;
% klist = klist(:,1) * b1b2(1,:) + klist(:,2) * b1b2(2,:) ; 
% q = qx(:) * b1b2(1,:) + qy(:) * b1b2(2,:) ; qx = q(:,1); qy = q(:,2);
% 
% for iq = 1:qxsus
%     kq = [klist(:,1)+qx(iq),klist(:,2)+qy(iq)];
%     kq_k_ind = kq/b1b2 ;
%     
%     logic = (kq_k_ind(:,1)-xup) > etol ; 
%     kq_k_ind(:,1) = kq_k_ind(:,1)-BZedge*logic;
%     logic = (kq_k_ind(:,2)-xup) > etol ; 
%     kq_k_ind(:,2) = kq_k_ind(:,2)-BZedge*logic;
%     for ik = 1:kxsus^2
%         [~,nx(ik,iq)] = min( abs( kx - kq_k_ind(ik,1) )) ;
%         [~,ny(ik,iq)] = min( abs( kx - kq_k_ind(ik,2) )) ;
%     end
% end

end