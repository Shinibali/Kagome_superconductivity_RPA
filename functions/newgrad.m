function [dEdx,dEdy,gradE] = newgrad(kstep,kx,ky,orbitals,bandind,b1b2)

nx = max(size(kx));
dEdx = zeros(nx,1); dEdy = dEdx;
dx = [kstep,0]/b1b2 ; dy = [0,kstep]/b1b2 ;

[~,Ekpdx] = eigenThomale(kx+dx(1),ky+dx(2),orbitals);
[~,Ekmdx] = eigenThomale(kx-dx(1),ky-dx(2),orbitals);

[~,Ekpdy] = eigenThomale(kx+dy(1),ky+dy(2),orbitals);
[~,Ekmdy] = eigenThomale(kx-dy(1),ky-dy(2),orbitals);

kpdx = ( kx+dx(1) ) * b1b2(1,:) + ( ky+dx(2) ) * b1b2(2,:);
kmdx = ( kx-dx(1) ) * b1b2(1,:) + ( ky-dx(2) ) * b1b2(2,:);
kpdy = ( kx+dy(1) ) * b1b2(1,:) + ( ky+dy(2) ) * b1b2(2,:);
kmdy = ( kx-dy(1) ) * b1b2(1,:) + ( ky-dy(2) ) * b1b2(2,:);

for i = 1:nx
    ind = bandind(i);
    dEdx(i) = Ekpdx(ind,i) - Ekmdx(ind,i);
    dEdy(i) = Ekpdy(ind,i) - Ekmdy(ind,i);
end
len_dx = sqrt( sum( (kpdx-kmdx).^2 ,2 ) ) ; len_dy = sqrt( sum( (kpdy-kmdy).^2 ,2 ) )  ;
dEdx = dEdx./len_dx ; dEdy = dEdy./len_dy ;
gradE = sqrt( dEdx.^2 + dEdy.^2 );

end