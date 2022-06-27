%###########################################
% Function: getH
%###########################################
function [nwan,Ek,ak,HK] = eigenThomale(kx,ky,orbitals)

a1a2 = [1.0,  0.0; 
        1/2,  sqrt(3)/2];

b1b2 = 2. * pi * eye(2,2)/ transpose(a1a2) ;
              
txz = 1; tyz = 0.5; t0 = 0.002 ; Exz = 2.182 ; Eyz = 0.055;
a123 = [1/2,   0;
        1/4,  sqrt(3)/4;
       -1/4,  sqrt(3)/4];

kxy = kx * b1b2(1,:) + ky * b1b2(2,:) ; nkxy = size(kx,1) ; 

phiAB = @(x,y) 1 + exp( -2i *( x*a123(1,1) + y*a123(1,2) ) ) ;
phiAC = @(x,y) 1 + exp( -2i *( x*a123(2,1) + y*a123(2,2) ) ) ;
phiBC = @(x,y) 1 + exp( -2i *( x*a123(3,1) + y*a123(3,2) ) ) ;

if orbitals == 1
    nwan = 3 ;
    hamil = @(x,y) -[Exz                     txz*phiAB(x,y)        txz*phiAC(x,y);
                     txz*conj(phiAB(x,y))    Exz                   txz*phiBC(x,y);
                     txz*conj(phiAC(x,y))    txz*conj(phiBC(x,y))    Exz];
    chempot = -Exz;
else
    nwan = 6 ;
    hamil = @(x,y) -[-Exz  0                          txz*phiAB(x,y)  t0*phiAB(x,y)      txz*phiAC(x,y)  t0*phiAC(x,y);
                     0    Eyz                         t0*phiAB(x,y)  tyz*phiAB(x,y)      t0*phiAC(x,y)  tyz*phiAC(x,y);
                     txz*conj(phiAB(x,y))  t0*conj(phiAB(x,y))   -Exz    0                          txz*phiBC(x,y) t0*phiBC(x,y);
                     t0*conj(phiAB(x,y))   tyz*conj(phiAB(x,y))   0    Eyz                           t0*phiBC(x,y) tyz*phiBC(x,y);
                     txz*conj(phiAC(x,y))  t0*conj(phiAC(x,y))    txz*conj(phiBC(x,y)) t0*conj(phiBC(x,y))     -Exz   0;
                     t0*conj(phiAC(x,y))   tyz*conj(phiAC(x,y))   t0*conj(phiBC(x,y))  tyz*conj(phiBC(x,y))     0   Eyz];
    chempot = 0;%-Eyz;
end
Ek = zeros(nwan,nkxy);  ak = zeros(nwan,nwan,nkxy); HK = zeros(nwan,nwan,nkxy);
Efmatrix =  chempot * eye(nwan,nwan);
for ik = 1:nkxy
    xx = kxy(ik,1); yy = kxy(ik,2);
    HK(:,:,ik) = hamil(xx,yy) - Efmatrix ;    
    [vec,val] = eig(HK(:,:,ik)); val=real(diag(val));
    Ek(:,ik) = val; ak(:,:,ik) = vec ;
end

end