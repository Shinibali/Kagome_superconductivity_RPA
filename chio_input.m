%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name: write_text_in.m
% Shinibali, OCT 26,2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = chio_input(inputs)

kxsus = inputs(1); qxsus= inputs(2); orbitals = inputs(3); Ef = inputs(4); 
MatsuT = inputs(5); 

fin = './input_files/'; addpath('./functions/') ;

modify_hr(); extension = 1; % BZ extended scheme

fout = 'modified_hr.mat'; load([fin,fout],'nsite','Rvec','-mat');
fout = 'a1a2_b1b2.mat'; load([fin,fout],'a1a2','b1b2','-mat');
[Hamind, chiind] = basis_formation(nsite,orbitals) ;
Rvec = unique(Rvec,'rows'); nR = size(Rvec,1);
trans_Rvec = Rvec(:,1)*a1a2(1,:) + Rvec(:,2)*a1a2(2,:);

kx = (0 : extension/kxsus : extension) ; kx=kx(1:end-1) ; ky = kx; 
trans_kvec = zeros(2,kxsus,kxsus) ; [KX,KY]=meshgrid(kx,ky) ;
[nwan,Ek,ak] = eigenThomale(KX(:),KY(:),orbitals);
for ix = 1:kxsus
    for iy = 1:kxsus
        kxky = kx(ix)*b1b2(1,:) + ky(iy)*b1b2(2,:) ;
        trans_kvec(1:2,iy,ix) = kxky(:);
    end
end
Ek = reshape(Ek,nwan,kxsus,kxsus); ak = reshape(ak,nwan,nwan,kxsus,kxsus);
Ek = Ek-Ef ; kT = MatsuT * 8.617e-5; fermi = 1.0 ./ ( exp(Ek./kT) + 1.0 ) ;

qx = (0 : extension/qxsus : extension) ; qx=qx(1:end-1);
[qx,qy]=meshgrid(qx); regq=[qx(:),qy(:)];
[nx,ny] = find_index(regq(:,1),regq(:,2),kx,ky,b1b2);

fw = 'NsiteNorbNkxNqxNr.bin'; fnx = 'nx.bin'; fny = 'ny.bin'; fqp = 'qpoints.bin';
fek = 'ek.bin'; ffermi = 'fermi.bin'; fakre = 'akre.bin'; fakim = 'akim.bin';
fHi = 'Hamind.bin'; fIi = 'chiind.bin'; fRvec = 'Rvec.bin'; fkv = 'kvec.bin';

inputs = [nsite,orbitals,kxsus,qxsus,MatsuT, Ef,nR];
fid = fopen([fin,fw],'w+','n');fwrite(fid,inputs,'double'); fclose(fid);
fid = fopen([fin,fqp],'w+','n');fwrite(fid,regq,'double'); fclose(fid);
fid = fopen([fin,fnx],'w+','n');fwrite(fid,nx,'double'); fclose(fid);
fid = fopen([fin,fny],'w+','n');fwrite(fid,ny,'double'); fclose(fid);
fid = fopen([fin,fek],'w+','n');fwrite(fid,Ek,'double'); fclose(fid);
fid = fopen([fin,fRvec],'w+','n');fwrite(fid,trans_Rvec,'double'); fclose(fid);
fid = fopen([fin,fkv],'w+','n');fwrite(fid,trans_kvec,'double'); fclose(fid);
fid = fopen([fin,fHi],'w+','n');fwrite(fid,Hamind,'double'); fclose(fid);
fid = fopen([fin,fIi],'w+','n');fwrite(fid,chiind,'double'); fclose(fid);
fid = fopen([fin,ffermi],'w+','n');fwrite(fid,fermi,'double'); fclose(fid);
fid = fopen([fin,fakre],'w+','n');fwrite(fid,real(ak),'double'); fclose(fid);
fid = fopen([fin,fakim],'w+','n');fwrite(fid,imag(ak),'double'); fclose(fid); 

end