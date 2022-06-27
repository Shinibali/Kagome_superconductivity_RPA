function [] = evec_kkp(orbitals,klist)

fin = './input_files/'; ftotk = 'totk_onFS.bin';
fpre = 'peig_re.bin'; fpim = 'peig_im.bin';
fmre = 'meig_re.bin'; fmim = 'meig_im.bin';

[~,~,peigjk] = eigenThomale( klist(:,1),  klist(:,2),orbitals);
[~,~,meigjk] = eigenThomale(-klist(:,1), -klist(:,2),orbitals);   

totk = size(klist,1);
fid = fopen([fin,ftotk],'w+','n'); fwrite(fid,totk,'double'); fclose(fid);

fid = fopen([fin,fpre],'w+','n'); fwrite(fid,real(peigjk),'double'); fclose(fid);
fid = fopen([fin,fpim],'w+','n'); fwrite(fid,imag(peigjk),'double'); fclose(fid);

fid = fopen([fin,fmre],'w+','n'); fwrite(fid,real(meigjk),'double'); fclose(fid);
fid = fopen([fin,fmim],'w+','n'); fwrite(fid,imag(meigjk),'double'); fclose(fid);

end