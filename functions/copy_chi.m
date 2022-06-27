function [fullchio,lindhardt_chio] = copy_chi(nsite,norb)

fin = './input_files/'; fwkq = 'NsiteNorbNkxNqxNr.bin';
fid = fopen([fin,fwkq],'r','n'); ins = fread(fid,'double'); fclose(fid);
nqgrid = ins(4); nR = ins(7); totq = nqgrid^2; dimen = norb*nsite;

fres = './bareresults/';fout = 'baresus.mat'; load([fres,fout],'chio','-mat');
chio = reshape( chio,dimen,dimen,dimen,dimen,nR,nR,[]) ;

fin = './input_files/'; fout = 'modified_hr.mat'; load([fin,fout],'Rvec','-mat');
[uniR] = unique(Rvec,'rows');
R0ind = find( sqrt( uniR(:,1).^2 + uniR(:,2).^2 ) == 0 );
lindhardt_chio = reshape( chio(:,:,:,:,R0ind,R0ind,:), ...
    dimen,dimen,dimen,dimen,nqgrid,nqgrid);

% [~,~,~,short_Int_ind] = basis_formation(nsite,norb);
% d1 = (norb*nsite*2)^2 * nR ; fullchio = zeros(d1,d1,totq);
% for iq = 1:totq
%     for ir2 = 1:nR
%         for ir1 = 1:nR
%             for so4 = 1:nsite
%                 for so3 = 1:nsite
%                     for so2 = 1:nsite
%                         for so1 = 1:nsite
%                             block = chio(so1,so2,so3,so4,ir1,ir2,iq);
%                             for isp2 = 1:2
%                                 for isp1 = 1:2
%                                     row = short_Int_ind( so1,so2,isp1,isp2,ir1 );
%                                     col = short_Int_ind( so3,so4,isp2,isp1,ir2 );
%                                     fullchio(row,col,iq) = block;
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

d1 = (norb*nsite*2)^2 * nR ;
% % % chio = reshape( chio,dimen,dimen,dimen,dimen,nR,nR,totq) ;
chio = permute(chio,[1,2,5,3,4,6,7]);
fullchio = zeros(dimen,dimen,2,2,nR,dimen,dimen,2,2,nR,totq) ;
for isp2 = 1:2
    for isp1 = 1:2
        fullchio(:,:,isp1,isp2,:,:,:,isp2,isp1,:,:) = chio(:,:,:,:,:,:,:);
    end
end
fullchio = reshape(fullchio,d1,d1,totq);

end