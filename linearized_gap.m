%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for linearized gap evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = linearized_gap(U,J,Vnn)

fres = './bareresults/'; fin = './input_files/'; addpath('./functions')

fw = 'NsiteNorbNkxNqxNr.bin'; fid = fopen([fin,fw],'r','n'); 
ins = fread(fid,'double'); fclose(fid); orbitals = ins(2);

fn = 'scatvec.mat'; load([fres,fn],'rearrklist','rearrbdin','len_element','-mat');
fout = 'a1a2_b1b2.mat'; load([fin,fout],'b1b2','-mat');
totk = max(size(rearrklist)); 

kstep = 1/99; % k-step size for finding gradient in Ek
[~,~,grad] = newgrad(kstep,rearrklist(:,1),rearrklist(:,2),orbitals,rearrbdin,b1b2 );

noU = numel(U); GammaU = zeros(totk,totk,4,4,noU); 
for uu = 1:noU
    fgre = ['re_gamma_U',num2str(uu),'.bin']; fid = fopen([fres,fgre],'r','n');
    re = fread(fid,'double'); fclose(fid); delete([fres,fgre]);
    re = reshape(re,totk,totk,4,4);
    
    fgim = ['im_gamma_U',num2str(uu),'.bin']; fid = fopen([fres,fgim],'r','n');
    im = fread(fid,'double'); fclose(fid); delete([fres,fgim]);
    im = reshape(im,totk,totk,4,4);
    
    gn1n2 = complex(re,im);
    for kp = 1:totk
        gn1n2(:,kp,:,:) = gn1n2(:,kp,:,:) * (len_element(kp)/grad(kp)) ;
    end
    GammaU(:,:,:,:,uu) = gn1n2; 
end

fn = 'gapstructure.mat'; save([fres,fn],'U','J','Vnn','GammaU','-mat');

end