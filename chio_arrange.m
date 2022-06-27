%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns the two orbital bare chi(p.q.s.t) for a specific q point
% Shinibali, NOV 4,2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = chio_arrange(qxsus)

fres = './bareresults/';  addpath('./functions/')
fout = 'baresus.mat'; fchiore = 'chio_re.bin'; fchioim = 'chio_im.bin';  

fid = fopen([fres,fchiore],'r','n'); chio_re = fread(fid,'double'); fclose(fid); delete([fres,fchiore]);
fid = fopen([fres,fchioim],'r','n'); chio_im = fread(fid,'double'); fclose(fid); delete([fres,fchioim]); 
chio = complex(chio_re,chio_im); clear chio_re chio_im
totq = qxsus^2; chio = reshape(chio,[],totq); save([fres,fout],'chio','-mat');

end