%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns the two orbital bare chi(p.q.s.t) for a specific q point
% Shinibali, NOV 4,2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = pairing_bare_interactions(inputs)

fdid = '/data/shinibali/kagome/oneorb/0p08/';

ncont = inputs(1); U = inputs(2); J = inputs(3); Vnn = inputs(4) ;

fres = './bareresults/'; fin = './input_files/'; addpath('./functions/')
fw = 'NsiteNorbNkxNqxNr.bin'; fid = fopen([fin,fw],'r','n');
ins = fread(fid,'double'); fclose(fid); nsite=ins(1); orbitals = ins(2); Ef= ins(6);

fout = 'a1a2_b1b2.mat'; load([fin,fout],'b1b2','hexagon_vert','-mat');

[klist,bandind,len_dkf] = kpcontour(orbitals,Ef,ncont,b1b2,hexagon_vert) ;
fscat = 'scatvec.mat'; FScontour_postproc(klist,bandind,len_dkf,fres,fscat,b1b2);
load([fres,fscat],'rearrbdin','rearrklist','PQ','MQ','-mat');

fb='bandind.bin'; fid = fopen([fin,fb],'w+','n');
fwrite(fid,rearrbdin,'double'); fclose(fid);

evec_kkp(orbitals,rearrklist);
vec_emkR = exponential_kR(rearrklist) ;
smalldim = orbitals*nsite; Rblock_size = (2* smalldim)^2 ; totk = max(size(rearrklist));

%% bubble and ladder contributions
MQ_or_PQ = 1; vkkp = kkp_pairing(PQ,vec_emkR,U,J,Vnn,fin,fres,MQ_or_PQ) ; % k+kp
MQ_or_PQ = 2; vkmkp = kkp_pairing(MQ,vec_emkR,U,J,Vnn,fin,fres,MQ_or_PQ) ; % k-kp
% % change site-orbital-spin ordering [1234] to [1432] for (k-kp) subset, spin-flip processes
vkmkp = reshape(vkmkp,smalldim,smalldim,2,2,smalldim,smalldim,2,2,totk^2);
vkmkp = permute(vkmkp,[1,6,3,8,5,2,7,4,9]);
% % v(k+kp) subset - v(k-kp) subset
vkmkp = reshape(vkmkp,Rblock_size,Rblock_size,totk,totk);
vkkp = 0.5 .* (vkkp - vkmkp);  % 1/2 factor for symmetrization

%% only in case of LGE only from bare int matrix:
if ~exist('vkkp','var')
    vkkp = zeros(Rblock_size,Rblock_size,totk,totk);
end

%% adding bare_int(MQ,PQ) and bubble,ladder contributions from above
UJ_or_V=1; % adding only U,J terms
k_p_or_m_kp = 0; vkkp = add_bare_int_vkkp(vkkp,nsite,orbitals,U,J,Vnn,PQ,UJ_or_V,k_p_or_m_kp);
UJ_or_V=2; % adding V terms for (1) 0.5*V(k+k') and (2) 0.5*V(k-k')
k_p_or_m_kp = 1; vkkp = add_bare_int_vkkp(vkkp,nsite,orbitals,U,J,Vnn,PQ,UJ_or_V,k_p_or_m_kp);
k_p_or_m_kp = 2; vkkp = add_bare_int_vkkp(vkkp,nsite,orbitals,U,J,Vnn,MQ,UJ_or_V,k_p_or_m_kp);

fb='re_vkkp.bin'; fid = fopen([fdid,fb],'w+','n');
fwrite(fid,real(vkkp),'double'); fclose(fid);
fb='im_vkkp.bin'; fid = fopen([fdid,fb],'w+','n');
fwrite(fid,imag(vkkp),'double'); fclose(fid);

end