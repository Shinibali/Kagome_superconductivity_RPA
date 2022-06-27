%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: k,k' pairing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ vkkp ] = kkp_pairing(MPQ,expmkR,U,J,Vnn,fin,fres,MQ_or_PQ)

fwkq = 'NsiteNorbNkxNqxNr.bin';
fid = fopen([fin,fwkq],'r','n'); ins = fread(fid,'double'); fclose(fid);
orbitals = ins(2); nsite=ins(1); nqgrid = ins(4); nR = ins(7); totq = nqgrid^2;

% finding unique k+-k' values
bzwidth = 1; totk = size(MPQ,1);
MPQ = reshape(MPQ,[],2); [uniq,~,ium] = unique(MPQ,'rows');
newq = mod(uniq, bzwidth); [finalq,~,ifm] = unique(newq,'rows');
clear newq uniq; nq = max(size(finalq));

% To generate back original MPQ, the following will do:
kkp_uniind = zeros(totk^2,1); kkp_uniind(:) = ifm(ium(:));
kkp_uniind = reshape(kkp_uniind,totk,totk) ;

% evaluating U,J,V interaction matrix
[~,~,Int_ind] = basis_formation(nsite,orbitals); dimen = numel(Int_ind);
[Intmat,Identity,R0ind] = nn_interaction(nsite,orbitals,U,J,Vnn);
[rows,cols,val_Vnn_fullq] = nn_interact_R0block(finalq(:,1),finalq(:,2),...
    nsite,orbitals,Vnn,Int_ind,R0ind);

% initializing vkkp(:,:,totk,totk) block
Rblock_size = (2* nsite * orbitals)^2 ;
vkkp = zeros(Rblock_size,Rblock_size,totk,totk);

% preparing for the interpolation scheme
fqp = 'qpoints.bin'; fid = fopen([fin,fqp],'r','n');
regq = fread(fid,[totq 2],'double'); fclose(fid); qpts = sqrt(totq);
X = reshape(regq(:,1),qpts,qpts); Y = reshape(regq(:,2),qpts,qpts);
P = [2,1]; X = permute(X,P) ; Y = permute(Y,P) ; clear regq;
chidimen = (orbitals*nsite)^4 * nR^2 ;
fout = 'baresus.mat'; load([fres,fout],'chio','-mat');
chio = reshape(chio,[],qpts,qpts) ; chio = permute(chio,[1,3,2]) ;
typ='cubic'; F = cell(chidimen,1) ;
for t = 1 : chidimen
    chi(:,:) = chio(t,:,:);
    F{t} = griddedInterpolant(X,Y,chi,typ);
end

intrchi = zeros(chidimen,1); smalldim = orbitals*nsite; d1 = (smalldim*2)^2 * nR ;
if MQ_or_PQ == 1 % i.e. k+kp
    exppkR = conj(expmkR);
elseif MQ_or_PQ == 2 % i.e. k-kp 
    exppkR = expmkR;
end
for iq = 1:nq
    qx = finalq(iq,1); qy = finalq(iq,2) ;
    % interpolate chio for each finalq
    for t = 1 : chidimen
        intrchi(t) = F{t}(qx,qy);
    end
    intrchi = reshape(intrchi, smalldim,smalldim,smalldim,smalldim,nR,nR ) ;
    
    % for each finalq, copy chio to full spin-basis
    intrchi = permute(intrchi,[1,2,5,3,4,6]);
    Achi = zeros(smalldim,smalldim,2,2,nR,smalldim,smalldim,2,2,nR) ;
    for isp2 = 1:2
        for isp1 = 1:2
            Achi(:,:,isp1,isp2,:,:,:,isp2,isp1,:) = intrchi(:,:,:,:,:,:);
        end
    end
    Achi = reshape(Achi,d1,d1);
    
    % for each finalq, find vpair = Uint*chirpa*Uint
    int_nn_R0 = sparse(rows,cols,val_Vnn_fullq(:,iq),dimen,dimen);
    Interaction = Intmat + int_nn_R0 ;
    vpair = Identity/( Identity - Achi * Interaction ) * Achi ;
    vpair = Interaction * vpair * Interaction ;
    
    % for k,kp pair, evaluate vkkp summed over exp(ik.R)
    [krow,kcol] = find(kkp_uniind==iq);
    
    for ikkp = 1:numel(krow)
        k = krow(ikkp); kp = kcol(ikkp) ;
        for r1 = 1:nR
            start_col = (r1-1) * Rblock_size + 1 ;
            end_col = r1 * Rblock_size ;
            for r2 = 1:nR
                start_row = (r2-1) * Rblock_size + 1 ;
                end_row = r2 * Rblock_size ;
                vkkp(:,:,k,kp) = vkkp(:,:,k,kp) + ( expmkR(r2,k) * exppkR(r1,kp) ) ...
                    .* vpair( start_row:end_row, start_col:end_col ) ;
            end
        end
    end
end

end