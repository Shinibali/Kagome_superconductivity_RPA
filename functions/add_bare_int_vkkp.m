function vkkp = add_bare_int_vkkp(vkkp,nsite,orbitals,U,J,Vnn,MPQ,UJ_or_V,k_p_or_m_kp)

% finding unique k-+k' values
bzwidth = 1; totk = size(MPQ,1);
MPQ = reshape(MPQ,[],2); [uniq,~,ium] = unique(MPQ,'rows');
newq = mod(uniq, bzwidth); [finalq,~,ifm] = unique(newq,'rows');
clear newq uniq; nq = max(size(finalq));

% To generate back original MPQ, the following will do:
kkp_uniind = zeros(totk^2,1); kkp_uniind(:) = ifm(ium(:));
kkp_uniind = reshape(kkp_uniind,totk,totk) ;

% evaluating U,J,V bare interaction matrix
[~,~,~,~,bare_int_ind] = basis_formation(nsite,orbitals); dimen = numel(bare_int_ind);

if UJ_or_V == 1
    [bareIntmat] = full_bare_int(nsite,orbitals,U,J) ;
    for kp = 1:totk
        for k = 1:totk
            vkkp(:,:,k,kp) = vkkp(:,:,k,kp) + bareIntmat;
        end
    end
elseif UJ_or_V == 2
    [rows,cols,bare_Vnn_fullq] = full_bare_nn_int(finalq(:,1),finalq(:,2),...
        nsite,orbitals,Vnn,bare_int_ind,k_p_or_m_kp);
    for iq = 1:nq
        int_nn_R0 = sparse(rows,cols,bare_Vnn_fullq(:,iq),dimen,dimen);
        % for k,kp pair, evaluate vkkp summed over exp(ik.R)
        [krow,kcol] = find(kkp_uniind==iq);
        for ikkp = 1:numel(krow)
            k = krow(ikkp); kp = kcol(ikkp) ;
            % adding V terms as 0.5*( V(k+k') + V(k-k') )
            vkkp(:,:,k,kp) = vkkp(:,:,k,kp) + int_nn_R0*0.5;
        end
    end
else
    
end