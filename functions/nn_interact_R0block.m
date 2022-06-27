%###########################################
% Function: getH
%###########################################
function [rows,cols,val_Vnn_fullq] = nn_interact_R0block(qx,qy,nsite,norb,Vnn,Int_ind,R0ind)

a1a2 = [1.0,  0.0;
    1/2,  sqrt(3)/2];

b1b2 = 2. * pi * eye(2,2)/ transpose(a1a2) ;

a123 = [1/2,   0;
    1/4,  sqrt(3)/4;
    -1/4,  sqrt(3)/4];

phiAB = @(x,y) 1 + exp( 2i *( x*a123(1,1) + y*a123(1,2) ) ) ;
phiAC = @(x,y) 1 + exp( 2i *( x*a123(2,1) + y*a123(2,2) ) ) ;
phiBC = @(x,y) 1 + exp( 2i *( x*a123(3,1) + y*a123(3,2) ) ) ;

qxy = qx * b1b2(1,:) + qy * b1b2(2,:) ; nqxy = max(size(qxy)) ;

% For density-density interaction, fermion legs 1=4 and 3=2 and matrix is
% constructed with joint indices (1&2, 3&4). Next, finding the NN interaction
% terms for R0 block:
for iq= 1:nqxy
    rows=zeros(1,1) ; cols=rows ; vals=rows ;
    xx = qxy(iq,1) ; yy = qxy(iq,2);
    for isp2 = 1:2
        for isp1 = 1:2
            for o2 = 1:norb
                for o1 = 1:norb
                    for is1= 1 : nsite-1
                        for is2 = is1+1 : nsite
                            if is1 == 1 && is2 == 2
                                value1 = Vnn*phiAB(xx,yy) ;
                            elseif is1 == 1 && is2 == 3
                                value1 = Vnn*phiAC(xx,yy) ;
                            elseif is1 == 2 && is2 == 3
                                value1 = Vnn*phiBC(xx,yy) ;
                            end
                            value2 = conj(value1) ;
                            % setting -V(1+exp(iq.2a)) from bubble diagrams
                            rw1 = Int_ind( o1,is1,o1,is1,isp1,isp1,R0ind );
                            cl1 = Int_ind( o2,is2,o2,is2,isp2,isp2,R0ind );
                            
                            rw2 = Int_ind( o2,is2,o2,is2,isp2,isp2,R0ind );
                            cl2 = Int_ind( o1,is1,o1,is1,isp1,isp1,R0ind );
                            
                            rows = vertcat(rows, rw1(:) ); % joint 1&2
                            cols = vertcat(cols, cl1(:) ); % joint 3&4
                            vals = vertcat(vals, -value1 .* ones(size(rw1(:))) );
                            rows = vertcat(rows, rw2(:) ); % joint 1&2
                            cols = vertcat(cols, cl2(:) ); % joint 3&4
                            vals = vertcat(vals, -value2 .* ones(size(rw2(:))) );
                            
                            % setting +V for intra-cell ladder diagrams
                            rw1 = Int_ind( o1,is1,o2,is2,isp1,isp2,R0ind );
                            cl1 = Int_ind( o2,is2,o1,is1,isp2,isp1,R0ind );
                            
                            rw2 = Int_ind( o2,is2,o1,is1,isp2,isp1,R0ind );
                            cl2 = Int_ind( o1,is1,o2,is2,isp1,isp2,R0ind );
                            
                            rows = vertcat(rows, rw1(:) ); % joint 1&2
                            cols = vertcat(cols, cl1(:) ); % joint 3&4
                            vals = vertcat(vals, Vnn .* ones(size(rw1(:))) );
                            rows = vertcat(rows, rw2(:) ); % joint 1&2
                            cols = vertcat(cols, cl2(:) ); % joint 3&4
                            vals = vertcat(vals, Vnn .* ones(size(rw2(:))) );
                        end
                    end
                end
            end
        end
    end
    rows=rows(2:end) ; cols=cols(2:end) ; vals=vals(2:end) ;
    val_Vnn_fullq(:,iq) = vals ;
end

end