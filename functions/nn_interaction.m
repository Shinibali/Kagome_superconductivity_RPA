function [Intmat,Identity,R0ind,rows,cols,vals] = nn_interaction(nsite,norb,U,J,Vnn)

Up = U; Jp = J;
[~,~,Int_ind] = basis_formation(nsite,norb); dimen = numel(Int_ind);
Identity = sparse(eye( dimen ));

fin = './input_files/'; fout = 'modified_hr.mat';
load([fin,fout],'Rvec','sisj','oioj','-mat');

if norb == 1
    indices = (oioj(:,1)+oioj(:,2)) == 2 ; % selecting the dxz orbital with orbitals number 1
    Rvec = Rvec(indices,:); sisj = sisj(indices,:);
end

indices = sisj(:,1) ~= sisj(:,2) ;
Rvec = Rvec(indices,:); sisj = sisj(indices,:);
[uniR,~,iRv] = unique(Rvec,'rows'); nR = size(uniR,1);
R0ind = find( sqrt( uniR(:,1).^2 + uniR(:,2).^2 ) == 0 ); %checking for R=(0,0) point

% For density-density interaction, fermion legs 1=4 and 3=2 and matrix is
% constructed with joint indices (1&2, 3&4). The int matrix has block
% diagonal form only. The Rvec list has hopping elements from (R=0,si,oi) to
% (R=R,sj,oj). Next, finding the NN and on-site interaction terms:
rows=zeros(1,1) ; cols=rows ; vals=rows ;
for ir = 1:nR
    subset_sisj = sisj(iRv == ir,:);
    
    for isp2 = 1:2
        for isp1 = 1:2
            for o2 = 1:norb
                for o1 = 1:norb
                    if ir ~= R0ind % R~= 0
                        for isij = 1:size(subset_sisj,1)
                            s1 = subset_sisj(isij,1); s2 = subset_sisj(isij,2) ;
                            rw = Int_ind( o1,s1,o2,s2,isp2,isp1,ir );
                            cl = Int_ind( o2,s2,o1,s1,isp1,isp2,ir );
                            rows = vertcat(rows, rw(:) ); % joint 1&2
                            cols = vertcat(cols, cl(:) ); % joint 3&4
                            vals = vertcat(vals, Vnn .* ones(size(rw(:))) );
                            
%                             rw = Int_ind( o1,s1,o1,s1,isp1,isp1,ir );
%                             cl = Int_ind( o2,s2,o2,s2,isp2,isp2,ir );
%                             rows = vertcat(rows, rw(:) ); % joint 1&2
%                             cols = vertcat(cols, cl(:) ); % joint 3&4
%                             vals = vertcat(vals, -Vnn .* ones(size(rw(:))) );
                        end
                    elseif ir == R0ind % same site, R0 block
                        
                        for s1 = 1:nsite
                            if o1 == o2 && isp2 ~= isp1
                                value = U;
                                
                                rw = Int_ind( o1,s1,o1,s1,isp2,isp1,R0ind );
                                cl = Int_ind( o1,s1,o1,s1,isp1,isp2,R0ind );
                                rows = vertcat(rows, rw(:) ); % joint 1&2
                                cols = vertcat(cols, cl(:) ); % joint 3&4
                                vals = vertcat(vals, value .* ones(size(rw(:))) );
                                
                                rw = Int_ind( o1,s1,o1,s1,isp2,isp2,R0ind );
                                cl = Int_ind( o1,s1,o1,s1,isp1,isp1,R0ind );
                                rows = vertcat(rows, rw(:) ); % joint 1&2
                                cols = vertcat(cols, cl(:) ); % joint 3&4
                                vals = vertcat(vals, -value .* ones(size(rw(:))) );
                                
                            elseif o1 ~= o2 && isp2 == isp1
                                value = Up-J ;
                                
                                rw = Int_ind( o1,s1,o2,s1,isp2,isp1,R0ind );
                                cl = Int_ind( o2,s1,o1,s1,isp1,isp2,R0ind );
                                rows = vertcat(rows, rw(:) ); % joint 1&2
                                cols = vertcat(cols, cl(:) ); % joint 3&4
                                vals = vertcat(vals, value .* ones(size(rw(:))) );
                                
                                rw = Int_ind( o1,s1,o1,s1,isp2,isp2,R0ind );
                                cl = Int_ind( o2,s1,o2,s1,isp1,isp1,R0ind );
                                rows = vertcat(rows, rw(:) ); % joint 1&2
                                cols = vertcat(cols, cl(:) ); % joint 3&4
                                vals = vertcat(vals, -value .* ones(size(rw(:))) );
                                
                            elseif o1 ~= o2 && isp2 ~= isp1
                                value = Up ;
                                
                                rw = Int_ind( o1,s1,o2,s1,isp2,isp1,R0ind );
                                cl = Int_ind( o2,s1,o1,s1,isp1,isp2,R0ind );
                                rows = vertcat(rows, rw(:) ); % joint 1&2
                                cols = vertcat(cols, cl(:) ); % joint 3&4
                                vals = vertcat(vals, value .* ones(size(rw(:))) );
                                
                                rw = Int_ind( o1,s1,o1,s1,isp2,isp2,R0ind );
                                cl = Int_ind( o2,s1,o2,s1,isp1,isp1,R0ind );
                                rows = vertcat(rows, rw(:) ); % joint 1&2
                                cols = vertcat(cols, cl(:) ); % joint 3&4
                                vals = vertcat(vals, -value .* ones(size(rw(:))) );
                                
                                value = J ;
                                
                                rw = Int_ind( o1,s1,o1,s1,isp2,isp1,R0ind );
                                cl = Int_ind( o2,s1,o2,s1,isp1,isp2,R0ind );
                                rows = vertcat(rows, rw(:) ); % joint 1&2
                                cols = vertcat(cols, cl(:) ); % joint 3&4
                                vals = vertcat(vals, value .* ones(size(rw(:))) );
                                
                                rw = Int_ind( o1,s1,o2,s1,isp2,isp2,R0ind );
                                cl = Int_ind( o2,s1,o1,s1,isp1,isp1,R0ind );
                                rows = vertcat(rows, rw(:) ); % joint 1&2
                                cols = vertcat(cols, cl(:) ); % joint 3&4
                                vals = vertcat(vals, -value .* ones(size(rw(:))) );
                                
                                value = Jp ;
                                
                                rw = Int_ind( o1,s1,o2,s1,isp2,isp1,R0ind );
                                cl = Int_ind( o1,s1,o2,s1,isp1,isp2,R0ind );
                                rows = vertcat(rows, rw(:) ); % joint 1&2
                                cols = vertcat(cols, cl(:) ); % joint 3&4
                                vals = vertcat(vals, value .* ones(size(rw(:))) );
                                
                                rw = Int_ind( o1,s1,o2,s1,isp2,isp2,R0ind );
                                cl = Int_ind( o1,s1,o2,s1,isp1,isp1,R0ind );
                                rows = vertcat(rows, rw(:) ); % joint 1&2
                                cols = vertcat(cols, cl(:) ); % joint 3&4
                                vals = vertcat(vals, -value .* ones(size(rw(:))) );
                            end
                        end
                    end
                end
            end
        end
    end
end
rows=rows(2:end) ; cols=cols(2:end) ; vals=vals(2:end) ;
Intmat = sparse(rows,cols,vals,dimen,dimen);

end