function [chixyz] = physical_chi(nsite,norb,R0ind,Int_ind,chiblock)

pauli_0 = [1    0;
    0    1];
pauli_x = [0    1;
    1    0];
pauli_y = [0    -1i;
    1i   0];
pauli_z = [1    0;
    0    -1];

chixyz = zeros(4,1);
for isp4 = 1:2
    for isp3 = 1:2
        for isp2 = 1:2
            for isp1 = 1:2
                for s2 = 1:nsite
                    for o2 = 1:norb
                        for s1 = 1:nsite
                            for o1 = 1:norb
                                row = Int_ind( o1,s1,o1,s1,isp1,isp2,R0ind );
                                col = Int_ind( o2,s2,o2,s2,isp3,isp4,R0ind );
                                % CDO
                                chixyz(1) = chixyz(1) + chiblock(row,col) * pauli_0(isp1,isp2) * pauli_0(isp3,isp4);
                                % Spin xx
                                chixyz(3) = chixyz(3) + chiblock(row,col) * pauli_x(isp1,isp2) * pauli_x(isp3,isp4);
                                % Spin zz
                                chixyz(4) = chixyz(4) + chiblock(row,col) * pauli_z(isp1,isp2) * pauli_z(isp3,isp4);
                                if s1 ~=s2
                                    row = Int_ind( o1,s1,o2,s2,isp1,isp2,R0ind );
                                    col = Int_ind( o2,s2,o1,s1,isp2,isp1,R0ind );
                                    % CBO
                                    chixyz(2) = chixyz(2) + chiblock(row,col)./16; % factor 16 for overcounting orbitals and spins
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
chixyz = chixyz./4;
end