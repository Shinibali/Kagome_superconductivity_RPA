function [Ham_ind,chi_ind,Int_ind,short_Int_ind,bare_int_ind] = basis_formation(nsite,norb)

Ham_ind = zeros(norb,nsite) ; count = 0;
for ns = 1:nsite % site
    for no = 1:norb % orbital
        count = count+1;
        Ham_ind( no,ns ) = count ;
    end
end

chi_ind = zeros(norb,nsite,norb,nsite) ; count2 = 0;
for ns1 = 1:nsite
    for no1 = 1:norb
        for ns2 = 1:nsite
            for no2 = 1:norb
                % ind = no2 + (no1-1) * norb + (ns-1) * norb * norb ;
                count2 = count2+1;
                chi_ind( no2,ns2,no1,ns1 ) = count2;
            end
        end
    end
end

fin = './input_files/'; fout = 'modified_hr.mat'; load([fin,fout],'Rvec','-mat');
Rvec = unique(Rvec,'rows'); nR = size(Rvec,1);
Int_ind = zeros(norb,nsite,norb,nsite,2,2,nR) ; count2 = 0;
for iR = 1:nR
    for isp1 = 1:2
        for isp2 = 1:2
            for ns2 = 1:nsite
                for no2 = 1:norb
                    for ns1 = 1:nsite
                        for no1 = 1:norb
                            count2 = count2+1;
                            Int_ind( no1,ns1,no2,ns2,isp2,isp1,iR ) = count2;
                        end
                    end
                end
            end
        end
    end
end

short_Int_ind = zeros(norb*nsite,norb*nsite,2,2,nR) ; count2 = 0;
for iR = 1:nR
    for isp1 = 1:2
        for isp2 = 1:2
            for nso2 = 1:nsite*norb
                for nso1 = 1:nsite*norb
                    count2 = count2+1;
                    short_Int_ind( nso1,nso2,isp2,isp1,iR ) = count2;
                end
            end
        end
    end
end

bare_int_ind = zeros(norb,nsite,norb,nsite,2,2) ; count2 = 0;
for isp1 = 1:2
    for isp2 = 1:2
        for ns2 = 1:nsite
            for no2 = 1:norb
                for ns1 = 1:nsite
                    for no1 = 1:norb
                        count2 = count2+1;
                        bare_int_ind( no1,ns1,no2,ns2,isp2,isp1 ) = count2;
                    end
                end
            end
        end
    end
end


end