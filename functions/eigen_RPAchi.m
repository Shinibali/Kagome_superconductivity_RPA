function [lead_eig_RPA,chiRPA] = eigen_RPAchi(nsite,norb,R0ind,Int_ind,chiblock)

chiRPA = zeros( (nsite*norb)^2,(nsite*norb)^2, 2 );

for isp2 = 1:2
    row = Int_ind( :,:,:,:,isp2,isp2,R0ind );
    col = Int_ind( :,:,:,:,isp2,isp2,R0ind );
    % Charge
    chiRPA(:,:, 1) = chiRPA(:,:, 1) + chiblock(row,col) ;
    % Spin
    chiRPA(:,:, 2) = chiRPA(:,:, 2) + chiblock(row,col) ;
    for isp1 = 1:2
        if isp1 ~= isp2
            row = Int_ind( :,:,:,:,isp2,isp2,R0ind );
            col = Int_ind( :,:,:,:,isp1,isp1,R0ind );
            % Charge
            chiRPA(:,:, 1) = chiRPA(:,:, 1) + chiblock(row,col) ;
            % Spin
            chiRPA(:,:, 2) = chiRPA(:,:, 2) - chiblock(row,col) ;
        end
    end
end

lead_eig_RPA = zeros(1,2);
for i=1:2
    [~,val] = eig(chiRPA(:,:,i));
    lead_eig_RPA(i) = max(real(diag(val)));
end

end