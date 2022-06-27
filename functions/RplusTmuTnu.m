function [RTmuTnu] = RplusTmuTnu(nsite,norb,a1a2,TABC,Rvec)

Ham_ind = basis_formation(nsite,norb) ;
nr = size(Rvec,1); nwan = nsite*norb ;
RTmuTnu = zeros(nwan,nwan,2,nr); % 2 for x and y components

Rvec = Rvec(:,1)*a1a2(1,:) + Rvec(:,2)*a1a2(2,:);

for ir = 1:nr
    for s2=1:nsite
        %         for o2 =1:norb
        for s1=1:nsite
            %                 for o1 =1:norb
            row = Ham_ind( :,s1 ) ;
            col = Ham_ind( :,s2 ) ;
            r = Rvec(ir,:) ;%+ TABC(s1,:) + TABC(s2,:) ;
            RTmuTnu(row,col,1,ir) = r(1);
            RTmuTnu(row,col,2,ir) = r(2);
            %                 end
        end
        %         end
    end
end

end