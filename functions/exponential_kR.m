function [vec_emkR] = exponential_kR(klist)

fin = './input_files/'; fout = 'modified_hr.mat'; load([fin,fout],'Rvec','-mat');
[uniR] = unique(Rvec,'rows'); nR = size(uniR,1); totk = size(klist,1);

fout = 'a1a2_b1b2.mat'; load([fin,fout],'a1a2','b1b2','-mat');

klist = klist(:,1) * b1b2(1,:) + klist(:,2) * b1b2(2,:) ;
uniR  =  uniR(:,1) * a1a2(1,:) + uniR(:,2)  * a1a2(2,:) ;

vec_emkR = zeros(nR,totk);
for ik = 1:totk
    kx = klist(ik,1); ky = klist(ik,2);
    vec_emkR(1:nR,ik) = exp( - 1i .* ( kx .* uniR(1:nR,1) + ky .* uniR(1:nR,2) ) ) ;
end

end