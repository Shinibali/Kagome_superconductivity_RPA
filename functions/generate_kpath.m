function [kpath,str_xt,vec_list,index] = generate_kpath(npkseg)

vec_list = [0   0;
            0.5 0.5;
            2/3 1/3;
            0   0];
str_xt={'\Gamma','M','K','\Gamma'};

s1 = size(vec_list,1);  index = zeros(s1,1);
kpath = zeros( npkseg * (s1-1) - (s1-2), 2 );  
% s1-2 for excluding the middle points from double counting

for i = 2:s1
    initvec = vec_list(i-1,:);
    segvec = ( vec_list(i,:) - initvec )./npkseg ;
    index(i-1) = 1 + (i-2)*npkseg; 
    for nk = 1: npkseg        
        kpath( nk + (i-2)*npkseg ,1:2) = initvec + segvec .* (nk-1) ; 
    end
end
index(s1) =  npkseg + (i-2)*npkseg; 
kpath = [kpath; vec_list(end,:)];

end