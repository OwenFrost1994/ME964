%================= C_{iqkp} ==================================
% Computes material elastic tensor C_{ijkl} 
function C = Ciqkp(u,F)
I=[[1,0];[0,1]];  
B=F*F';
C=zeros(2,2,2,2);

for i=1:1:2 
    for q=1:1:2 
        for k=1:1:2 
            for p=1:1:2 
                C(i,q,k,p)=(1/det(F))*(u*B(p,q)*I(i,k));
            end
        end
    end
end
end