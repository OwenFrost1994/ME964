clear
clc
close all

volf=[0.01, 0.012, 0.02, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]; %volume fraction
uc = 55; %shear contrast of fiber to matrix

lambda2s=1:-0.01:0.7;
k2s=[0.0001 0.01:0.01:5];


for n=1:1:size(volf,2)
    %% material and workpiece variables
    c2=volf(n);     %volume fraction of fiber
    c1=1-c2;     %volume fraction of matrix
    
    zerothres=1E-15;
    if c2>0.2
       zerothres=1.5E-8; 
    end

    %% calculate
    u1=1;        %shear modules of matrix
    u2=u1*uc;    %shear modules of fiber
    L0=1;
    L01=L0*c1;   %length of matrix
    L02=L0*c2;   %length of fiber
    
    for i=1:1:size(lambda2s,2)
        for j=1:1:size(k2s,2)
            %% deformation
            lambda2=lambda2s(i);
            lambda1=1/lambda2;
            k2=k2s(j);
            F=[[lambda1,0];[0,lambda2]];
            C1=Ciqkp(u1,F);
            C2=Ciqkp(u2,F);
            
            p1=u1*F(1,1)^2;
            p2=u2*F(1,1)^2;
            
            %% matrices
            V1=zeros(4,4);
            V2=zeros(4,4);
            V1(1,2)=-1i*k2;   V1(2,3)=1;   V1(3,2)=-k2^2*(C1(2,1,1,2)+C1(2,2,1,1)-C1(2,2,2,2))/C1(2,1,2,1);
            V2(1,2)=-1i*k2;   V2(2,3)=1;   V2(3,2)=-k2^2*(C2(2,1,1,2)+C2(2,2,1,1)-C2(2,2,2,2))/C2(2,1,2,1);
            V1(3,4)=1i*k2/C1(2,1,2,1);   V1(4,1)=-k2^2*C1(1,2,1,2);   V1(4,3)=-1i*k2*(C1(1,1,1,1)-C1(1,1,2,2)-C1(1,2,2,1));
            V2(3,4)=1i*k2/C2(2,1,2,1);   V2(4,1)=-k2^2*C2(1,2,1,2);   V2(4,3)=-1i*k2*(C2(1,1,1,1)-C2(1,1,2,2)-C2(1,2,2,1));
            
            [W1,Z1]=eig(V1);
            [W2,Z2]=eig(V2);
            
            Q1=zeros(4,4);
            Q2=zeros(4,4);
            Q1(1,1)=1;Q1(2,2)=1;Q1(3,4)=-1;Q1(3,2)=1i*k2*(C1(1,1,2,2)-C1(1,1,1,1)-p1);
            Q2(1,1)=1;Q2(2,2)=1;Q2(3,4)=-1;Q2(3,2)=1i*k2*(C2(1,1,2,2)-C2(1,1,1,1)-p2);
            Q1(4,1)=1i*k2*p1;Q1(4,3)=C1(2,1,2,1);
            Q2(4,1)=1i*k2*p2;Q2(4,3)=C2(2,1,2,1);
            
            G1=Q1*W1;
            G2=Q2*W2;
            
            L1=L01*lambda2;
            L2=L02*lambda2;
            L=L0*lambda1;
            
            K=pinv(G1)*G2*diag(exp(diag(Z1*L1)))*pinv(G2)*G1*diag(exp(diag(Z2*L2)));
            
            if det(G1)~=0 && det(G2)~=0
                I1=trace(K);
                I2=0.5*(trace(K)^2-trace(K^2));
                a=sqrt(I1^2-4*I2+8);

                A1=abs(abs(I1/4+a/4+sqrt(I1^2-2*I2-4+I1*a)/(2*sqrt(2)))-1)<zerothres;
                A2=abs(abs(I1/4+a/4-sqrt(I1^2-2*I2-4+I1*a)/(2*sqrt(2)))-1)<zerothres;
                A3=abs(abs(I1/4-a/4+sqrt(I1^2-2*I2-4-I1*a)/(2*sqrt(2)))-1)<zerothres;
                A4=abs(abs(I1/4-a/4-sqrt(I1^2-2*I2-4-I1*a)/(2*sqrt(2)))-1)<zerothres;
                Amax(j,i)=max([A1,A2,A3,A4]);
            else
                Amax(j,i)=0;
            end
            
        end
    end
    
    %================= identify instability in  ==================================
    % Computes onset of
    Insta=0;
    for i=1:1:size(lambda2s,2)
        for j=1:1:size(k2s,2)
            if Insta==0
                if Amax(j,i)==1 && lambda2s(i)~=1
                    lambda2cir(n)=lambda2s(i);
                    k2cir(n)=k2s(j);
                    Insta=1;
                end
            end
        end
    end

end

figure(size(volf,2)+1)
subplot(1,2,1)
plot(volf,lambda2cir,'rs-','Markersize',12,'LineWidth',3)
hold on;
xlabel({'c^{(F)}'},'FontSize',20);
ylabel({'\lambda^{cr}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(1,2,2)
plot(volf,k2cir,'r^-','Markersize',12,'LineWidth',3)
xlabel({'c^{(F)}'},'FontSize',20);
ylabel({'k_2^{cr}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)