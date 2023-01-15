clear
clc
close all

%% 1. TI fiber composites Relevant components of the effective stress 
%%tensor for various loading conditions
mu_M = 1;
shearcont = 55;
c_F = 0.1;
c_M = 1-c_F;
mu_F = shearcont*mu_M;
N_0=[0; 1; 0];
mu_w = mu_M*(((1+c_F)*mu_F+(1-c_F)*mu_M)/((1-c_F)*mu_F+(1+c_F)*mu_M));
mu_b = c_M*mu_M+c_F*mu_F;
%% uniaxial compression
lamb = 1.2:-0.01:0.6;

for i = 1:1:size(lamb,2)
    F = diag([1/sqrt(lamb(i)),lamb(i),1/sqrt(lamb(i))]);
    I4 = N_0'*F*F'*N_0;
    p = mu_w/lamb(i);
    sig = mu_w*F*F' + (mu_b-mu_w)*(1-I4^(1-3/2))*(F*N_0)*(F*N_0)' - p*eye(3);
    sig_11(i) = sig(1,1);
    sig_22(i) = sig(2,2);
    sig_33(i) = sig(3,3);
end

figure(1)
subplot(1,3,1)
plot(lamb,sig_11,'b-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\lambda_2'},'FontSize',20);
ylabel({'\sigma_{11}'},'FontSize',20);
ylim([-1 1])
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(1,3,2)
plot(lamb,sig_22,'r-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\lambda_2'},'FontSize',20);
ylabel({'\sigma_{22}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(1,3,3)
plot(lamb,sig_33,'k-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\lambda_2'},'FontSize',20);
ylabel({'\sigma_{33}'},'FontSize',20);
ylim([-1 1])
set(gca, 'FontName','Times New Roman','FontSize', 20)

clearvars sig_11 sig_22 sig_33

%% simple shear
gamma = 0:0.001:0.1;

for i = 1:1:size(gamma,2)
    F = [1,gamma(i),0;
        0,1,0;
        0,0,1];
    I4 = N_0'*F*F'*N_0;
    sig = mu_w*F*F' + (mu_b-mu_w)*(1-I4^(1-3/2))*(F*N_0)*(F*N_0)';
    sig_11(i) = sig(1,1)- sig(3,3);
    sig_22(i) = sig(2,2)- sig(3,3);
    sig_12(i) = sig(1,2);
end

figure(2)
subplot(1,3,1)
plot(gamma,sig_11,'b-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\gamma'},'FontSize',20);
ylabel({'\sigma_{11}-\sigma_{33}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(1,3,2)
plot(gamma,sig_22,'r-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\gamma'},'FontSize',20);
ylabel({'\sigma_{22}-\sigma_{33}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(1,3,3)
plot(gamma,sig_12,'k-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\gamma'},'FontSize',20);
ylabel({'\sigma_{12}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)

clearvars sig_11 sig_22 sig_12

%% pure shear
gamma = 0:0.001:0.1;

for i = 1:1:size(gamma,2)
    F = [1+gamma(i)^2,gamma(i),0;
        gamma(i),1,0;
        0,0,1];
    I4 = N_0'*F*F'*N_0;
    sig = mu_w*F*F' + (mu_b-mu_w)*(1-I4^(1-3/2))*(F*N_0)*(F*N_0)';
    sig_11(i) = sig(1,1) - sig(3,3);
    sig_22(i) = sig(2,2) - sig(3,3);
    sig_12(i) = sig(1,2);
end
figure(3)
subplot(1,3,1)
plot(gamma,sig_11,'b-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\gamma'},'FontSize',20);
ylabel({'\sigma_{11}-\sigma_{33}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(1,3,2)
plot(gamma,sig_22,'r-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\gamma'},'FontSize',20);
ylabel({'\sigma_{22}-\sigma_{33}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(1,3,3)
plot(gamma,sig_12,'k-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\gamma'},'FontSize',20);
ylabel({'\sigma_{12}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)

clearvars sig_11 sig_22 sig_12

%% 2. For neo-Hookean laminates, determine the stresses and deformation
%%developing in the phases (fiber and matrix) as functions of the macroscopically
%%applied stretch and shear
%% uniaxial compression
lamb = 1.2:-0.01:0.6;

for i = 1:1:size(lamb,2)
    F = diag([1/lamb(i),lamb(i)]);
    p_M = mu_M*(1/lamb(i))^2;
    p_F = mu_F*(1/lamb(i))^2;
    
    sigM = det(F)*(mu_M*F*F'-p_M*eye(2));
    sigM_11(i) = sigM(1,1);
    sigM_22(i) = sigM(2,2);
    
    sigF = det(F)*(mu_F*F*F'-p_F*eye(2));
    sigF_11(i) = sigF(1,1);
    sigF_22(i) = sigF(2,2);
end

figure(4)
subplot(2,2,1)
plot(lamb,sigM_11,'b-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\lambda_2'},'FontSize',20);
ylabel({'\sigma^{(M)}_{11}'},'FontSize',20);
ylim([-1 1])
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(2,2,2)
plot(lamb,sigM_22,'r-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\lambda_2'},'FontSize',20);
ylabel({'\sigma^{(M)}_{22}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(2,2,3)
plot(lamb,sigF_11,'b-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\lambda_2'},'FontSize',20);
ylabel({'\sigma^{(F)}_{11}'},'FontSize',20);
ylim([-1 1])
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(2,2,4)
plot(lamb,sigF_22,'r-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\lambda_2'},'FontSize',20);
ylabel({'\sigma^{(F)}_{22}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)

clearvars sigM_11 sigM_22 sigF_11 sigF_22

%% simple shear
gamma = 0:0.001:0.1;

for i = 1:1:size(gamma,2)
    F = [1,gamma(i);
        0,1;];
    
    p_M = mu_M;
    p_F = mu_F;
    
    sigM = det(F)*(mu_M*F*F'-p_M*eye(2));
    sigM_11(i) = sigM(1,1);
    sigM_22(i) = sigM(2,2);
    sigM_12(i) = sigM(1,2);
    
    sigF = det(F)*(mu_F*F*F'-p_F*eye(2));
    sigF_11(i) = sigF(1,1);
    sigF_22(i) = sigF(2,2);
    sigF_12(i) = sigF(1,2);
end

figure(5)
subplot(2,3,1)
plot(gamma,sigM_11,'b-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\gamma'},'FontSize',20);
ylabel({'\sigma^{(M)}_{11}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(2,3,2)
plot(gamma,sigM_22,'r-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\gamma'},'FontSize',20);
ylabel({'\sigma^{(M)}_{22}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(2,3,3)
plot(gamma,sigM_12,'k-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\gamma'},'FontSize',20);
ylabel({'\sigma^{(M)}_{12}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(2,3,4)
plot(gamma,sigF_11,'b-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\gamma'},'FontSize',20);
ylabel({'\sigma^{(F)}_{11}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(2,3,5)
plot(gamma,sigF_22,'r-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\gamma'},'FontSize',20);
ylabel({'\sigma^{(F)}_{22}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)
subplot(2,3,6)
plot(gamma,sigF_12,'k-','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'\gamma'},'FontSize',20);
ylabel({'\sigma^{(F)}_{12}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)

clearvars sigM_11 sigM_22 sigM_12 sigF_11 sigF_22 sigF_12

%% 3. Plot the critical stretch ratio vs. fiber volume fraction for
%the assigned matrix-to-fiber initial shear moduli contrast.
%Critical stretch ratio versus fiber volume fraction
mu_M = 1;
shearcont1 = 20;
shearcont2 = 55;
shearcont3 = 80;
c_Fl = 0.05:0.01:0.95;
for i = 1:1:size(c_Fl,2)
    c_F = c_Fl(i);
    c_M = 1-c_F;
    
    mu_F = shearcont1*mu_M;
    mu_w = mu_M*(((1+c_F)*mu_F+(1-c_F)*mu_M)/((1-c_F)*mu_F+(1+c_F)*mu_M));
    mu_b = c_M*mu_M+c_F*mu_F;
    lamb_cr1(i) = (1-mu_w/mu_b)^(1/3);
    
    mu_F = shearcont2*mu_M;
    mu_w = mu_M*(((1+c_F)*mu_F+(1-c_F)*mu_M)/((1-c_F)*mu_F+(1+c_F)*mu_M));
    mu_b = c_M*mu_M+c_F*mu_F;
    lamb_cr2(i)=(1-mu_w/mu_b)^(1/3);
    
    mu_F = shearcont3*mu_M;
    mu_w = mu_M*(((1+c_F)*mu_F+(1-c_F)*mu_M)/((1-c_F)*mu_F+(1+c_F)*mu_M));
    mu_b = c_M*mu_M+c_F*mu_F;
    lamb_cr3(i)=(1-mu_w/mu_b)^(1/3);
end

figure(6)
plot(c_Fl,lamb_cr1,'b-','Markersize',6,'LineWidth',2)
hold on
plot(c_Fl,lamb_cr2,'r--','Markersize',6,'LineWidth',2)
hold on
plot(c_Fl,lamb_cr3,'k-.','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'c^{(F)}'},'FontSize',20);
ylabel({'\lambda_{cr}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)
legend(strcat('\mu^{(F)}/\mu^{(M)}=',num2str(shearcont1)),strcat('\mu^{(F)}/\mu^{(M)}=',num2str(shearcont2)),strcat('\mu^{(F)}/\mu^{(M)}=',num2str(shearcont3)))

shearcontl = 5:1:100;
c_F1 = 0.05;
c_F2 = 0.1;
c_F3 = 0.2;
for i = 1:1:size(shearcontl,2)
    mu_F = shearcontl(i)*mu_M;
    
    c_F = c_F1;
    c_M = 1-c_F;
    mu_w = mu_M*(((1+c_F)*mu_F+(1-c_F)*mu_M)/((1-c_F)*mu_F+(1+c_F)*mu_M));
    mu_b = c_M*mu_M+c_F*mu_F;
    lamb_cr4(i)=(1-mu_w/mu_b)^(1/3);
    
    c_F = c_F2;
    c_M = 1-c_F;
    mu_w = mu_M*(((1+c_F)*mu_F+(1-c_F)*mu_M)/((1-c_F)*mu_F+(1+c_F)*mu_M));
    mu_b = c_M*mu_M+c_F*mu_F;
    lamb_cr5(i)=(1-mu_w/mu_b)^(1/3);
    
    c_F = c_F3;
    c_M = 1-c_F;
    mu_w = mu_M*(((1+c_F)*mu_F+(1-c_F)*mu_M)/((1-c_F)*mu_F+(1+c_F)*mu_M));
    mu_b = c_M*mu_M+c_F*mu_F;
    lamb_cr6(i)=(1-mu_w/mu_b)^(1/3);
end

figure(7)
plot(shearcontl,lamb_cr4,'b-','Markersize',6,'LineWidth',2)
hold on
plot(shearcontl,lamb_cr5,'r--','Markersize',6,'LineWidth',2)
hold on
plot(shearcontl,lamb_cr6,'k-.','Markersize',6,'LineWidth',2)
hold on
grid on
xlabel({'c^{(F)}'},'FontSize',20);
ylabel({'\lambda_{cr}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)
legend(strcat('c^{(F)}=',num2str(c_F1)),strcat('c^{(F)}=',num2str(c_F2)),strcat('c^{(F)}=',num2str(c_F3)))

%% 4. Plot the critical stretch and wavenumber vs. fiber volume fraction.
%Compare with the explicit results obtained through the loss of ellipticity
%analysis for macroscopic instabilities (from Sec. 3).

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

c_Fl = 0.05:0.01:0.95;
figure(8)%size(volf,2)+1
subplot(1,2,1)
plot(volf,lambda2cir,'rs-','Markersize',12,'LineWidth',3)
hold on;
plot(c_Fl,lamb_cr2,'r--','Markersize',6,'LineWidth',2)
hold on
xlabel({'c^{(F)}'},'FontSize',20);
ylabel({'\lambda^{cr}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)
legend('Micro-instability method','Loss of ellipticity')
subplot(1,2,2)
plot(volf,k2cir,'r^-','Markersize',12,'LineWidth',3)
xlabel({'c^{(F)}'},'FontSize',20);
ylabel({'k_2^{cr}'},'FontSize',20);
set(gca, 'FontName','Times New Roman','FontSize', 20)