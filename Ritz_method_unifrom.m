%/*************************************************************************
%* File Name     : Ritz_method_uniform.m
%* Code Title    : �\���U���_: Ritz�@�ɂ�闼�[���R���̐U�����
%* Programmer    : Shuji Ochi
%* Affiliation   : Dept. of Aeronautics & Astronautics
%* Creation Date : 2020/01/14
%* Language      : Matlab
%* Version       : 1.0.0
%**************************************************************************
%* NOTE:Reference
%*  [1]�u�U���_�v�ߓ�����
%*  [2]"Vibration of Rectangular Plates by the Ritz Method". Dana Young,et
%*      al.vol17,no4,1950,pp448-453
%*************************************************************************/

clear all; close all;  %�S�Ă̕ϐ���Figure���폜

%/*************************************************************************

%/**********************Uniform Beam***************************************
%���̏���
l = 1;   %����(m)
E = 70*10^9;     %�����O��(Pa)
h = 0.003;      %����(m)
rho = 2700;  %���x(kg/m^3)


N = 10000;
x = linspace(0,l,N);

EI = E*2*h*2*h*2*h/12;
mu = rho*2*h;


%���[�h�֐��̐ݒ�
n = 5;
phi = zeros(n,N);
d2phi = zeros(n,N);

lam = zeros(1,n);
lam(1,1) = 0.0;
lam(1,2) = 0.0;
lam(1,3) = 4.7300408;
lam(1,4) = 7.8532046;
lam(1,5) = 10.9956078;
%lam(1,6) = 14.1371655;
%lam(1,7) = 17.2787596;
%lam(1,8) = 20.4203522;
%lam(1,9) = 23.5619449;
%lam(1,10) = 26.7035375;

beta = zeros(1,n);

phi(1,:) = ones(1,N); %1��
phi(2,:) = sqrt(3)*(ones(1,N)-(2/l)*x); %2��
for i = 3:n   %3���ȏ�
    beta(1,i) = (sinh(lam(1,i))+sin(lam(1,i)))/(cosh(lam(1,i))-cos(lam(1,i)));
    phi(i,:) = cosh((lam(1,i)/l).*x)+cos((lam(1,i)/l).*x)-beta(1,i).*(sinh((lam(1,i)/l).*x)+sin((lam(1,i)/l).*x));
end

%���[�h�֐��̓�K����
d2phi(1,:) = zeros(1,N); %1��
d2phi(2,:) = zeros(1,N); %2��
for i = 3:n              %3���ȏ�
    d2phi(i,:) = (lam(1,i)/l)^2.*(cosh((lam(1,i)/l).*x)-cos((lam(1,i)/l).*x)-beta(1,i).*(sinh((lam(1,i)/l).*x)-sin((lam(1,i)/l).*x)));
end

%�����}�g���b�N�X�Ǝ��ʃ}�g���b�N�X�̌v�Z
k = zeros(n,n);
m = zeros(n,n);

for i = 1:n
    for j = 1:n
        f = EI*d2phi(i,:).*d2phi(j,:);
        k(i,j) = trapz(f);
        g = mu*phi(i,:).*phi(j,:);
        m(i,j) = trapz(g);
    end
end

%�ŗL�U�����̌v�Z
[C,ome2] = eig(k,m);   
[ome,ind] = sort(diag(sqrt(ome2)));  %�ŗL�p�U����
freq = ome./2./pi;    %�ŗL�U����
ome_0 = ome./(1/l^2*sqrt(EI/mu));     %�ŗL�p�U����(��������)
ome_t = (lam./l).*(lam./l).*sqrt(EI/mu);   %���_�l
W = zeros(n,N);       %�ŗL���[�h


%�ŗL���[�h�̐}��
for i = 1:n
    for j = 1:n
        W(i,:) = W(i,:)+C(j,i).*phi(j,:);
    end
    if i == ind(1+2)       %1�����[�h�̈ʑ���tapered�ɍ��킹��         
        W(i,:)= 1*W(i,:);
    end
    if i == ind(2+3)       %2�����[�h�̈ʑ���tapered�ɍ��킹�� 
        W(i,:)= -1*W(i,:);
    end
    figure(i);
    plot(x,W(i,:));
    title('mode'+string(find(ind==i)-2));
    xlabel("x[m]")
end

figure(n+1);
hold on;
for i = 1:3
    plot(x,W(ind(i+2),:));
end
legend("mode1","mode2","mode3");
title("mode1-3")
xlabel("x[m]")
