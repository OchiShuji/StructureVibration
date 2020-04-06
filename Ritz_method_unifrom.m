%/*************************************************************************
%* File Name     : Ritz_method_uniform.m
%* Code Title    : 構造振動論: Ritz法による両端自由梁の振動解析
%* Programmer    : Shuji Ochi
%* Affiliation   : Dept. of Aeronautics & Astronautics
%* Creation Date : 2020/01/14
%* Language      : Matlab
%* Version       : 1.0.0
%**************************************************************************
%* NOTE:Reference
%*  [1]「振動論」近藤恭平
%*  [2]"Vibration of Rectangular Plates by the Ritz Method". Dana Young,et
%*      al.vol17,no4,1950,pp448-453
%*************************************************************************/

clear all; close all;  %全ての変数とFigureを削除

%/*************************************************************************

%/**********************Uniform Beam***************************************
%梁の諸元
l = 1;   %長さ(m)
E = 70*10^9;     %ヤング率(Pa)
h = 0.003;      %高さ(m)
rho = 2700;  %密度(kg/m^3)


N = 10000;
x = linspace(0,l,N);

EI = E*2*h*2*h*2*h/12;
mu = rho*2*h;


%モード関数の設定
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

phi(1,:) = ones(1,N); %1次
phi(2,:) = sqrt(3)*(ones(1,N)-(2/l)*x); %2次
for i = 3:n   %3次以上
    beta(1,i) = (sinh(lam(1,i))+sin(lam(1,i)))/(cosh(lam(1,i))-cos(lam(1,i)));
    phi(i,:) = cosh((lam(1,i)/l).*x)+cos((lam(1,i)/l).*x)-beta(1,i).*(sinh((lam(1,i)/l).*x)+sin((lam(1,i)/l).*x));
end

%モード関数の二階微分
d2phi(1,:) = zeros(1,N); %1次
d2phi(2,:) = zeros(1,N); %2次
for i = 3:n              %3次以上
    d2phi(i,:) = (lam(1,i)/l)^2.*(cosh((lam(1,i)/l).*x)-cos((lam(1,i)/l).*x)-beta(1,i).*(sinh((lam(1,i)/l).*x)-sin((lam(1,i)/l).*x)));
end

%剛性マトリックスと質量マトリックスの計算
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

%固有振動数の計算
[C,ome2] = eig(k,m);   
[ome,ind] = sort(diag(sqrt(ome2)));  %固有角振動数
freq = ome./2./pi;    %固有振動数
ome_0 = ome./(1/l^2*sqrt(EI/mu));     %固有角振動数(無次元化)
ome_t = (lam./l).*(lam./l).*sqrt(EI/mu);   %理論値
W = zeros(n,N);       %固有モード


%固有モードの図示
for i = 1:n
    for j = 1:n
        W(i,:) = W(i,:)+C(j,i).*phi(j,:);
    end
    if i == ind(1+2)       %1次モードの位相をtaperedに合わせる         
        W(i,:)= 1*W(i,:);
    end
    if i == ind(2+3)       %2次モードの位相をtaperedに合わせる 
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
