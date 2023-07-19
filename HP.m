clc
clear all
load Cp.txt
temperaturas = [133:20:233]+273.15;
Piniciais = ones(1,length(temperaturas))*1e5;
Pfinais = ones(1,length(temperaturas))*200*1e5;
hP = 1e5; 
T0 = 298.15%temperatura de formação padrão
H0 = -46100;%entalpia padrão de formação
N = 1000;%numero de passos na etapa isobárica
for l = 1:length(temperaturas)
dT = (temperaturas(l)-T0)/N
vetorH = H0;
i = 1;
for k = 1:N
  T = T0+dT;
  T0 = T;
  while T >= Cp(i,1)
    i++;
  end
  cp = (Cp(i-1,2)*(Cp(i,1)-T)+Cp(i,2)*(T-Cp(i-1,1)))/5;
  H = H0+dT*cp;
  H0 = H;
end
Hinicial(l) = H;
function f = f(V,P,T)
  Tc = 132.41+273.15;
  Pc = 11363400;
  R = 8.3145;
  a = 0.42748*R^2*Tc^2.5/Pc; %L^2*atm/mol
  b = 0.08664*R*Tc/Pc; %L/mol
  %f = V^3*P+V^2*(-R*T-P*b)+V*(a)-a*b;
  f = P*(V-b)*(T^0.5*(V+b)*V)-R*T*(T^0.5*(V+b)*V)+a*(V-b);
endfunction
function df = df(V,P,T)
  a = 0.4225; %L^2*atm/mol
  b = 3.71*10^-5;  %L/mol
  R = 8.3145;
  df = (f(V+1e-9,P,T)-f(V-1e-9,P,T))/2e-9;
endfunction
function Vf = redlich(P,T)
  V = 1e0;%L
  for i = 1:50
    F = f(V,P,T);
    dF = df(V,P,T);
    V = V-F/dF;
  end
  Vf = V;
endfunction
P = Piniciais(l)
T = temperaturas(l)
N = (Pfinais(l)-Piniciais(l))/hP;
hT = 0.0001;
Href = Hinicial(l);
vetorH(1) = Href;
vetorP(1) = P;
xGauss = [-0.9602899 -0.7966664 -0.5255324 -0.1834346];
xGauss = [xGauss -flip(xGauss)]/2+0.5;
pGauss = [0.1012285 0.2223810 0.3137666 0.3626838];
pGauss = 0.5*[pGauss flip(pGauss)]/1.0001;
for i = 1:N
  printf("i = %i (l = %i)\n", i, l)
  quadratura = 0;  
  for k = 1:length(xGauss)%oito pontos
    hx = xGauss(k)*(hP);
    x = P+hx;
    Vteste = redlich(x,T);
    dVteste = (redlich(x,T+hT)-redlich(x,T-hT))/(2*hT);
    Y(k) = (-T*dVteste+Vteste);
    quadratura = quadratura+pGauss(k)*Y(k);
  endfor
  H = Href+hP*(quadratura);
  Href = H;
  P = P+hP;
  vetorH(i+1) = Href;
  vetorP(i+1) = P;
end
P
Pbar = P/100000
H
figure 1
switch l
  case 1
    plot(vetorH,vetorP/100000,'r')
    hold on   
  case 2
    plot(vetorH,vetorP/100000,'y')
    hold on
  case 3
    plot(vetorH,vetorP/100000,'g')
    hold on
  case 4
    plot(vetorH,vetorP/100000,'c')
    hold on
  case 5
    plot(vetorH,vetorP/100000,'b')
    hold on
  case 6
    plot(vetorH,vetorP/100000,'m')
    hold on
endswitch
endfor
xlabel("Entalpia (J/mol)")
ylabel("P (Pa)")
title("Diagrama PH da amônia")
legend('133ºC','153ºC','173ºC','193ºC','213ºC','233ºC')
figure 1
hold off