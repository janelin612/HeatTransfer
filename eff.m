%Known
w=1;    %in 2-D case, w=1
t=0.04; %thickness of fin
L=0.10; %length of fin

%Definitions
Lc=L+0.5*t;
thetaB=39.8126-TempInfinite;
Ac=w*t;
P=2*w+2*t;
M=(h*P*k*Ac*thetaB)^0.5;
m=(h*P/(k*Ac))^0.5;

Qf=M*(sinh(m*L)+(h*cosh(m*L))/(m*k))/(cosh(m*L)+(h*sinh(m*L))/(m*k))
Ee=(tanh(m*Lc))/(m*Lc)
Ef=Qf/h*Ac*thetaB
Qt=1*Ee*h*2*w*Lc+h*Ac*thetaB
