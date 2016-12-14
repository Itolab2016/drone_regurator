
g= 9.80665

Ixx = 0.082
Iyy = 0.083
Izz = 0.062
Ixz = 0.000
mass = 0.712
l = 0.137

den = (Ixx*Izz - Ixz^2)
ApL = Izz/den
ApN = Ixz/den
AqM = 1/Iyy
ArL = Ixz/den
ArN = Ixx/den

K1 = 3.321*g
K2 = 3.426*g
K3 = 3.622*g
K4 = 3.454*g

KT1 = 0.06329
KT2 = 0.04349
KT3 = 0.02104
KT4 = 0.05735


B31 = -K1/mass
B32 = -K2/mass
B33 = -K3/mass
B34 = -K4/mass

B41 = -ApL*l*K1 + ApN*KT1
B42 =  ApL*l*K2 + ApN*KT2
B43 = -ApN*KT3 
B44 = -ApN*KT4

B53 =  AqM*l*K3
B54 = -AqM*l*K4

B61 = -ArL*l*K1 + ArN*KT1
B62 =  ArL*l*K2 + ArN*KT2
B63 = -ArN*KT3
B64 = -ArN*KT4

A=[ 0 0 0 0 0 0 ;
    0 0 0 0 0 0 ;
    0 0 0 0 0 0 ;
    1 0 0 0 0 0 ;
    0 1 0 0 0 0 ;
    0 0 1 0 0 0 ]
    
B=[ B41 B42 B43 B44;
    0   0   B53 B54;
    B61 B62 B63 B64;
    0 0 0 0;
    0 0 0 0;
    0 0 0 0 ]
    
Csim = eye(6,6)
C=[ 0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1]

D = zeros(6,4)

// Contoroller Desigen
Q=diag([0.3,0.3,0.3,3,3,0.3]);R=diag([1,1,1,1]);     //Usual notations x'Qx + u'Ru
Big=sysdiag(Q,R);    //Now we calculate C1 and D12
[w,wp]=fullrf(Big);C1=wp(:,1:6);D12=wp(:,7:$);   //[C1,D12]'*[C1,D12]=Big
P=syslin('c',A,B,C1,D12);    //The plant (continuous-time)
[K,X]=lqr(P)
disp(spec(A+B*K))    //check stability
disp(norm(A'*X+X*A-X*B*inv(R)*B'*X+Q,1))  //Riccati check


L=ppol(A',C',[-1 -1 -1 -1 -1 -1 ])

printf("{%10f, %10f, %10f, %10f, %10f, %10f},\n",K(1,1),K(1,2),K(1,3),K(1,4),K(1,5),K(1,6))
printf("{%10f, %10f, %10f, %10f, %10f, %10f},\n",K(2,1),K(2,2),K(2,3),K(2,4),K(2,5),K(2,6))
printf("{%10f, %10f, %10f, %10f, %10f, %10f},\n",K(3,1),K(3,2),K(3,3),K(3,4),K(3,5),K(3,6))
printf("{%10f, %10f, %10f, %10f, %10f, %10f}\n",K(4,1),K(4,2),K(4,3),K(4,4),K(4,5),K(4,6))

