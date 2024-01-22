function [M_delta,M_R,M_D,trad,psi,delt,theta,A] = Luchipman(M)
%% MD
D = (M(1, 2:4))';
D_mag = norm(D);
D_unit = D ./ D_mag;
I = eye(3);
Z=[0 0 0];
if(D_mag>1e-15) 
    if(D_mag>1)
        D_mag = 1;
    end
    m_D = (sqrt(1 - D_mag^2) .* I) + ((1 - sqrt(1 - D_mag^2)) .* (D_unit*D_unit'));
    M_D =  [1, D'; D, m_D];
else
  M_D = [1, Z; Z', I];
end
%% MR and Mdelta
Mp=  M * pinv(M_D) ;
Mp = Mp/Mp(1,1);
PP=[Mp(2,1) Mp(3,1) Mp(4,1)]';
m0=[Mp(2,2) Mp(2,3) Mp(2,4);Mp(3,2) Mp(3,3) Mp(3,4);Mp(4,2) Mp(4,3) Mp(4,4)];
mprimedeterminant = abs(det(m0));
E=m0*m0'; 
[V,D]=eig(E);
L1=D(1,1);
L2=D(2,2);
L3=D(3,3) ;
i=mprimedeterminant/abs(mprimedeterminant);
m1=i*pinv(E+(sqrt(L1*L2)+sqrt(L2*L3)+sqrt(L1*L3)).* I)*((sqrt(L1)+sqrt(L2)+sqrt(L3)).*E+sqrt(L1*L2*L3).*I);
mR=pinv(m1) * m0; 
M_delta=[1 Z;PP m1];
M_R = [1, Z; Z', mR];
%% Parameter
M_delta11=M_delta(2,2);
M_delta22=M_delta(3,3);
trad=1-(trace( abs(M_delta))-1)/(3*1);% depolarization
psi=1/2*atan2((M_R(2,3)-M_R(3,2)),(M_R(2,2)+M_R(3,3)));% optical rotation
delt =acos((sqrt((M_R(2,2)+M_R(3,3))^2+(M_R(2,3)-M_R(3,2))^2))-1);% linear retardance
theta=90*atan2(-M_R(2,4),M_R(3,4))/pi;% orientation angle
A=(M_delta11-M_delta22)./(M_delta11+M_delta22);% linear depolarization anisotropy
end