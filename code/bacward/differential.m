function [trad,psi,delt,theta,A] = differential(M)
L=real(logm(M));
G =diag([1 -1 -1 -1]);
Lm=1/2*(L-G*L'*G);
Lu=1/2*(L+G*L'*G);
theta= 90*atan2(-Lm(2,4),Lm(3,4))/pi;% orientation angle
delt =sqrt(Lm(2,4)^2+Lm(3,4)^2);% linear retardance
psi =1/2*Lm(2,3);% optical rotation
Lu2=Lu-Lu(1,1)*eye(4);
trad= 1-1/3*(exp(Lu2(2,2))+exp(Lu2(3,3))+exp(Lu2(4,4)));% depolarization
A=(exp(Lu2(2,2))-exp(Lu2(3,3)))/(exp(Lu2(2,2))+exp(Lu2(3,3)));% linear depolarization anisotropy
end