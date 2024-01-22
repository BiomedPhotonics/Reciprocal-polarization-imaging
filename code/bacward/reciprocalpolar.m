function [M_delta,M_R,MD,trad,psi,delt,theta,A] = reciprocalpolar(M)
Q = diag([1,1,-1,1],0);
G =diag([1,-1,-1,-1],0);

%% Symmetric matrix QM
QM = Q * M;
QM = (QM + QM')/2;
QMG=QM*G;

%% MD
[V,D]=eig(QMG);
W1 = real(V);
U1 = real(D);
U=diag(U1);
W=W1 /diag(W1(1,:),0);
J= find(U>=0); 
multiply=W.*W;
indictor = 1-sum(multiply(2:4, J));
% There should be one unique eigenvector.The eigenvector should first come from positive eigenvalues
[~,i] = max(indictor);
Dv = W(2:4,J(i));
D=norm(Dv,'fro');
if (D>1 || D<1e-6)
    D = 1;
    MD=diag([1,1,1,1],0);
    d0=1;
    MDinv = MD;
else
    Dv0 = Dv/D;
    I = eye(3);
    mD = (sqrt(1 - D^2) * I) + ((1 - sqrt(1 - D^2)) * Dv0'.*Dv0);
    MD =  [1, Dv'; Dv, mD];
    d0 = U(J(i))/(1-D^2);
    MDinv = G * MD * G / (1-D^2); 
end

%% N
N = MDinv * QM * MDinv;
N1=N(2:4,2:4);
%% The order of W3 and U3 is arbitary. We need to reorder it. It means the decomposition is not unique!
[V1, D1] = eig(N1);
W3=diag(D1);
U3 = V1';
rmin=100;
cur_ind=[];
cur_sign=[];
minW= min(W3);
%% Try out all permutation and +- to find the least R
    for i =perms(1:3)'
        if( W3(i(3))==minW) %d3<0
            for j1=[-1,1]
                for j2=[-1,1]
                    for j3=[-1,1]
                        tmp=U3(i,:);
                        tmp1(1,:) = j1 * tmp(1,:) ;
                        tmp1(2,:) = j2 * tmp(2,:);
                        tmp1(3,:) = j3 * tmp(3,:);
                        R=(trace(tmp1)+1)/2-1;
                        if(det(tmp1)>0  && tmp1(3,3)<=0 )%det(tmp1)>0,and Linear retardance>pi/2 (optional)
                            p = R;
                            if(abs(p)-rmin<1e-5)
                                cur_ind = i;
                                cur_sign = [j1, j2, j3];
                                Rmin = abs(p);
                                rmin=Rmin;
                            end

                        end
                    end
                end
            end
        end
    end
W4 = W3(cur_ind);
U4 = U3(cur_ind,:);
for i=1:3
    U4(i,:) = U4(i,:) * cur_sign(i)  ;
end
Z=[0 0 0];
M_R=[1, Z; Z', U4] ;
M_delta = diag([d0, W4(1), -W4(2), W4(3)],0);

%% Parameter
trad=1-(trace( abs(M_delta))-d0)/(3*d0); % depolarization
psi=1/2*atan2((M_R(2,3)-M_R(3,2)),(M_R(2,2)+M_R(3,3)));% optical rotation
delt =acos((sqrt((M_R(2,2)+M_R(3,3))^2+(M_R(2,3)-M_R(3,2))^2))-1);% linear retardance
theta=atan2((M_R(4,2)),(-M_R(4,3)))*90/pi;% orientation angle
A=(abs(M_delta(2,2))-abs(M_delta(3,3)))/(abs(M_delta(2,2))+abs(M_delta(3,3)));% linear depolarization anisotropy
end