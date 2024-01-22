%% Lu-chipman polar,Differential and Reciprocal polar decomposition in backward geometry
%% Measure Mueller matrix (see Supplementary documents in https://doi.org/10.48550/arXiv.2305.12053)

clc;clear ;
%% Import File
sample_path = 'C:\Users\dell\Desktop\PRX\github\Original polarization image\target-backward\';
Readfiles = dir(fullfile(sample_path,'*.tiff'));
%% Separate the four polarization directions of a polarization camera:  0¡ã, 45¡ã, 90¡ã, and 135¡ãdirections 
wx=1;
for i=1:8
    s= double(imread(strcat(sample_path,Readfiles(i).name)));
    for x=1:2:2047
        for y=1:2:2447
            a11((x+1)/2,(y+1)/2)=s(x,y);
        end
    end
    H = fspecial('average', [wx, wx]);
    A11(:,:,i) = imfilter(a11, H);
    for x=2:2:2048
        for y=1:2:2447
            a21(x/2,(y+1)/2)=s(x,y);
        end
    end
    A21(:,:,i) = imfilter(a21, H);
    for x=1:2:2047
        for y=2:2:2448
            a12((x+1)/2,y/2)=s(x,y);
        end
    end
    A12(:,:,i) = imfilter(a12, H);

    for x=2:2:2048
        for y=2:2:2448
            a22(x/2,y/2)=s(x,y);
        end
    end
    A22(:,:,i) = imfilter(a22, H);
end
b11=A11(:,:,1);a11=A11(:,:,2);d11=A11(:,:,3);c11=A11(:,:,4);f11=A11(:,:,5);e11=A11(:,:,6);k11=A11(:,:,7);h11=A11(:,:,8);
b12=A12(:,:,1);a12=A12(:,:,2);d12=A12(:,:,3);c12=A12(:,:,4);f12=A12(:,:,5);e12=A12(:,:,6);k12=A12(:,:,7);h12=A12(:,:,8);
b21=A21(:,:,1);a21=A21(:,:,2);d21=A21(:,:,3);c21=A21(:,:,4);f21=A21(:,:,5);e21=A21(:,:,6);k21=A21(:,:,7);h21=A21(:,:,8);
b22=A22(:,:,1);a22=A22(:,:,2);d22=A22(:,:,3);c22=A22(:,:,4);f22=A22(:,:,5);e22=A22(:,:,6);k22=A22(:,:,7);h22=A22(:,:,8);

%% Mueller matrix 
M11=1/8*(a11+a21+a22+a12+b11+b21+b22+b12);
M21=1/2*(b22-b11);
M31=1/2*(a12-a21);
M41=1/4*(b12-b21-a22+a11);
M12=1/8*(c11+c21+c22+c12+d11+d21+d22+d12);
M22=1/2*(d22-d11);
M32=1/2*(c12-c21);
M42=1/4*(d12-d21-c22+c11);
M13=1/8*(e11+e21+e22+e12+f11+f21+f22+f12);
M23=1/2*(f22-f11);
M33=1/2*(e12-e21);
M43=1/4*(f12-f21-e22+e11);
M14=1/8*(h11+h21+h22+h12+k11+k21+k22+k12);
M24=1/2*(k22-k11);
M34=1/2*(h12-h21);
M44=1/4*(k12-k21-h22+h11);
%% Mueller matrix of optical element
Rin =[1.0000    0.4967    0.7353    0.7400;
    0.9712   -0.4824    0.2027    0.2577;
    0.0040    0.0028   -0.6757   -0.2150;
    0.0027   -0.0036    0.1683   -0.6524];% R * Sin

M2=[1.0000   -0.0532   -0.0268   -0.0090;
    -0.0524    0.9956   -0.0247   -0.0205;
    0.0383   -0.0018   -0.9419    0.3706;
    -0.0023   -0.0725   -0.3802   -0.8836];

T = [1.0000   -0.3952    0.0149   -0.0002;
    -0.3776    0.9969    0.0056   -0.0230;
    0.0248   -0.0497    0.9066    0.1027;
    0.0318   -0.0009   -0.0447    0.8740];
%% ROI
n1 = 475;
m1 =380;
n2= 730;
m2=555;
%% Lu-chipman polar decomposition in backward geometry
for i=n1:n2
    for j=m1:m2
        M(1,1)=M11(i ,j );
        M(1,2)=M12(i ,j );
        M(2,1)=M21(i ,j );
        M(2,2)=M22(i ,j );
        M(1,3)=M13(i ,j );
        M(1,4)=M14(i ,j );
        M(2,3)=M23(i ,j );
        M(2,4)=M24(i ,j );
        M(3,1)=M31(i ,j );
        M(3,2)=M32(i ,j );
        M(3,3)=M33(i ,j );
        M(3,4)=M34(i ,j );
        M(4,1)=M41(i ,j );
        M(4,2)=M42(i ,j );
        M(4,3)=M43(i ,j );
        M(4,4)=M44(i ,j );
        M =inv(T) * inv(M2)  * M * inv(Rin) ;
        M=M/M(1,1);
        M=diag([1,1,-1,-1],0)*M;
        [M_delta,M_R,MD,trad,psi,delt,theta,A] = Luchipman(M);
        a1(i-n1+1,j-m1+1)=A;% linear depolarization anisotropy
        trad1(i-n1+1,j-m1+1)=trad;% depolarization
        psi1(i-n1+1,j-m1+1)=psi;% optical rotation
        delt1(i-n1+1,j-m1+1)=delt;% linear retardance
        theta1(i-n1+1,j-m1+1)=theta;% orientation angle
    end
end

%% differential decomposition in backward geometry
for i=n1:n2
    for j=m1:m2
        M(1,1)=M11(i ,j );
        M(1,2)=M12(i ,j );
        M(2,1)=M21(i ,j );
        M(2,2)=M22(i ,j );
        M(1,3)=M13(i ,j );
        M(1,4)=M14(i ,j );
        M(2,3)=M23(i ,j );
        M(2,4)=M24(i ,j );
        M(3,1)=M31(i ,j );
        M(3,2)=M32(i ,j );
        M(3,3)=M33(i ,j );
        M(3,4)=M34(i ,j );
        M(4,1)=M41(i ,j );
        M(4,2)=M42(i ,j );
        M(4,3)=M43(i ,j );
        M(4,4)=M44(i ,j );
            M =inv(T) * inv(M2)  * M * inv(Rin) ;
            M=M/M(1,1);
            M=diag([1,1,-1,-1],0)*M;
            [trad,psi,delt,theta,A] = differential(M);       
            a2(i-n1+1,j-m1+1)=A;% linear depolarization anisotropy 
            trad2(i-n1+1,j-m1+1)=trad;% depolarization
            psi2(i-n1+1,j-m1+1)=psi;% optical rotation
            delt2(i-n1+1,j-m1+1)=delt;% linear retardance
            theta2(i-n1+1,j-m1+1)=theta;% orientation angle
    end
end

%% Reciprocal polar decomposition in backward geometry
for i=n1:n2
    for j=m1:m2
        M(1,1)=M11(i ,j );
        M(1,2)=M12(i ,j );
        M(2,1)=M21(i ,j );
        M(2,2)=M22(i ,j );
        M(1,3)=M13(i ,j );
        M(1,4)=M14(i ,j );
        M(2,3)=M23(i ,j );
        M(2,4)=M24(i ,j );
        M(3,1)=M31(i ,j );
        M(3,2)=M32(i ,j );
        M(3,3)=M33(i ,j );
        M(3,4)=M34(i ,j );
        M(4,1)=M41(i ,j );
        M(4,2)=M42(i ,j );
        M(4,3)=M43(i ,j );
        M(4,4)=M44(i ,j );
        M =inv(T) * inv(M2)  * M * inv(Rin) ;
        M=M/M(1,1);
        [M_delta,M_R,MD,trad,psi,delt,theta,A] = reciprocalpolar(M);
        a3(i-n1+1,j-m1+1)=A;% linear depolarization anisotropy
        trad3(i-n1+1,j-m1+1)=trad;% depolarization
        psi3(i-n1+1,j-m1+1)=psi;% optical rotation
        delt3(i-n1+1,j-m1+1)=delt;% linear retardance
        theta3(i-n1+1,j-m1+1)=theta;% orientation angle
    end
end