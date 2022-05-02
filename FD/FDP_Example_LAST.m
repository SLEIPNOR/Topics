clear all 
close all
clc

%%
%define matrices
A = -3;
B = 2;
Bf = 1;
Bd = [5 -5];
C = [4; -6];
Df = [7; 2];
Dd = [-5 0;4 3];
[U,sigma,V] = svd(Df);
U1 = U(:,1);
U2 = U(:,2);
sigma1 = sigma(1);
% a = inv(S1);
% a1 = S1^(-1);
% a2 = 1/S1;

%provide dimensions of matrices
n = 1;
nu = 1;
nf = 1;
nd = 2;
ny = 2;

%%
%define Df_dagger and Df_perp
Df_dagger = V*inv(sigma1)*U1';
Df_perp = U2';
%%
%deifne A1,B1,C1,C2,D1,D2
A1 = A - Bf*Df_dagger*C;
B1 = Bd - Bf*Df_dagger*Dd;
C1 = Df_dagger*C;
C2 = Df_perp*C;
D1 = Df_dagger*Dd;
D2 = Df_perp*Dd;
matrices = [A1 B1;C1 D1;C2 D2];

%%
%main calculation: use LMI solver to obtain P, Z, R, S, and gamma
setlmis([]);
g2 = lmivar(1,[1,1]);
P = lmivar(1,[n,1]);
Z = lmivar(2,[n,ny-nf]);
S = lmivar(2,[nf,ny-nf]);
L = newlmi;
lmiterm([L,1,1,P],1,A1,'s'); lmiterm([L,1,1,Z],1,C2,'s'); %%should not be A, should be A1
lmiterm([L,1,2,P],1,B1); lmiterm([L,1,2,Z],1,D2);
lmiterm([L,1,3,0],C1'); lmiterm([L,1,3,-S],C2',1);
lmiterm([L,2,2,g2],-0.5,1,'s');
lmiterm([L,2,3,0],D1'); lmiterm([L,2,3,-S],D2',1);
lmiterm([L,3,3,0],-1);
LMI = getlmis;
[g2,x] = mincx(LMI,eye(decnbr(LMI),1));
P = dec2mat(LMI,x,P);
Z = dec2mat(LMI,x,Z);
S = dec2mat(LMI,x,S);
gamma = sqrt(g2);
R = inv(P)*Z;
%get P, Z, R, S, and gamma
a = [P Z R S];
% gamma = sqrt(g2);
%%
%calculate L and H
L = -Bf*Df_dagger + R*Df_perp;
H = Df_dagger + S*Df_perp;

eigs(A+L*C)
% Q1 = obsv(A+L*C,H*C);
% rank(Q1)
% H_new = inv(H*Df)*H;
% L_new = L - (Bf+L*Df)*H_new;
G = -(H*C)'*(H*C);
[X,K1,L1] = icare((A+L*C)',[],[],[],[],[],G);
L_new = L - X*(H*C)'*H;
eigs(A+L_new*C)

%%
%matrices
D = [0; 0];
G = [A B;C D]
G0 = G;
Gf = [A Bf;C Df]
Gd = [A Bd;C Dd]
F = [A+L_new*C L_new; H*C H]



