clear;
A=-2;B=1;Bf=1;Bd=[4,-4];C=[4;-5];Df=[6;1];Dd=[-4,0;3,2];
n = 1;
nu = 1;
nf = 1;
nd = 2;
ny = 2;
[U,S,V] = svd(Df);
Dfd = V*(S(1)^-1)*U(:,1)';
Dfv = U(:,2)';
A1=A-Bf*Dfd*C;
B1=Bd-Bf*Dfd*Dd;
C1=Dfd*C;
C2=Dfv*C;
D1=Dfd*Dd;
D2=Dfv*Dd;
%% main calculation: use LMI solver to obtain P, Z, R, S, and gamma
setlmis([]);
g2 = lmivar(1,[1,1]);
P = lmivar(1,[n,1]);
Z = lmivar(2,[n,ny-nf]);
S = lmivar(2,[nf,ny-nf]);
L = newlmi;
lmiterm([L,1,1,P],1,A1,'s'); lmiterm([L,1,1,Z],1,C2,'s'); 
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
R = (P^-1)*Z;
%% L&H
H=Dfd+S*Dfv;
L=-Bf*Dfd+R*Dfv;
%% stable evaluation
if eig(A+L*C)<0
else
    [X,Kp,Lp] = icare((A+L*C)',[],[],[],[],[],-(H*C)'*H*C);
    LN=L-X*(H*C)'*H;
end
%% matrices
D = [0; 0];
G = [A B;C D]
G0 = G;
Gf = [A Bf;C Df]
Gd = [A Bd;C Dd]
F = [A+LN*C LN; H*C H]