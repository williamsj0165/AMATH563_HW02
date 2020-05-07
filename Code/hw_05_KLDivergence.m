%% 5. KL Divergence

X = [
    hare'
    lynx'
    ];

m = length(X);

% Does it matter if these are different?
x_hare = 1 : 150;
x_lynx = 1 : 100;

f_hare = hist(X(1,:),x_hare)+0.01;
f_lynx = hist(X(2,:),x_lynx)+0.01;

f_hare = f_hare / trapz(x_hare,f_hare);
f_lynx = f_lynx / trapz(x_lynx,f_lynx);

%% 1. DMD

g1_hare = hist(x(1,1:m),x_hare)+0.01;
g1_lynx = hist(x(2,1:m),x_lynx)+0.01;

g1_hare = g1_hare / trapz(x_hare,g1_hare);
g1_lynx = g1_lynx / trapz(x_lynx,g1_lynx);

Int1_hare = f_hare.*log(f_hare./g1_hare);
Int1_lynx = f_lynx.*log(f_lynx./g1_lynx);

I1_hare = trapz(x_hare,Int1_hare)
I1_lynx = trapz(x_lynx,Int1_lynx)


%% 2. DMD w/ Time-Delay Embedding

g2_hare = hist(h(1,1:m),x_hare)+0.01;
g2_lynx = hist(h(1,1:m),x_lynx)+0.01;

g2_hare = g2_hare / trapz(x_hare,g2_hare);
g2_lynx = g2_lynx / trapz(x_lynx,g2_lynx);

Int2_hare = f_hare.*log(f_hare./g2_hare);
Int2_lynx = f_lynx.*log(f_lynx./g2_lynx);

I2_hare = trapz(x_hare,Int2_hare)
I2_lynx = trapz(x_lynx,Int2_lynx)


%% 3. Lotka-Voltera

g3_hare = hist(y_KL(1,:),x_hare)+0.01;
g3_lynx = hist(y_KL(2,:),x_lynx)+0.01;

g3_hare = g3_hare / trapz(x_hare,g3_hare);
g3_lynx = g3_lynx / trapz(x_lynx,g3_lynx);

Int3_hare = f_hare.*log(f_hare./g3_hare);
Int3_lynx = f_lynx.*log(f_lynx./g3_lynx);

I3_hare = trapz(x_hare,Int3_hare)
I3_lynx = trapz(x_lynx,Int3_lynx)


%% 4. SINDy

g4_hare = hist(x_pred_stl(1:30,1),x_hare)+0.01;
g4_lynx = hist(x_pred_stl(1:30,2),x_lynx)+0.01;

g4_hare = g4_hare / trapz(x_hare,g4_hare);
g4_lynx = g4_lynx / trapz(x_lynx,g4_lynx);

Int4_hare = f_hare.*log(f_hare./g4_hare);
Int4_lynx = f_lynx.*log(f_lynx./g4_lynx);

I4_hare = trapz(x_hare,Int4_hare)
I4_lynx = trapz(x_lynx,Int4_lynx)


%% 5.

KL_Divs = [
    I1_hare
    I1_lynx
    I2_hare
    I2_lynx
    I3_hare
    I3_lynx
    I4_hare
    I4_lynx
    ]







% 
% 
% 
% 
% Int1 = Int1=f.*log(f./g1);
% 
% I1 = trapz(x,Int1)
