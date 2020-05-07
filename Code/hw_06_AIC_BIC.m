%% 6. AIC and BIC Scores



%% lambda = 0.0025
k_4a_hare = 16;
k_4a_lynx = 16;

L_4a_hare = max(mle(g4_hare));
L_4a_lynx = max(mle(g4_lynx));

AIC_4a_hare = 2 * k_4a_hare - 2 * log(L_4a_hare)
AIC_4a_lynx = 2 * k_4a_lynx - 2 * log(L_4a_lynx)

BIC_4a_hare = log(length(year)) * k_4a_hare - 2 * log(L_4a_hare)
BIC_4a_lynx = log(length(year)) * k_4a_lynx - 2 * log(L_4a_lynx)


%% lambda = 0.025
k_4b_hare = 15;
k_4b_lynx = 15;

L_4b_hare = max(mle(g4_hare));
L_4b_lynx = max(mle(g4_lynx));


AIC_4b_hare = 2 * k_4b_hare - 2 * log(L_4b_hare)
AIC_4b_lynx = 2 * k_4b_lynx - 2 * log(L_4b_lynx)

BIC_4b_hare = log(length(year)) * k_4b_hare - 2 * log(L_4b_hare)
BIC_4b_lynx = log(length(year)) * k_4b_lynx - 2 * log(L_4b_lynx)


%% lambda = 2.5

k_4c_hare = 6;
k_4c_lynx = 4;

L_4c_hare = max(mle(g4_hare));
L_4c_lynx = max(mle(g4_lynx));

AIC_4c_hare = 2 * k_4c_hare - 2 * log(L_4c_hare)
AIC_4c_lynx = 2 * k_4c_lynx - 2 * log(L_4c_lynx)

BIC_4c_hare = log(length(year)) * k_4c_hare - 2 * log(L_4c_hare)
BIC_4c_lynx = log(length(year)) * k_4c_lynx - 2 * log(L_4c_lynx)

