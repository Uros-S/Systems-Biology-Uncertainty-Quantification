function dYdt = autoactivating_feedback_loop(t,Y,parameters)

%parameters setting
n = parameters(1);
K = parameters(2);
a = parameters(3);

% Hindmarsh-Rose model equations
x = Y;
dYdt = K + a*(x^n)/(1+x^n) - x;

end