function dYdt = Epileptor(t,Y,parameters)

dYdt = zeros(3,1);

%parameters setting
r_1 = parameters(1);
r_2 = parameters(2);
I_1 = parameters(3);
I_2 = parameters(4);
tau_0 = parameters(5);
m = parameters(6); 
tau_2 = parameters(7);
gamma = parameters(8); 

% Epileptor model equations
x_1 = Y(1); x_2 = Y(2); x_3 = Y(3); x_4 = Y(4); x_5 = Y(5); u = Y(6);

if x_1 < 0
    f_1 = x_1^3-3*x_1^2;
else
    f_1 = (x_4-0.6*(x_3-4)^2-m)*x_1;
end

if x_4 < -0.25
    f_2 = 0;
else
    f_2 = 6*(x_4+0.25);
end

dYdt(1) = x_2 - f_1 - x_3 + I_1;
dYdt(2) = r_2 - 5*x_1^2 - x_2;
dYdt(3) = 1/tau_0*(4*(x_1 - r_1) - x_3);
dYdt(4) = -x_5 + x_4 -x_4^3 + I_2 + 2*u - 0.3*(x_3 - 3.5);
dYdt(5) = 1/tau_2*(-x_5 + f_2);
dYdt(6) = -gamma*(u-0.1*x_1);

end