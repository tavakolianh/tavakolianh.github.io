var 
    c
    ii
    pi
    y 
    k 
    m 
    tau 
    b 
    a;

varexo 
    eps_a;

parameters
    betta
    siggma
    etta
    chi
    n
    deltta
    allpha
    thetta
    rho_a;

    betta   = 0.99;
    siggma  = 1.5;
    etta    = 2;
    chi     = 2;
    n       = 0;
    deltta  = 0.01;
    allpha  = 1/3;
    thetta  = 0.1;
    rho_a   = 0.8;

model;
// Demand for money
    chi * m^(-etta)/ c^(-siggma) = ii / (1 + ii);
// Euler Equation
    (c(+1) / c)^(siggma) = betta * (allpha * a * k(-1)^(allpha - 1) + 1 - deltta);
// Production function
    y = a * (k(-1)/(1 + n))^allpha;
// Fisher Equation
   allpha * a * k(-1)^(allpha - 1) + 1 - deltta = (1 + ii)/ (1 + pi(+1)); 
//Household Budget
   y + tau + (1 - deltta) * k(-1)/(1 + n) + ((1 + ii(-1)) * b(-1) + m(-1))/((1 + pi) * (1 + n))
   = c + k + m + b ;
// government debt
   b = 0; 
// Money Growth
   m *(1+ pi) / m(-1) = 1 + thetta;
// Government Budget
   tau + ((1 + ii(-1)) * b(-1) + m(-1))/((1 + pi) * (1 + n)) = m + b;
// Technology
   log(a) = rho_a * log(a(-1)) + eps_a;
end;


steady_state_model;
    a = 1;
    b = 0;
    k = (allpha * betta / (1 + betta*(deltta - 1)))^(1/ (1 - allpha));
    y = a * k^allpha;
    c = y - deltta * k;
    pi = thetta;
    ii = (allpha * a * k^(allpha - 1) + 1 - deltta)*(1 + pi) - 1;
    m = (ii * c^(-siggma) / ((1 + ii) * chi))^(-1/etta);
    tau = pi * m / (1 +pi);    
end;


steady;
check; 

shocks;
var eps_a; stderr 0.01;
end;

stoch_simul(periods=1000000, order=2);



