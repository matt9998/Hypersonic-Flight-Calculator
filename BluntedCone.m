function[pinf,Tinf,rhoinf,p2,h2,rho2,T2,u2,eps,iterationHist,CL,CD,qstag,qref] = BluntedCone(z,M1,x,dc,rn,rc)

% 1b) Freestream Conditions 
r0 = 6378.15; %Earth's radius, in km
g0 = 9.806; %Gravitational constant for Earth
mw0 = 28.9644; %Molecular Weight
RBar = 8.314E3; %Gas Cosntant
b = 3.31E-7; %Accounts for change in molecular weight for regions 7+
h = (z*r0)/(r0+z); %Obtaining geopotential altitude,used for calculation if in regions 1-6
A = [0,11.0102,-6.5,288.15,101325;11,20.0631,0,216.65,22631.95;20,32.1619,1,216.65,5474.79;32,47.3501,2.8,228.65,868.01;47,51.4125,0,270.65,110.9;51,71.8020,-2.8,270.65,66.94;71,86,-2,214.65,3.956]; %Array of atmospheric table, regions 0-6
B = [86,1.6481,186.945,186.946,28.9644,.34418;100,5,210.65,210.02,28.88,.29073;110,10,260.65,257,28.56,.006801;120,20,360.65,349.49,28.08,.002247]; %Array of atmosperic table, regions 7-10
if h <= 11
    pinf = A(1,5)*((A(1,4)/(A(1,4)+A(1,3)*(h-A(1,1))))^((g0*mw0)/(8.314*A(1,3)))); %Eq 8 for pressure, non-isothermal layer
    Tinf = A(1,4)+(A(1,3)*(h-A(1,1))); %Eq 7, linearly varying temperature
    rhoinf = (pinf*mw0)/(RBar*Tinf); %Eq 10, Computes freestream density based on pinf and Tinf
elseif h > 11 && h <= 20
    pinf = A(2,5)*exp((-1*g0*mw0*(h-A(2,1)))/(8.314*A(2,4))); %Eq 9 for pressure, isothermal layer
    Tinf = A(2,4); %No change in temperature across layer
    rhoinf = (pinf*mw0)/(RBar*Tinf); %Eq 10, Computes freestream density based on pinf and Tinf
elseif h > 20 && h <= 32
    pinf = A(3,5)*((A(3,4)/(A(3,4)+A(3,3)*(h-A(3,1))))^((g0*mw0)/(8.314*A(3,3)))); %Eq 8 for pressure, non-isothermal layer
    Tinf = A(3,4)+(A(3,3)*(h-A(3,1))); %Eq 7, linearly varying temperature
    rhoinf = (pinf*mw0)/(RBar*Tinf); %Eq 10, Computes freestream density based on pinf and Tinf
elseif h > 32 && h <= 47
    pinf = A(4,5)*((A(4,4)/(A(4,4)+A(4,3)*(h-A(4,1))))^((g0*mw0)/(8.314*A(4,3)))); %Eq 8 for pressure, non-isothermal layer
    Tinf = A(4,4)+(A(4,3)*(h-A(4,1))); %Eq 7, linearly varying temperature 
    rhoinf = (pinf*mw0)/(RBar*Tinf); %Eq 10, Computes freestream density based on pinf and Tinf
elseif h > 47 && h <= 51
    pinf = A(5,5)*exp((-1*g0*mw0*(h-A(5,1)))/(8.314*A(5,4))); %Eq 9 for pressure, isothermal layer
    Tinf = A(5,4); %No change in temperature across layer
    rhoinf = (pinf*mw0)/(RBar*Tinf); %Eq 10, Computes freestream density based on pinf and Tinf
elseif h > 51 && h <= 71
    pinf = A(6,5)*((A(6,4)/(A(6,4)+A(6,3)*(h-A(6,1))))^((g0*mw0)/(8.314*A(6,3)))); %Eq 8 for pressure, non-isothermal layer
    Tinf = A(6,4)+(A(6,3)*(h-A(6,1))); %Eq 7, linearly varying temperature
    rhoinf = (pinf*mw0)/(RBar*Tinf); %Eq 10, Computes freestream density based on pinf and Tinf
elseif h > 71 && h <= 84.852
    pinf = A(7,5)*((A(7,4)/(A(7,4)+A(7,3)*(h-A(7,1))))^((g0*mw0)/(8.314*A(7,3)))); %Eq 8 for pressure, non-isothermal layer
    Tinf = A(7,4)+(A(7,3)*(h-A(7,1))); %Eq 7, linearly varying temperature
    rhoinf = (pinf*mw0)/(RBar*Tinf); %Eq 10, Computes freestream density based on pinf and Tinf
elseif z > 86 && z <= 100
    R = RBar/(B(1,5)); %Computes gas constant for layer
    Ti = (B(1,3)*B(1,5))/28.9664; %Computes temperature of bottom of layer
    rhoi = (B(1,6)*B(1,5))/(RBar*Ti); %Computes density at bottom of layer
    POW = -1*((g0/(R*B(1,2)*.001))*(1+b*((B(1,3)/(B(1,2)*.001))-(B(1,1)*1000)))); %Computes power term for pressure
    Term1 = (B(1,3)+B(1,2)*(z-B(1,1)))/B(1,3); %Part of pressure computation
    Term2 = exp(((g0*b)/(R*B(1,2)*.001))*(z-B(1,1))); %Part of pressure computation
    pinf = B(1,6)*(Term1^POW)*Term2; %Computes freestream pressure
    POWD = -1*((g0/(R*B(1,2)*.001))*(((R*B(1,2)*.001)/g0)+1+(b*((B(1,3)/(B(1,2)*.001))-(B(1,1)*1000))))); %Exponent term for density
    rhoinf = rhoi*(Term1^POWD)*Term2; %Computes freestream density
    Tinf = (pinf*B(1,5))/(RBar*rhoinf); %Computes freestream temperature
elseif z > 100 && z <= 110
    R = RBar/(B(2,5)); %Computes gas constant for layer
    Ti = (B(2,3)*B(2,5))/28.9664; %Computes temperature of bottom of layer
    rhoi = (B(2,6)*B(2,5))/(RBar*Ti); %Computes density at bottom of layer
    POW = -1*((g0/(R*B(2,2)*.001))*(1+b*((B(2,3)/(B(2,2)*.001))-(B(2,1)*1000)))); %Computes power term for pressure
    Term1 = (B(2,3)+B(2,2)*(z-B(2,1)))/B(2,3); %Part of pressure computation
    Term2 = exp(((g0*b)/(R*B(2,2)*.001))*(z-B(2,1))); %Part of pressure computation
    pinf = B(2,6)*(Term1^POW)*Term2; %Computes freestream pressure
    POWD = -1*((g0/(R*B(2,2)*.001))*(((R*B(2,2)*.001)/g0)+1+(b*((B(2,3)/(B(2,2)*.001))-(B(2,1)*1000))))); %Exponent term for density
    rhoinf = rhoi*(Term1^POWD)*Term2; %Computes freestream density
    Tinf = (pinf*B(2,5))/(RBar*rhoinf); %Computes freestream temperature
elseif z > 110 && z <= 120
    R = RBar/(B(3,5)); %Computes gas constant for layer
    Ti = (B(3,3)*B(3,5))/28.9664; %Computes temperature of bottom of layer
    rhoi = (B(3,6)*B(3,5))/(RBar*Ti); %Computes density at bottom of layer
    POW = -1*((g0/(R*B(3,2)*.001))*(1+b*((B(3,3)/(B(3,2)*.001))-(B(3,1)*1000)))); %Computes power term for pressure
    Term1 = (B(3,3)+B(3,2)*(z-B(3,1)))/B(3,3); %Part of pressure computation
    Term2 = exp(((g0*b)/(R*B(3,2)*.001))*(z-B(3,1))); %Part of pressure computation
    pinf = B(3,6)*(Term1^POW)*Term2; %Computes freestream pressure
    POWD = -1*((g0/(R*B(3,2)*.001))*(((R*B(3,2)*.001)/g0)+1+(b*((B(3,3)/(B(3,2)*.001))-(B(3,1)*1000))))); %Exponent term for density
    rhoinf = rhoi*(Term1^POWD)*Term2; %Computes freestream density
    Tinf = (pinf*B(3,5))/(RBar*rhoinf); %Computes freestream temperature    
elseif z > 120 && z <= 150
    R = RBar/(B(4,5)); %Computes gas constant for layer
    Ti = (B(4,3)*B(4,5))/28.9664; %Computes temperature of bottom of layer
    rhoi = (B(4,6)*B(4,5))/(RBar*Ti); %Computes density at bottom of layer
    POW = -1*((g0/(R*B(4,2)*.001))*(1+b*((B(4,3)/(B(4,2)*.001))-(B(4,1)*1000)))); %Computes power term for pressure
    Term1 = (B(4,3)+B(4,2)*(z-B(4,1)))/B(4,3); %Part of pressure computation
    Term2 = exp(((g0*b)/(R*B(4,2)*.001))*(z-B(4,1))); %Part of pressure computation
    pinf = B(4,6)*(Term1^POW)*Term2; %Computes freestream pressure
    POWD = -1*((g0/(R*B(4,2)*.001))*(((R*B(4,2)*.001)/g0)+1+(b*((B(4,3)/(B(4,2)*.001))-(B(4,1)*1000))))); %Exponent term for density
    rhoinf = rhoi*(Term1^POWD)*Term2; %Computes freestream density
    Tinf = (pinf*B(4,5))/(RBar*rhoinf); %Computes freestream temperature
end

% 1c) Computing Conditions Downstream of Normal Shock
Cp = 1004.5; %Defines Cp upstream
if M1 < 3 %Using CPG Assumption/Shock equations for low speeds under Mach 3
    h1 = Cp*Tinf; %Computes upstream enthalpy
    a = sqrt(1.4*287*Tinf); %Defines speed of sound, useful for calculating u2
    rho2 = rhoinf*((2.4*M1^2)/(2+(.4*M1^2))); %Eq 17 for downstream density
    u1 = M1*a; %Gives upstream velocity in m/s
    u2 = (u1*rhoinf)/rho2; %Eq 17 for downstream velocity
    p2 = pinf*(1+(1.167*(M1^2-1))); %Eq 18 for downstream pressure
    T2 = Tinf*((1.167*(M1^2-1))*((2+(.4*M1^2))/(2.4*M1^2))); %Eq 19 for downstream temperature
    h2 = (T2*h1)/Tinf; %Eq 19 for downstream enthalpy
elseif M1 >= 3 %Using equilibrium (iterative) solution procedure for chemically reacting gas
    h1 = Cp*Tinf; %Computes upstream enthalpy
    a = sqrt(1.4*287*Tinf); %Defines speed of sound, useful for calculating u2    
    u1 = M1*a; %Gives upstream velocity in m/s
    eps = .1; %Initial guess for epsilon
    Deps = 1; %Arbitrary starting value dor Deps, which will be driven to 0
    func = 1; %Sets arbitary starting value for h2 - h2j, which will be driven to 0 in a nested while loop, used below
    num = 0; %Starts iteration counter
    G = [1.4	0	0	0	0	0	0	0	0.0	0	0.0000
1.42598	0.000918	-0.092209	-0.002226	0.019772	-0.0366	-0.077469	0.043878	-15.0	-1.000	-1.0400
1.64689	-0.062133	-0.334994	0.063612	-0.038332	-0.014468	0.073421	-0.002442	-15.0	-1.000	-1.3600
1.48558	-0.453562	0.152096	0.30335	-0.459282	0.448395	0.220546	-0.292293	-10.0	-1.000	-1.6000
1.4	0	0	0	0	0	0	0	0.0	0.000	0.0000
1.42176	-0.000366	-0.083614	0.000675	0.005272	-0.115853	-0.007363	0.146179	-20.0	-1.000	-0.8600
1.74436	-0.035354	-0.415045	0.061921	0.018536	0.043582	0.044353	-0.04975	-20.0	-1.040	-1.3360
1.49674	-0.021583	-0.197008	0.030886	-0.157738	-0.009158	0.123213	-0.006553	-10.0	-1.050	-1.8950
1.10421	-0.033664	0.031768	0.024335	-0.176802	-0.017456	0.080373	0.002511	-15.0	-1.080	-2.6500
1.4	0	0	0	0	0	0	0	0.0	0.000	0.0000
1.47003	0.007939	-0.244205	-0.025607	0.872248	0.049452	-0.764158	0.000147	-20.0	-1.000	-0.7420
3.18652	0.13793	-1.89529	-0.10349	-2.14572	-0.272717	2.06586	0.223046	-15.0	-1.000	-1.0410
1.63963	-0.001004	-0.303549	0.016464	-0.852169	-0.101237	0.503123	0.04358	-10.0	-1.000	-1.5440
1.55889	0.055932	-0.211764	-0.023548	-0.549041	-0.101758	0.276732	0.046031	-15.0	-1.000	-2.2500
]; %Copying gamma tilda regions from Excel spreadsheet
    D = [0.27407	0	1.00082	0	0	0	0	0	0	0	0.000	0
0.235869	-0.043304	1.17619	0.046498	-0.143721	-1.3767	0.160465	1.08988	-0.083489	-0.217748	-10.000	-1.78
0.281611	0.001267	0.990406	0	0	0	0	0	0	0	0.000	0
0.457643	-0.034272	0.819119	0.046471	0	-0.073233	-0.169816	0.043264	0.111854	0	-15.000	-1.28
1.04172	0.041961	0.412752	-0.009329	0	-0.434074	-0.196914	0.264883	0.100599	0	-15.000	-1.778
0.418298	-0.2521	0.784048	0.144576	0	-2.00015	-0.639022	0.716053	0.206457	0	-10.000	-2.4
2.72964	0.003725	0.938851	-0.01192	0	0.682406	0.089153	-0.646541	-0.070769	0	-20.000	-0.82
2.50246	-0.042827	1.12924	0.041517	0	1.72067	0.268008	-1.25038	-0.179711	0	-20.000	-1.33
2.44531	-0.047722	1.00488	0.034349	0	1.95893	0.316244	-1.012	-0.151561	0	-20.000	-1.88
2.50342	0.026825	0.83886	-0.009819	0	3.58284	0.533853	-1.36147	-0.195436	0	-20.000	-2.47
]; %Copying temperature regions from Excel spreadsheet
    p2 = pinf + (rhoinf*(u1^2)*(1-eps)); %Eq 14, Computes value for p2
    rhostart = rhoinf*(1/eps); %Initial guess for rho2, used to locate gamme tilda region
    Y = log10(rhostart/1.292); %Initial value for Y
    X = log10(p2/101300); %Initial value for X
    Z = X - Y; %Initial value for Z
    if Y > -.5 && Z <= .3 %This if statement obtains the region in which to perform gamma tilda calculations
        n=1;
    elseif Y > -.5 && Z > .3 && Z <= 1.15
        n=2;
    elseif Y > -.5 && Z > 1.15 && Z <= 1.6
        n=3;
    elseif Y > -.5 && Z > 1.6
        n=4;
    elseif Y > -4.5 && Y <= -.5 && Z <= .3
        n=5;
    elseif Y > -4.5 && Y <= -.5 && Z > .3 && Z <=.98
        n=6;
    elseif Y > -4.5 && Y <= -.5 && Z > .98 && Z <= 1.38
        n=7;
    elseif Y > -4.5 && Y <= -.5 && Z > 1.38 && Z <= 2.04
        n=8;
    elseif Y > -4.5 && Y <= -.5 && Z > 2.04
        n=9;
    elseif Y > -7 && Y < -4.5 && Z <= .398
        n=10;
    elseif Y > -7 && Y < -4.5 && Z > .398 && Z <= .87
        n=11;
    elseif Y > -7 && Y < -4.5 && Z > .87 && Z <= 1.27
        n=12;
    elseif Y > -7 && Y < -4.5 && Z > 1.27 && Z <= 1.863
        n=13;
    elseif Y > -7 && Y < -4.5 && Z > 1.863
        n=14;
    end
    epsvec = [.1]; %Creates intial vectors, which will be use to store iteration history for p, eps, rho, etc
    p2vec = [];
    h2vec = []; 
    rho2vec = [];
    Depsvec = [];
    num = 0; %Iteration counter, outside loop
    while abs(Deps) > .001 %Will continue to loop through until epsilon stops changing
        p2 = pinf + (rhoinf*(u1^2)*(1-eps)); %Eq 14, Computes value for p2
        h2 = h1 + ((.5*u1^2)*(1-eps^2)); %Eq 15, Computes value for h2
        h2j = h2 - 1; %Arbitrary starting value for h2j, which will converge to h2 in the nested loop below
        rhoguess = rhoinf*(1/eps); %Initial guess for rho2, using current epsilon value
        while h2j < h2
            Y = log10(rhoguess/1.292); %Current value for Y
            X = log10(p2/101300); %Current value for X
            Z = X - Y; %Current value for Z
            gamma = G(n,1)+G(n,2)*Y+G(n,3)*Z+G(n,4)*Y*Z+((G(n,5)+G(n,6)*Y+G(n,7)*Z+G(n,8)*Y*Z)/(1+exp(G(n,9)*(X+G(n,10)*Y+G(n,11))))); %Computes Gamma tilda
            h2j = (p2/rhoguess)*(gamma/(gamma-1)); %Computes value for h2j
            rhoguess = rhoguess - .001; %Updates new value for rhoguess
        end
        eps = rhoinf/rhoguess; %Updates new value of epsilon
        num = num + 1; %Updates iteration counter
        epsvec = [epsvec;eps]; %Updates vectors of stored epsilon values
        Deps = epsvec(num+1,1)-epsvec(num,1); %Computes new value for change in epsilon, which is being driven to 0
        p2vec = [p2vec; p2]; %Updates remaining vectors
        h2vec = [h2vec; h2];
        rho2vec = [rho2vec; rhoguess];
        Depsvec = [Depsvec; Deps];
    end
    iteration = [0:(num-1)];
    iteration = transpose(iteration);
    epsvec(1)=[];
    iterationHist = table(iteration,epsvec,p2vec,h2vec,rho2vec,Depsvec);
    iterationHist.Properties.VariableNames = {'i','eps','p2','h2','rho2','Deps'};
    rho2 = rhoguess; %Assigns downstream density, rho2, to final updated value of rhoguess
    Y = log10(rho2/1.225); %Assigns X value for temperature region
    X = log10(p2/101340); %Assigns y value for temperature region
    Z = X - Y; %Assigns Z value for temperature region
    if Y > -.5 && Z > .48 && Z <= .90 %Locates the region for temperature d constants
        m=1;
    elseif Y > -.5 && Z > .90
        m=2;
    elseif Y > -4.5 && Y <= -.5 && Z > .48 && Z <= .9165
        m=3;
    elseif Y > -4.5 && Y <= -.5 && Z > .9165 && Z <=1.478
        m=4;
    elseif Y > -4.5 && Y <= -.5 && Z > 1.478 && Z <= 2.176
        m=5;
    elseif Y > -4.5 && Y <= -.5 && Z > 2.176
        m=6;
    elseif Y > -7 && Y < -4.5 && Z > .3 && Z <= 1.07
        m=7;
    elseif Y > -7 && Y < -4.5 && Z > 1.07 && Z <= 1.57
        m=8;
    elseif Y > -7 && Y < -4.5 && Z > 1.57 && Z <= 2.24
        m=9;
    elseif Y > -7 && Y < -4.5 && Z > 2.24
        m=10;
    end
    eqn23 = D(m,1)+D(m,2)*Y+D(m,3)*Z+D(m,4)*Y*Z+D(m,5)*Z^2+((D(m,6)+D(m,7)*Y+D(m,8)*Z+D(m,9)*Y*Z+D(m,10)*Z^2)/(1+exp(D(m,11)*(Z+D(m,12))))); %Eqn 23 for log10(T2/T0)
    TRatio = 10^eqn23; %Obtains ratio T2/T0
    T2 = TRatio*151.78; %Obtains downstream temperature T2
    u2 = u1*eps; %Obtains downstream velocity, u2
end

% 1a) Computing Lift and Drag Coeffs Using Modified Newtonian
r = rn/rc; %Ratio of Nose radius to Cone radius
CL = [];
CD = [];
Cpmax = (2-eps)*sind(90)^2; %Eqn 27, obtains Cpmax value based on compute epsilon value from shock calculations
for a = 0:1:10 %Using for loop to iterate through all 10 AOA's of Interest and Creating Vector Outputs for CL and CD
    CN = Cpmax*(1-(r^2*(cosd(dc)^2)))*((cosd(dc)^2)*sind(a)*cosd(a)); %Eq 36 - Normal Force Coefficient for Blunted Cone
    CA = Cpmax*((.5*(1-(sind(dc))^4)*r^2)+(((sind(dc)^2*(cosd(a)^2))+(.5*(cosd(dc)^2*(sind(a)^2))))*(1-(r^2*(cosd(dc)^2))))); %Eq 38 - Axial Force Coefficient for Blunted Cone
    CL1 = CN*cosd(a)-CA*sind(a); %Eq 28b - Lift Coefficient Using Modified Newtonian 
    CD1 = CN*sind(a)+CA*cosd(a); %Eq 28b - Drag Coefficient Using Modified Newtonian
    CL = [CL CL1]; %Adds newest value of lift to the output vector
    CD = [CD CD1];
end

% 2a) Finding heat flux at stagnation point using Fay & Riddell
% First, will compute the transport properties using soph. model
pnew = p2/101325; %Converts p2 from Pa to atm
B = 1.38E-23; %Boltzmann Constant
mf = [.78578 .205749 0 .0001626 .008304]; %Given mole fractions for 5 species model
mw = [28.014 31.998 14.007 15.999 30.006]; %Given molecular weights
s = [3.667 3.433 2.94 2.33 3.47]; %Sigma values (In Angstroms) from Bird
ek = [99.8 113 78 154 119]; %eps/k values from Bird
ep = ek.*B; %eps values for each species
ket = T2.*(1./ek); %(k/eps)*T
om = (1.16145./ket.^0.14874)+(0.52487./exp(0.7732.*ket))+(2.16178./exp(2.43787.*ket)); %omega values for mu and k
mui = 0.000026693.*sqrt(mw.*T2)./(s.^2.*om); %mu,i values in g/cm-s
xmu = mui.*mf; %Obtains Xmu,i values
phiarr = [];
for var = 1:5 %This double for loop iterates through species i and j to create 5x5 array of phi i,j values
    phivec = [];
    for var2 = 1:5
        phi = (1./sqrt(8)).*(1+(mw(var)./mw(var2))).^-0.5.*(1+(mui(var)./mui(var2)).^0.5.*(mw(var2)./mw(var)).^0.25).^2;
        phivec = [phivec phi];
    end
    phiarr = [phiarr;phivec];
end
sumvec = [];
for var3 = 1:5; %This single for loop iterates to obtain the sum(Xj*phi) vector
    summ = mf(1).*phiarr(var3,1) + mf(2).*phiarr(var3,2) + mf(3).*phiarr(var3,3) + mf(4).*phiarr(var3,4) + mf(5).*phiarr(var3,5);
    sumvec = [sumvec summ];
end
mu = xmu./sumvec; %Obtains mu value for each of the 5 species
mu = sum(mu); %mu value for mixture in g/cm-s
mum = (mu./1000).*100; %mu value for mixture in kg/m-s
xkwmk = 4.75.*mui.*(8.314./mw).*100.*mf; %Obtains X*k in W/m-K (Used for N2 O2 and NO)
kcal = .00019891.*(sqrt(T2./mw)./(s.^2.*om)).*mf.*418.634; %Obtains X*k in cal/cm-s-k, then converts to W/m-K (Used for N and O)
k = [xkwmk(1) xkwmk(2) kcal(3) kcal(4) xkwmk(5)]; %Obtains k values for each species, using the respective 2 equations above
kmix = sum(k);

% 2a) Using Fay-Riddell method for stagnation point heat flux
pw = p2; %Assuming pressure at the edge of the boundary layer is the same as pressure at the wall
Tw = 300; %Using cold wall assumption
rhow = pw/(287*Tw); %Obtains density at the wall
muw = (1.462E-6)*((Tw^1.5)/(Tw+110.4)); %Computes mu at the wall
kw = (1.993E-3)*((Tw^1.5)/(Tw+110.4)); %Computes k at the wall
he = h2/1000; %Converts edge (Downstream) enthalpy to kJ/kg
Cpwall = he/T2; %Computes specific heat at the wall
hw = Cpwall*Tw; %Computes the wall enthalpy
dudx = (1/rn)*((2*(p2-pinf))/rho2)^.5; %Computes due/dx in units 1/sec
Prandtl = (muw*Cpwall*1000)/kw; %Computes Prandtl Number at the wall
qstag = .76*(Prandtl)^-.6*(rho2*mum)^.4*(rhow*muw)^.1*sqrt(dudx)*(he-hw)*1000; %Eqn 89,Computes stagnation point heat flux, in W/m^2
qstag = qstag/10000; %Converts final answer to W/cm^2

% 2a) Using Reference Temperature Method for heat flux at specified x value
%First, entering the provided conditions downstream of the oblique shock
To = 1195.78;
po = 215619.6;
rhoo = .62567;
uo = 1464.59;
ho = 1279.1482;
ho0 = 2351.66;
cpo = 1.1754;
mwo = 28.85;
geq = (2.39683E-19)*To^5 - (3.0436E-15)*To^4 + (8.89216E-12)*To^3 + (2.77835E-8)*To^2 - (1.46939E-4)*To + 1.4517; %Using curve fit, eqn 107, for gamma eq
aeq = sqrt(geq*(po/rhoo)); %Obtains equivalent speed of sound
Me = uo/aeq; %Uses speed of sound and given velocity to obtain mach number at the boundary layer edge
Tref = To*(1+.032*Me^2+.58*((Tw/To)-1)); %Obtains reference temperature
muref = (1.462E-6)*((Tref^1.5)/(Tref+110.4)); %Using simple model for mu, assuming reference temperature is below 2000K
kref = (1.993E-3)*((Tref^1.5)/(Tref+112)); %Using simple model for k
cpref = (8.31935E-12)*Tref^3 - (8.6249E-8)*Tref^2 + (3.14166E-4)*Tref + .901439; %Using curve fit, eqn 109, for c eq
Pr = (muref*cpref*1000)/kref; %Computes Prandtl number
x1 = rn*(1-sind(dc)); %Computes initial x location, in meters
y1 = rn*cosd(dc); %Computes initial y location, in meters
y = y1 + tand(dc)*(x-x1); %Computes y location at specified x value
s = rn*deg2rad(90-dc) + sqrt((x-x1)^2+(y-y1)^2); %Obtains running lentgh, s
rhos = po/(287*Tref); %Obtains reference density
Res = (rhos*uo*s)/muref; %Computes reference Reynold's Number to see if flow is laminar or turbulent at point of interest
if Res >= 10^5 %Defines recovery factor and mf, depending on whether flow is turbulent or laminar
    mf = sqrt(2);
    rec = Pr^.33;
elseif Res < 10^5
    mf = sqrt(3);
    rec = sqrt(Pr);
end
Ch = (.332/sqrt(Res))*Pr^-.67*mf; %Obtains value for CH
Hr = ho + rec*(uo^2/2); %Computes recovery enthalpy
qref = (rhoo*uo*Ch)*(Hr-(cpo*Tw)); %Computes heat flux at point of interest, in W/m^2
qref = qref/10000; %Converts final answer to W/cm^2
  

























































end
