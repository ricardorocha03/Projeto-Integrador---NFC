mtom=13.82 ; %kg
mwing=2.07 ;
mtail=0.42+0.285 ;
mbooms=0.21*3 ;
mmotors=1.28 + 0.8 ;
meletro=2.786 + 0.896 ;
mtrem=0.26 ;

%available mass for fuselage+payload
mavai=mtom-mwing-mtail-mbooms-mmotors-meletro-mtrem;


bestp=0;
bestd=0;
bestVpay=0.01;
bestMpay=0.01;

dmin = 0.1;
dmax = 0.35;

for d = dmin:0.005:dmax

    l = d/0.3;
    
    %volume and surface area calculation
    Vfus = pi*(d/2)^3/(3*tan(20)) + pi*(d/2)^2*(l-d/2*(1+1/tan(20))) + 2/3*pi*(d/2)^3;
    Sfus = pi*(d/2)^2/sin(20) +pi*d*(l-d/2*(1+1/tan(20))) + 2*pi*(d/2)^2;

    %fuselage mass calculation
    Mfus = 0.82*Sfus + 700*0.003*pi*((d/2)^2-(d/2-0.1*d/2)^2);

    if Mfus >= mavai
        break
    end
    
    % mass and volume of payload
    Mpay = mavai - Mfus;
    Vpay = Vfus*0.7;

    %comparative points calculation
    points = 18*tanh(Vpay/bestVpay-1) + 10*tanh(Mpay/bestMpay-1);

    if points > 0
        bestd = d;
        bestVpay=Vpay;
        bestMpay=Mpay;
    end

end

disp(bestd)
disp(bestVpay)
disp(bestMpay)



