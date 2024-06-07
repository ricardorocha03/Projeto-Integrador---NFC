    

%% Asa    

AR=9.917;
S=0.898;
sections=4;
    dens_casca = 0.45; % [Kg/m^2]    usar 0.820 se for necessário ou 1.64
    n_oli = 9; 
    c_oli = 0.972; % [m]

    mass_casca = dens_casca * S * 2;  % duplicated to take into account the upper and lower parts of the wing

    envergadura = sqrt(AR * S); % [m]
    chord_raiz = (S/envergadura); %[m]
    
    n = (envergadura/2)*(n_oli)/(c_oli);  %calculation of the number of stringers/spars using a simple rule of three (based on data from OLISSIPO)
    espacamento = (envergadura/2/n);  %spacing between ribs [m]

    num_longarinas= 2; % this variable can either be assigned the value 1 or 2, with 2 being more realistic
    mass_longarinas = (30/34)* 10^(-1) * envergadura/2 * 2 * 1.4 *2^(min(num_longarinas,2)-1);
    
    if sections == 4
        mass_longarinas = mass_longarinas + (30/34)* 10^(-1)*2*(3*2*espacamento*10^-3) + (0.0115*8); %Assumindo configuração PEGASUS, com o segundo termo a massa das peças 3D 
    end
    
    %Obtaining and summing the rib areas - 1 panel
    area_nervuras = 0;    
    flecha= 3.72; %average value obtained through data for the geometry in the transportation box
    [x,y] = get_coordinates();
    for dist = 0:espacamento:envergadura/2      
        chord = chord_raiz - (tand(flecha)*dist);      %calculates the chord of each rib [mm]
    
        x_t=x.*chord;        %airfoil adapted to the chord length
        y_t=y.*chord;
        
        [k,A]=boundary(x_t,y_t);        %calculates the airfoil area
        
        area_nervuras=area_nervuras + A; % [m^2]
   
    end
    
    mass_nervuras = 2*(0.105/(0.305*0.330))*(area_nervuras);    %[kg] result for 2 panels 
    
    %data for calculating mass: plate measuring 305mm x330mm weighs 105g
    
    mass_resina = 0.030;
    
    mass_asa = (mass_nervuras + mass_longarinas + mass_casca + mass_resina)*1.15^(0.5*min(sections,4)-1);
    
    
    
function [x_airfoil, y_airfoil] = get_coordinates()
    % Define the matrix
    matrix = [
        1.000000	0.008475
        0.998993	0.009618
        0.995977	0.010918
        0.990964	0.012397
        0.983974	0.014074
        0.975036	0.015966
        0.964184	0.018088
        0.951463	0.020452
        0.936925	0.023066
        0.920627	0.025929
        0.902635	0.029031
        0.883022	0.032353
        0.861867	0.035866
        0.839255	0.039530
        0.815276	0.043298
        0.790028	0.047115
        0.763613	0.051021
        0.736136	0.054832
        0.707708	0.058494
        0.678443	0.061955
        0.648460	0.065165
        0.617879	0.068080
        0.586824	0.070660
        0.555419	0.072869
        0.523791	0.074678
        0.492067	0.076063
        0.460375	0.077008
        0.428843	0.077499
        0.397597	0.077527
        0.366763	0.077090
        0.336466	0.076192
        0.306827	0.074845
        0.277967	0.073063
        0.250000	0.070869
        0.223040	0.068289
        0.197195	0.065352
        0.172570	0.062092
        0.149263	0.058539
        0.127368	0.054726
        0.106973	0.050685
        0.088162	0.046446
        0.071008	0.042043
        0.055582	0.037506
        0.041946	0.032867
        0.030154	0.028160
        0.020254	0.023418
        0.012285	0.018676
        0.006281	0.013972
        0.002264	0.009342
        0.000252	0.004824
        0.000252	0.000458
        0.002264	-0.003714
        0.006281	-0.007650
        0.012285	-0.011309
        0.020254	-0.014651
        0.030154	-0.017642
        0.041946	-0.020251
        0.055582	-0.022457
        0.071008	-0.024243
        0.088162	-0.025606
        0.106973	-0.026552
        0.127368	-0.027097
        0.149263	-0.027265
        0.172570	-0.027090
        0.197195	-0.026609
        0.223040	-0.025864
        0.250000	-0.024896
        0.277967	-0.023747
        0.306827	-0.022454
        0.336466	-0.021053
        0.366763	-0.019576
        0.397597	-0.018050
        0.428843	-0.016498
        0.460375	-0.014941
        0.492067	-0.013394
        0.523791	-0.011871
        0.555419	-0.010383
        0.586824	-0.008939
        0.617879	-0.007545
        0.648460	-0.006210
        0.678443	-0.004940
        0.707708	-0.003742
        0.736136	-0.002624
        0.763613	-0.001593
        0.790028	-0.000657
        0.815276	0.000178
        0.839255	0.000911
        0.861867	0.001542
        0.883022	0.002073
        0.902635	0.002505
        0.920627	0.002844
        0.936925	0.003095
        0.951463	0.003269
        0.964184	0.003375
        0.975036	0.003423
        0.983974	0.003423
        0.990964	0.003384
        0.995977	0.003316
        0.998993	0.003225
        1.000000	0.003118

    ];

    % Extract the first and second columns
    x_airfoil = matrix(:, 1);
    y_airfoil = matrix(:, 2);
end
