clc;
clear all;
% diretoria
dir = '.\CaudasVersaoTesteMassas2';
mkdir(dir);
for Vht = 0.7
    for Vvt = 0.05
        for AR_ht = 3.5
            for AR_vt = 3
                
wing_span = 3.0;%em metros
S = 0.897;
S1=0.325;
span_1=wing_span*0.3
span_2=wing_span*0.7
S2=S-S1
AR=wing_span^2/S
Cmed = S/wing_span; % em metros

AR1=span_1^2/S1
AR2=span_2^2/S2
c_m = (4 * S2) / (pi * span_2)
c_r=2*S1/span_1-c_m
% TR2=(2*span_2)/(c_m*AR2)-1;
% c_t=c_m*TR2;

%Y_asa = @(t)((t/n_sections_asa)*(wing_span/2));            %função que define a posição de cada secção
%Y_asa = [0 0.308 0.399 0.639 0.709 0.774 0.835 0.892 0.945 0.976 0.993 1.008 1.022]*(wing_span/2)/1.022;
n_sections_asa_1 = 2; 
n_sections_asa_2 = 20;
n_sections_asa=n_sections_asa_1+n_sections_asa_2;
%definir o número de secções pretendidas
Y_max = wing_span/2;
Y_asa_1 = @(t)(((t-1)/(n_sections_asa_1-1))*(span_1/2)); 
growth_rate = 1.045;
Y_asa_2 = @(t)((((t-2)*growth_rate^(n_sections_asa_2-t+2))/n_sections_asa_2)*(span_2/2)+span_1/2);
%função que define a posição de cada secção
%n_sections_asa = length(Y_asa); 

corda_asa_1 = @(y)(c_r+2*y*(c_m-c_r)/span_1); %função que define a corda em função da distância à raíz
corda_asa_2 = @(y)(c_m * sqrt(1 - ((2 * (y-0.45)) / span_2)^2)+0.001);                
                
                
R1=0.30;
R2=0.02;
lt = sqrt(2 * S* (Vht * Cmed + Vvt * wing_span) / (pi * (R1 + R2)));
Sht = Vht * S * Cmed / lt;
Svt = Vvt * S * wing_span / lt;
bht = sqrt(AR_ht * Sht);
bvt = sqrt(AR_vt * Svt);
cht = Sht/bht;
cvt = Svt/bvt;

Vol_htail = Sht*0.08*cht;  
Vol_vtail = Svt*0.08*cvt;
m_htail= (Vol_htail*38.45)+(Sht*4.5*10^-2)+(bht*0.117809*10^-1)+0.4690846*2.04*Sht;
m_vtail = (Vol_vtail*38.45)+(Svt*4.5*10^-2)+(bvt*0.117809*10^-1)+0.4690846*2.04*Svt;  
m_pontual=11-2.07-m_htail-m_vtail;
             
nome = sprintf('V1_cauda_%0.4d_%0.4d_%0.4d_%0.4d', int64(Vht*1000),int64(Vvt*1000),int64(AR_ht*100),int64(AR_vt*100));                
        filename = strcat(nome, '.xml');               
        fid = fopen(fullfile(dir, filename),'wt');  
if fid<0
   fprintf('erro ao abrir o ficheiro\n');
   return;
end
format long;
fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid,'<!DOCTYPE explane>\n');
fprintf(fid,'<explane version="1.0">\n');
fprintf(fid,'<Units>\n');
fprintf(fid,'<length_unit_to_meter>1</length_unit_to_meter>\n');
fprintf(fid,'<mass_unit_to_kg>1</mass_unit_to_kg>\n');
fprintf(fid,'</Units>\n');

%-------------------------------------PLANE-------------------------------------%
%---------------Definições do avião---------------%
plane_name = nome;
plane_description = 'miau';
body = 'false';       %definir se há corpo ou não
%---------------Definições do avião---------------%

fprintf(fid,'<Plane>\n');
fprintf(fid,'<Name> %s </Name>\n',plane_name);
fprintf(fid,'<Description>\n');
fprintf(fid,'%s\n',plane_description);
fprintf(fid,'</Description>\n');
%------------Definir inércia------------%
fprintf(fid,'<Inertia>\n');
%------------Definir massas------------%
fprintf(fid,'<Point_Mass>\n');
fprintf(fid,'<Tag></Tag>\n');
fprintf(fid,strcat('<Mass> ',num2str(m_pontual),' </Mass>\n')); 
% fprintf(fid,'<Mass> ' , num2str(m_tail),' </Mass>\n');
fprintf(fid,'<coordinates>  0.05,   0, -0.07 </coordinates>\n');
fprintf(fid,'</Point_Mass>\n');
%----------------Massas----------------%
fprintf(fid,'</Inertia>\n');
%----------------Inércia----------------%
fprintf(fid,'<has_body> %s </has_body>\n',body);         


%--------------------------------------ASA--------------------------------------%
%---------------Definições do asa---------------%
name_wing = 'Asa PI';    %Nomear a asa
wing_description = '';      
tilt_angle_asa = 0;             %ângulo em graus  
pos_asa = [0 0 0];          %definir posição
%---------------Definições do asa---------------%

fprintf(fid,'<wing>\n');
fprintf(fid,'<Name>%s</Name>\n',name_wing);
fprintf(fid,'<Type>MAINWING</Type>\n');
%-----random definições de cor-----------%
fprintf(fid,'<Color>\n');
fprintf(fid,'<red>208</red>\n');                
fprintf(fid,'<green>145</green>\n');                
fprintf(fid,'<blue>140</blue>\n');                
fprintf(fid,'<alpha>255</alpha>\n');                
fprintf(fid,'</Color>\n'); 
%-----alterar valores numéricos para alterar a cor---%
fprintf(fid,'<Description>%s</Description>\n',wing_description); 
fprintf(fid,'<Position> %0.3f, %0.3f, %0.3f</Position>\n',pos_asa(1),pos_asa(2),pos_asa(3)); 
fprintf(fid,' <Tilt_angle> %0.3f </Tilt_angle>\n',tilt_angle_asa);                      
fprintf(fid,'<Symetric>true</Symetric>\n');          %convém né
fprintf(fid,'<isFin>false</isFin>\n');               
fprintf(fid,'<isDoubleFin>false</isDoubleFin>\n'); 
fprintf(fid,'<isSymFin>false</isSymFin>\n'); 
%-----Definir inércia da asa-----%
fprintf(fid,'<Inertia>\n'); 
fprintf(fid,'<Volume_Mass>  2.07</Volume_Mass>\n'); 
fprintf(fid,'</Inertia>\n'); 
%--------------Inércia-----------%

%------------Início das sections------------%
%------------------------------------Definições das sections------------------------------------------%
                                      %definir o número de secções pretendidas
                                    
% wing_span = 3.0;%em metros
% S = 0.897;
% S1=0.325;
% span_1=wing_span*0.3
% span_2=wing_span*0.7
% S2=S-S1
% AR=wing_span^2/S
% Cmed = S/wing_span; % em metros
% 
% AR1=span_1^2/S1
% AR2=span_2^2/S2
% c_m = (4 * S2) / (pi * span_2)
% c_r=2*S1/span_1-c_m
% % TR2=(2*span_2)/(c_m*AR2)-1;
% % c_t=c_m*TR2;
% 
% %Y_asa = @(t)((t/n_sections_asa)*(wing_span/2));            %função que define a posição de cada secção
% %Y_asa = [0 0.308 0.399 0.639 0.709 0.774 0.835 0.892 0.945 0.976 0.993 1.008 1.022]*(wing_span/2)/1.022;
% n_sections_asa_1 = 2; 
% n_sections_asa_2 = 20;
% n_sections_asa=n_sections_asa_1+n_sections_asa_2;
% %definir o número de secções pretendidas
% Y_max = wing_span/2;
% Y_asa_1 = @(t)(((t-1)/(n_sections_asa_1-1))*(span_1/2)); 
% growth_rate = 1.045;
% Y_asa_2 = @(t)((((t-2)*growth_rate^(n_sections_asa_2-t+2))/n_sections_asa_2)*(span_2/2)+span_1/2);
% %função que define a posição de cada secção
% %n_sections_asa = length(Y_asa); 
% 
% corda_asa_1 = @(y)(c_r+2*y*(c_m-c_r)/span_1); %função que define a corda em função da distância à raíz
% corda_asa_2 = @(y)(c_m * sqrt(1 - ((2 * (y-0.45)) / span_2)^2)+0.001); %função que define a corda em função da distância à raíz

    %corda_asa_2 = @(y)(c_m+2*(y-wing_span*0.3)*(c_t-c_m)/span_2); %função que define a corda em função da distância à raíz

    
%theta = acos(Y_asa);
%C_l = 0.39*(1-Y_asa.^2).^(3/2);
%twist_asa = @(y)(-((3*(C_l(1)*corda_asa(Y_asa(1))))/(8*wing_span))*cos(2*acos(y))*(180/pi)*0.5);
%twist_asa = @(y)(-sin((pi*y)/wing_span));                  %função que define a twist em função da distância à raíz
twist_asa = @(y)(0);
diedro_1_asa = 2;
%diedro_2 = ; caso haja mais que um diedro ao longo da asa;
%diedro_change = ; Y a partir do qual o diedro muda de valor;
x_paineis_asa = 10;                                        %número de paineis em x
x_distribuicao_asa = 'COSINE';                             %escolher entre 'COSINE', 'UNIFORM', 'SINE' ou '-SINE'
y_paineis_asa = 8;                                         %número de paineis em x
y_distribuicao_asa = 'COSINE';                             %escolher entre 'COSINE', 'UNIFORM', 'SINE' ou '-SINE'
%editar isto com base nas espessuras
vetor_espessuras=[0.3003,0.3003,0.3003,0.3003,0.3003,0.3003,0.3013,0.3093,0.3163,0.3163,0.3333];
 
% 
% % Foil_Mixer
% %Foil Mixer Proj Integrador
% format long
% 
% % Input
% n_section = 11;
% 
% % Distribution ratio (bigger = more tip influence)
% R = 0.3;
% % T parameter for bezier
% distribuition_density = 0.00001; 
% 
% % Foils (type XFLR5 .dat)
% fsroot = 'Root_final.dat';
% fstip = 'Tip_final.dat';
% 
% % Output directory
% dir = '.\PERFIS_PI_R0.3';
% 
% % \Input
% 
% % Read input foils
% froot = fopen(fsroot,'r');
% ftip = fopen(fstip,'r');
% % Foil name
% sroot = string(textscan(froot, '%s\n', 1));
% stip = string(textscan(ftip, '%s\n', 1));
% % Foil (x, y)
% croot = textscan(froot, '%f %f');
% ctip = textscan(ftip, '%f %f');
% ogFoils = {cell2mat(croot); cell2mat(ctip)};
% % Cleanup
% fclose(froot);
% fclose(ftip);
% 
% % foils: extraRoot{1,1}, intraRoot{1,2};
% %         extraTip{2,1},  intraTip{2,2}
% foils = {0,0;0,0};
% for i = 1 : 2
%     for j = 1 : length(ogFoils{i})
%         if ogFoils{i}(j, 2) < 0
%             foils{i, 1} = ogFoils{i}(1:j-1, :);
%             foils{i, 2} = ogFoils{i}(j:end, :);
%             break
%         end
%     end
% end
% 
% % Root Ratio along sections
% 
% % Quadratic Bezier curve
% B = @(t, a, b, c) (1-t).^2.*a + 2*(1-t).*t.*b + t.^2.*c;
% t = 0 : distribuition_density : 1;
% 
% % Define [starting control ending] (point)
% curve1 = {  [0; 0]   0 [1 - R; R]};
% curve2 = {[1 - R; R] 0   [1; 1]  };
% % Control dependent on ratio
% if R > 0.5
%     curve1{2} = [   0   ; 2*R-1];
%     curve2{2} = [2*(1-R);   1  ];
% else
%     curve1{2} = [1-2*R;  0 ];
%     curve2{2} = [  1  ; 2*R];
% end
% 
% % Calculate points on the curve
% bezier = transpose([B(t, curve1{1}, curve1{2}, curve1{3}) B(t, curve2{1}, curve2{2}, curve2{3})]);
% 
% % Search for points to fill in sectionRatio array
% sectionRatio = zeros(1, n_section);
% for i = 1 : n_section
%     x = (i - 1)/(n_section - 1);
%     for j = 1 : length(bezier) - 1 
%         y = linInterp(bezier(j, :), x, bezier(j + 1, :));       
%         if y ~= 0 
%             sectionRatio(i) = y; 
%         end
%     end
% end
% 
% % foil calc & write
% mkdir(dir);
% for i = 1 : n_section    
%     filename = sprintf('Perfil_PI%03i.dat', i);
%     fout = fopen(fullfile(dir, filename), 'wt');
%     fprintf(fout, 'Perfil_PI%03i\n', i);
%     for j = 1 : 2                
%         half = newHalfFoil(foils{1,j}, foils{2, j}, sectionRatio(i));
%         for k = 1 : length(half)
%             fprintf(fout, '%f %f\n', half(k, 1), half(k, 2));
%         end       
%     end
%     fclose(fout);
% end
% 
% % Given 2 matrices(half foil), 2 columns (x, y), and a (tip) ratio, create new (half) foil
% function hfoil = newHalfFoil(root, tip, ratio)
%     hfoil = root;
%     for i = 1 : length(root)
%         x = root(i, 1);        
%         for j = 1 : length(tip) - 1
%             y_tip = linInterp(tip(j, :), x, tip(j + 1, :));
%             if y_tip ~= 0
%                 y = (1 - ratio) * root(i, 2) + ratio * y_tip;
%                 hfoil(i, 1) = x;
%                 hfoil(i, 2) = y;
%                 break
%             end            
%         end
%     end
% end
% 
% % Linear Interpolation between 2 points, return 0 for invalid x
% function y = linInterp(p1, x, p2)
%     y = 0;
%     x1 = p1(1); y1 = p1(2);
%     x2 = p2(1); y2 = p2(2);
%     if (x1 <= x && x <= x2) || (x2 <= x && x <= x1)
%         m = (y2 - y1) / (x2 - x1);
%         y = y1 + m * (x - x1);
%     end
% end
% 




    %------------------------------------Definições das sections------------------------------------------%

    
    




fprintf(fid,'<Sections>\n');

for i = 1 :  21 
        fprintf(fid,'<Section>\n');
        if i==1
            fprintf(fid,'<y_position> %.3f </y_position>\n',Y_asa_1(i));
            fprintf(fid,'<Chord> %.3f </Chord>\n',corda_asa_1(Y_asa_1(i)));
            fprintf(fid,'<xOffset> %.3f </xOffset>\n',0);
        end
%         if i<=(n_sections_asa*0.6+1)
%             %fprintf(fid1,'<xOffset> %.3f </xOffset>\n',((0.2823*c_r) - (0.2432*c_t)))
%             fprintf(fid,'<Chord> %.3f </Chord>\n',corda_asa_1(Y_asa(i)));
%             fprintf(fid,'<xOffset> %.3f </xOffset>\n',(vetor_espessuras(1)*corda_asa_1(Y_asa(1))-vetor_espessuras(i)*corda_asa_1(Y_asa(i))));
%         end
%         if i>(n_sections_asa*0.6+1)
%             fprintf(fid,'<Chord> %.3f </Chord>\n',corda_asa_2(Y_asa(i)));
%             fprintf(fid,'<xOffset> %.3f </xOffset>\n',(vetor_espessuras(1)*corda_asa_1(Y_asa(1))-vetor_espessuras(i)*corda_asa_2(Y_asa(i))));
%         end
        if i<=2
            fprintf(fid,'<y_position> %.3f </y_position>\n',Y_asa_1(i));
            %fprintf(fid1,'<xOffset> %.3f </xOffset>\n',((0.2823*c_r) - (0.2432*c_t)))
            fprintf(fid,'<Chord> %.3f </Chord>\n',corda_asa_1(Y_asa_1(i)));
            fprintf(fid,'<xOffset> %.3f </xOffset>\n',((corda_asa_1(Y_asa_1(1))-corda_asa_1(Y_asa_1(i))))/2);
        end
        if i>2
            fprintf(fid,'<y_position> %.3f </y_position>\n',Y_asa_2(i));
            fprintf(fid,'<Chord> %.3f </Chord>\n',corda_asa_2(Y_asa_2(i)));
            fprintf(fid,'<xOffset> %.3f </xOffset>\n',((corda_asa_1(Y_asa_1(1))-corda_asa_2(Y_asa_2(i))))/2);
        end
      
    %-----------Caso haja dois diedros diferentes-----------%    
        sprintf('perfil_PI_r0.3%03i.dat', i);
            fprintf(fid,'<Dihedral> %.3f </Dihedral>\n',diedro_1_asa);
        foil_asa = sprintf('perfil_PI_r0.3%03i.dat', i);
    %-----------Diedros-----------%     
        fprintf(fid,'<Twist> %.3f </Twist>\n',twist_asa(Y_asa_1(i)));
        fprintf(fid,'<x_number_of_panels> %d </x_number_of_panels>\n',x_paineis_asa);
        fprintf(fid,'<x_panel_distribution>%s</x_panel_distribution>\n',x_distribuicao_asa);
        fprintf(fid,'<y_number_of_panels> %d </y_number_of_panels>\n',y_paineis_asa);
        fprintf(fid,'<y_panel_distribution>%s</y_panel_distribution>\n',y_distribuicao_asa);
        if i==1
        foil_asa = sprintf('Perfil_PI_r0.3001');    
        fprintf(fid,'<Left_Side_FoilName>%s</Left_Side_FoilName>\n',foil_asa);
        fprintf(fid,'<Right_Side_FoilName>%s</Right_Side_FoilName>\n',foil_asa);
        end
        if i==2
        foil_asa = sprintf('Perfil_PI_r0.3004');    
        fprintf(fid,'<Left_Side_FoilName>%s</Left_Side_FoilName>\n',foil_asa);
        fprintf(fid,'<Right_Side_FoilName>%s</Right_Side_FoilName>\n',foil_asa);
        end
        if i==3
        foil_asa = sprintf('Perfil_PI_r0.3005');    
        fprintf(fid,'<Left_Side_FoilName>%s</Left_Side_FoilName>\n',foil_asa);
        fprintf(fid,'<Right_Side_FoilName>%s</Right_Side_FoilName>\n',foil_asa);
        end
        if i>3 && i<=5
        foil_asa = sprintf('Perfil_PI_r0.3006');    
        fprintf(fid,'<Left_Side_FoilName>%s</Left_Side_FoilName>\n',foil_asa);
        fprintf(fid,'<Right_Side_FoilName>%s</Right_Side_FoilName>\n',foil_asa);
        end
        if i==6
        foil_asa = sprintf('Perfil_PI_r0.3007');    
        fprintf(fid,'<Left_Side_FoilName>%s</Left_Side_FoilName>\n',foil_asa);
        fprintf(fid,'<Right_Side_FoilName>%s</Right_Side_FoilName>\n',foil_asa);
        end
        if i>6 && i<=9
        foil_asa = sprintf('Perfil_PI_r0.3008');    
        fprintf(fid,'<Left_Side_FoilName>%s</Left_Side_FoilName>\n',foil_asa);
        fprintf(fid,'<Right_Side_FoilName>%s</Right_Side_FoilName>\n',foil_asa);
        end
        if i>9 && i<=12
        foil_asa = sprintf('Perfil_PI_r0.3009');    
        fprintf(fid,'<Left_Side_FoilName>%s</Left_Side_FoilName>\n',foil_asa);
        fprintf(fid,'<Right_Side_FoilName>%s</Right_Side_FoilName>\n',foil_asa);
        end
        if i>12 && i<=17
        foil_asa = sprintf('Perfil_PI_r0.3010');    
        fprintf(fid,'<Left_Side_FoilName>%s</Left_Side_FoilName>\n',foil_asa);
        fprintf(fid,'<Right_Side_FoilName>%s</Right_Side_FoilName>\n',foil_asa);
        end
        if i>17 
        foil_asa = sprintf('Perfil_PI_r0.3011');    
        fprintf(fid,'<Left_Side_FoilName>%s</Left_Side_FoilName>\n',foil_asa);
        fprintf(fid,'<Right_Side_FoilName>%s</Right_Side_FoilName>\n',foil_asa);
        end
        fprintf(fid,'</Section>\n');
end

   
    fprintf(fid,'</Sections>\n');
    %------------Fim das sections------------%
    fprintf(fid,'</wing>\n');
    


%{
%------------------------------------Definições das sections------------------------------------------%
fprintf(fid,'<Sections>\n'); 

for i = 0 : n_sections_asa
    fprintf(fid,'<Section>\n');
    fprintf(fid,'<y_position> %.3f </y_position>\n',Y_asa(i));
    fprintf(fid,'<Chord> %.3f </Chord>\n',corda_asa(Y_asa(i)));
    fprintf(fid,'<xOffset> %.3f </xOffset>\n',(corda_asa(0)-corda_asa(Y_asa(i)))/2);
%-----------Caso haja dois diedros diferentes-----------%    
    %if Y(i) < diedro_change
        fprintf(fid,'<Dihedral> %.3f </Dihedral>\n',diedro_1_asa);
    %else 
    %    fprintf(fid,'<Dihedral> %.3f </Dihedral>\n',diedro_2);
    %end
%-----------Diedros-----------%     
    fprintf(fid,'<Twist> %.3f </Twist>\n',twist_asa(Y_asa(i)));
    fprintf(fid,'<x_number_of_panels> %d </x_number_of_panels>\n',x_paineis_asa);
    fprintf(fid,'<x_panel_distribution>%s</x_panel_distribution>\n',x_distribuicao_asa);
    fprintf(fid,'<y_number_of_panels> %d </y_number_of_panels>\n',y_paineis_asa);
    fprintf(fid,'<y_panel_distribution>%s</y_panel_distribution>\n',y_distribuicao_asa);
    fprintf(fid,'<Left_Side_FoilName>%s</Left_Side_FoilName>\n',foil_asa);
    fprintf(fid,'<Right_Side_FoilName>%s</Right_Side_FoilName>\n',foil_asa);
    fprintf(fid,'</Section>\n');
end

fprintf(fid,'</Sections>\n');
%------------Fim das sections------------%
fprintf(fid,'</wing>\n');
%--------------------------------------Fim da ASA--------------------------------------%
%}

%--------------------------------------ELEVATOR--------------------------------------%
%---------------Definições do elevator---------------%
% R1=0.30;
% R2=0.02;
% lt = sqrt(2 * S* (Vht * Cmed + Vvt * wing_span) / (pi * (R1 + R2)));
% Sht = Vht * S * Cmed / lt;
% Svt = Vvt * S * wing_span / lt;
% bht = sqrt(AR_ht * Sht);
% bvt = sqrt(AR_vt * Svt);
% cht = Sht/bht;
% cvt = Svt/bvt;


name_elevator = 'Elevator';    %Nomear o elevator
elevator_description = '';      
tilt_angle_elev = 0;          %                    <--------ângulo em graus  
pos_elev = [lt+0.110-cht*0.25 0 -0.07];          %      <--------definir posição
%---------------Definições do elevator---------------%

fprintf(fid,'<wing>\n');
fprintf(fid,'<Name>%s</Name>\n',name_elevator);
fprintf(fid,'<Type>ELEVATOR</Type>\n');
%-----random definições de cor-----------%
fprintf(fid,'<Color>\n');
fprintf(fid,'<red>138</red>\n');                
fprintf(fid,'<green>216</green>\n');                
fprintf(fid,'<blue>140</blue>\n');                
fprintf(fid,'<alpha>255</alpha>\n');                
fprintf(fid,'</Color>\n'); 
%-----alterar valores numéricos para alterar a cor---%
fprintf(fid,'<Description>%s</Description>\n',elevator_description); 
fprintf(fid,'<Position> %0.3f, %0.3f, %0.3f</Position>\n',pos_elev(1),pos_elev(2),pos_elev(3)); 
fprintf(fid,' <Tilt_angle> %0.3f </Tilt_angle>\n',tilt_angle_elev);                      
fprintf(fid,'<Symetric>true</Symetric>\n');          %convém né
fprintf(fid,'<isFin>false</isFin>\n');               
fprintf(fid,'<isDoubleFin>false</isDoubleFin>\n'); 
fprintf(fid,'<isSymFin>false</isSymFin>\n'); 
%-----Definir inércia da asa-----%
fprintf(fid,'<Inertia>\n'); 
fprintf(fid,strcat('<Volume_Mass> ',num2str(m_htail),' </Volume_Mass>\n')); 
fprintf(fid,'</Inertia>\n'); 
%--------------Inércia-----------%

%------------Início das sections------------%
%------------------------------------Definições das sections------------------------------------------%
n_sections_elev = 2;      %definir o número de secções pretendidas 
%         elev_span = sqrt(ARvt(b)*V_area(c));                                           %em metros
         Y_elev = [0 bht/2];          %função que define a posição de cada secção
%         
        corda_elev = @(y)(cht); %função que define a corda em função da distância à raíz
        twist_elev = @(y)(0);                  %função que define a twist em função da distância à raíz
        diedro_1_elev = 0;
        offset_elev =0;
        %diedro_2 = ; caso haja mais que um diedro ao longo da asa;
        %diedro_change = ; Y a partir do qual o diedro muda de valor;
        x_paineis_elev = 10;                                        %número de paineis em x
        x_distribuicao_elev = 'COSINE';                               %escolher entre 'COSINE', 'UNIFORM', 'SINE' ou '-SINE'
        y_paineis_elev = 7;                                         %número de paineis em x
        y_distribuicao_elev = 'COSINE';                            %escolher entre 'COSINE', 'UNIFORM', 'SINE' ou '-SINE'
        foil_elev = 'NACA 0008';
%------------------------------------Definições das sections------------------------------------------%
        fprintf(fid,'<Sections>\n'); 

        for i = 1 : n_sections_elev 
            fprintf(fid,'<Section>\n');
            fprintf(fid,'<y_position> %.3f </y_position>\n',Y_elev(i));
            fprintf(fid,'<Chord> %.3f </Chord>\n',corda_elev(Y_elev(i)));
            fprintf(fid,'<xOffset> %.3f </xOffset>\n',offset_elev);
%-----------Caso haja dois diedros diferentes-----------%    
            %if Y(i) < diedro_change
             fprintf(fid,'<Dihedral> %.3f </Dihedral>\n',diedro_1_elev);
            %else 
            %fprintf(fid,'<Dihedral> %.3f </Dihedral>\n',diedro_2);
            %end
%-----------Diedros-----------%     
            fprintf(fid,'<Twist> %.3f </Twist>\n',twist_elev(Y_elev(i)));
            fprintf(fid,'<x_number_of_panels> %d </x_number_of_panels>\n',x_paineis_elev);
            fprintf(fid,'<x_panel_distribution> %s </x_panel_distribution>\n',x_distribuicao_elev);
            fprintf(fid,'<y_number_of_panels> %d </y_number_of_panels>\n',y_paineis_elev);
            fprintf(fid,'<y_panel_distribution> %s </y_panel_distribution>\n',y_distribuicao_elev);
            fprintf(fid,'<Left_Side_FoilName>%s</Left_Side_FoilName>\n',foil_elev);
            fprintf(fid,'<Right_Side_FoilName>%s</Right_Side_FoilName>\n',foil_elev);
            fprintf(fid,'</Section>\n');
        end
        fprintf(fid,'</Sections>\n');
%------------Fim das sections------------%
        fprintf(fid,'</wing>\n');
        
   
%--------------------------------------Fim de ELEVATOR--------------------------------------%


%--------------------------------------FIN--------------------------------------%
%---------------Definições do fin---------------%
name_fin = 'Fin';    %Nomear o fin
fin_description = '';      
tilt_angle_fin = 0;               %ângulo em graus  
pos_fin = [lt+0.110-cvt*0.25 0 -0.07];          %definir posição
%---------------Definições do fin---------------%

fprintf(fid,'<wing>\n');
fprintf(fid,'<Name>%s</Name>\n',name_fin);
fprintf(fid,'<Type>FIN</Type>\n');
%-----random definições de cor-----------%
fprintf(fid,'<Color>\n');
fprintf(fid,'<red>233</red>\n');                
fprintf(fid,'<green>201</green>\n');                
fprintf(fid,'<blue>110</blue>\n');                
fprintf(fid,'<alpha>255</alpha>\n');                
fprintf(fid,'</Color>\n'); 
%-----alterar valores numéricos para alterar a cor---%
fprintf(fid,'<Description>%s</Description>\n',fin_description); 
fprintf(fid,'<Position> %0.3f, %0.3f, %0.3f</Position>\n',pos_fin(1),pos_fin(2),pos_fin(3)); 
fprintf(fid,' <Tilt_angle> %0.3f </Tilt_angle>\n',tilt_angle_fin);                      
fprintf(fid,'<Symetric>true</Symetric>\n');          %convém né
fprintf(fid,'<isFin>true</isFin>\n');                
fprintf(fid,'<isDoubleFin>false</isDoubleFin>\n'); 
fprintf(fid,'<isSymFin>false</isSymFin>\n'); 
%-----Definir inércia da asa-----%
fprintf(fid,'<Inertia>\n'); 
fprintf(fid,strcat('<Volume_Mass> ',num2str(m_vtail),' </Volume_Mass>\n')); 
fprintf(fid,'</Inertia>\n'); 
%--------------Inércia-----------%

%------------Início das sections------------%
%------------------------------------Definições das sections------------------------------------------%
n_sections_fin = 2;                                       %definir o número de secções pretendidas
% fin_span = 0.44;                                           %em metros
Y_fin = [0 bvt/2];             %função que define a posição de cada secção
corda_fin = @(y)(cvt); %função que define a corda em função da distância à raíz
twist_fin = @(y)(0); 
offset_fin=0;            %função que define a twist em função da distância à raíz
diedro_1_fin = 0;
%diedro_2 = ; caso haja mais que um diedro ao longo da asa;
%diedro_change = ; Y a partir do qual o diedro muda de valor;
x_paineis_fin = 10;                                        %número de paineis em x
x_distribuicao_fin = 'COSINE';                            %escolher entre 'COSINE', 'UNIFORM', 'SINE' ou '-SINE'
y_paineis_fin = 7;                                         %número de paineis em x
y_distribuicao_fin = 'COSINE';                             %escolher entre 'COSINE', 'UNIFORM', 'SINE' ou '-SINE'
foil_fin = foil_elev;
%------------------------------------Definições das sections------------------------------------------%
fprintf(fid,'<Sections>\n'); 

for i = 1 : n_sections_fin
    fprintf(fid,'<Section>\n');
    fprintf(fid,'<y_position> %.3f </y_position>\n',Y_fin(i));
    fprintf(fid,'<Chord> %.3f </Chord>\n',corda_fin(Y_fin(i)));
    fprintf(fid,'<xOffset> %.3f </xOffset>\n',offset_fin);
%-----------Caso haja dois diedros diferentes-----------%    
    %if Y(i) < diedro_change
        fprintf(fid,'<Dihedral> %.3f </Dihedral>\n',diedro_1_fin);
    %else 
    %    fprintf(fid,'<Dihedral> %.3f </Dihedral>\n',diedro_2);
    %end
%-----------Diedros-----------%     
    fprintf(fid,'<Twist> %.3f </Twist>\n',twist_fin(Y_fin(i)));
    fprintf(fid,'<x_number_of_panels> %d </x_number_of_panels>\n',x_paineis_fin);
    fprintf(fid,'<x_panel_distribution> %s </x_panel_distribution>\n',x_distribuicao_fin);
    fprintf(fid,'<y_number_of_panels> %d </y_number_of_panels>\n',y_paineis_fin);
    fprintf(fid,'<y_panel_distribution> %s </y_panel_distribution>\n',y_distribuicao_fin);
    fprintf(fid,'<Left_Side_FoilName>%s</Left_Side_FoilName>\n',foil_elev);
    fprintf(fid,'<Right_Side_FoilName>%s</Right_Side_FoilName>\n',foil_elev);
    fprintf(fid,'</Section>\n');
end

fprintf(fid,'</Sections>\n');
%------------Fim das sections------------%
fprintf(fid,'</wing>\n');
%--------------------------------------Fim de FIN--------------------------------------%


fprintf(fid,'</Plane>\n');
fprintf(fid,'</explane>\n');
%--------------------------------------Fim do PLANE--------------------------------------%
fclose(fid);
            end
        end
    end
end