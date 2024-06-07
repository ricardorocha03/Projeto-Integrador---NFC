fid = fopen('elipse_PI_20_rate1.05.xml','wt');     %alterar nome do ficheiro aqui
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
perfil = 'NACA 0012';
plane_name = 'elipse_PI_20_rate1.05';
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
fprintf(fid,'<Mass>  6.200 </Mass>\n');
fprintf(fid,'<coordinates>  0.085,   0, -0.2</coordinates>\n');
fprintf(fid,'</Point_Mass>\n');
%----------------Massas----------------%
fprintf(fid,'</Inertia>\n');
%----------------Inércia----------------%
fprintf(fid,'<has_body> %s </has_body>\n',body);         


%--------------------------------------ASA--------------------------------------%
%---------------Definições do asa---------------%
name_wing = 'Main Wing';    %Nomear a asa
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
fprintf(fid,'<Volume_Mass>  1.500</Volume_Mass>\n'); 
fprintf(fid,'</Inertia>\n'); 
%--------------Inércia-----------%

%------------Início das sections------------%
%------------------------------------Definições das sections------------------------------------------%

%definir o número de secções pretendidas
n_section = 20;
S = 0.897;
AR = 10.033;
b = sqrt(S * AR);
growth_rate = 1.05;
Y_asa = @(t)(((t*growth_rate^(n_section-t))/n_section)*(b/2));
%Y_asa = @(t)((t/n_section)*(b/2));            %função que define a posição de cada secção
C_r = (4 * S) / (pi * b);
corda_asa = @(y)(C_r * sqrt(1 - ((2 * y) / b)^2)); %função que define a corda em função da distância à raíz
twist_asa = 0;                  %função que define a twist em função da distância à raíz
diedro_1_asa = 0;

x_paineis_asa = 16;                                        %número de paineis em x
x_distribuicao_asa = 'COSINE';                             %escolher entre 'COSINE', 'UNIFORM', 'SINE' ou '-SINE'
y_paineis_asa = 10;                                         %número de paineis em x
y_distribuicao_asa = 'COSINE';                             %escolher entre 'COSINE', 'UNIFORM', 'SINE' ou '-SINE'
foil_asa = perfil;
%------------------------------------Definições das sections------------------------------------------%
fprintf(fid,'<Sections>\n'); 

for i = 0 : 19
    fprintf(fid,'<Section>\n');
    fprintf(fid,'<y_position> %.3f </y_position>\n',Y_asa(i));
    fprintf(fid,'<Chord> %.3f </Chord>\n',corda_asa(Y_asa(i)));
    fprintf(fid,'<xOffset> %.3f </xOffset>\n', (corda_asa(0) - corda_asa(Y_asa(i))) / 2);    
%-----------Diedros-----------%     
    fprintf(fid,'<Dihedral> %.3f </Dihedral>\n',diedro_1_asa);
    fprintf(fid,'<Twist> %.3f </Twist>\n',twist_asa);
    fprintf(fid,'<x_number_of_panels> %d </x_number_of_panels>\n',x_paineis_asa);
    fprintf(fid,'<x_panel_distribution> %s </x_panel_distribution>\n',x_distribuicao_asa);
    fprintf(fid,'<y_number_of_panels> %d </y_number_of_panels>\n',y_paineis_asa);
    fprintf(fid,'<y_panel_distribution> %s </y_panel_distribution>\n',y_distribuicao_asa);
    fprintf(fid,'<Left_Side_FoilName>%s</Left_Side_FoilName>\n',foil_asa);
    fprintf(fid,'<Right_Side_FoilName>%s</Right_Side_FoilName>\n',foil_asa);
    fprintf(fid,'</Section>\n');
end

fprintf(fid,'</Sections>\n');
%------------Fim das sections------------%
fprintf(fid,'</wing>\n');

fprintf(fid,'</Plane>\n');
fprintf(fid,'</explane>\n');
%--------------------------------------Fim do PLANE--------------------------------------%
fclose(fid);