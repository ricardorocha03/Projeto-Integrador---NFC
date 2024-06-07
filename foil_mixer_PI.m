%Foil Mixer Proj Integrador
format long

% Input
n_section = 11;

% Distribution ratio (bigger = more tip influence)
R = 0.5;
% T parameter for bezier
distribuition_density = 0.00001; 

% Foils (type XFLR5 .dat)
fsroot = 'Root_final.dat';
fstip = 'Tip_final.dat';

% Output directory
dir = '.\11PERFIS_PI_R0.5_3';

% \Input

% Read input foils
froot = fopen(fsroot,'r');
ftip = fopen(fstip,'r');
% Foil name
sroot = string(textscan(froot, '%s\n', 1));
stip = string(textscan(ftip, '%s\n', 1));
% Foil (x, y)
croot = textscan(froot, '%f %f');
ctip = textscan(ftip, '%f %f');
ogFoils = {cell2mat(croot); cell2mat(ctip)};
% Cleanup
fclose(froot);
fclose(ftip);

% foils: extraRoot{1,1}, intraRoot{1,2};
%         extraTip{2,1},  intraTip{2,2}
foils = {0,0;0,0};
for i = 1 : 2
    for j = 1 : length(ogFoils{i})
        if ogFoils{i}(j, 2) < 0
            foils{i, 1} = ogFoils{i}(1:j-1, :);
            foils{i, 2} = ogFoils{i}(j:end, :);
            break
        end
    end
end

% Root Ratio along sections

% Quadratic Bezier curve
B = @(t, a, b, c) (1-t).^2.*a + 2*(1-t).*t.*b + t.^2.*c;
t = 0 : distribuition_density : 1;

% Define [starting control ending] (point)
curve1 = {  [0; 0]   0 [1 - R; R]};
curve2 = {[1 - R; R] 0   [1; 1]  };
% Control dependent on ratio
if R > 0.5
    curve1{2} = [   0   ; 2*R-1];
    curve2{2} = [2*(1-R);   1  ];
else
    curve1{2} = [1-2*R;  0 ];
    curve2{2} = [  1  ; 2*R];
end

% Calculate points on the curve
bezier = transpose([B(t, curve1{1}, curve1{2}, curve1{3}) B(t, curve2{1}, curve2{2}, curve2{3})]);

% Search for points to fill in sectionRatio array
sectionRatio = zeros(1, n_section);
for i = 1 : n_section
    x = (i - 1)/(n_section - 1);
    for j = 1 : length(bezier) - 1 
        y = linInterp(bezier(j, :), x, bezier(j + 1, :));       
        if y ~= 0 
            sectionRatio(i) = y; 
        end
    end
end

% foil calc & write
mkdir(dir);
for i = 1 : n_section    
    filename = sprintf('perfil_PI_r0.3%03i.dat', i);
    fout = fopen(fullfile(dir, filename), 'wt');
    fprintf(fout, 'perfil_PI_r0.3%03i\n', i);
    for j = 1 : 2                
        half = newHalfFoil(foils{1,j}, foils{2, j}, sectionRatio(i));
        for k = 1 : length(half)
            fprintf(fout, '%f %f\n', half(k, 1), half(k, 2));
        end       
    end
    fclose(fout);
end

% Given 2 matrices(half foil), 2 columns (x, y), and a (tip) ratio, create new (half) foil
function hfoil = newHalfFoil(root, tip, ratio)
    hfoil = root;
    for i = 1 : length(root)
        x = root(i, 1);        
        for j = 1 : length(tip) - 1
            y_tip = linInterp(tip(j, :), x, tip(j + 1, :));
            if y_tip ~= 0
                y = (1 - ratio) * root(i, 2) + ratio * y_tip;
                hfoil(i, 1) = x;
                hfoil(i, 2) = y;
                break
            end            
        end
    end
end

% Linear Interpolation between 2 points, return 0 for invalid x
function y = linInterp(p1, x, p2)
    y = 0;
    x1 = p1(1); y1 = p1(2);
    x2 = p2(1); y2 = p2(2);
    if (x1 <= x && x <= x2) || (x2 <= x && x <= x1)
        m = (y2 - y1) / (x2 - x1);
        y = y1 + m * (x - x1);
    end
end
