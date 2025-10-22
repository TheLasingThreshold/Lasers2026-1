% Ajuste para la obtención de los parámetros M^2, w0 y z0

%% Parámetros

% Longitud de onda en nanómetros
lambda = 800;

% Datos de directorio                
ruta = 'images\';         % Ruta local
nombre = 'image_';        % Nombre de la foto

% Tamaño de pixel
X_pixel = 5.3;
Y_pixel = X_pixel;

% Número de imágenes a ser analizados
n = 81;

% Intervalo de distancia de la platina que se tomó entre cada foto
res = 0.5; % mm

% Se define el mínimo valor como el valor más alto que se repite 10
% veces para eliminar el offset de los datos
repeatedValues = 10;

% Vector de posición en z (recorrido de la platina), en mm
Z = 0:res:(res*(n - 1));

% Arreglo que contiene los anchos del haz en x y en y
width = zeros(2,n);

%% Cálculo la desviación estándar de los marginales en x y y
for i = 1:n
    % Se carga la foto
    imagen = double(imread([ruta nombre sprintf('%d', i-1) '.jpg']));
    imagen = imagen(:,:,1);

    imagen = imagen - min(imagen(:));
    
    %{
    figure(1)
    surf(imagen', 'EdgeColor','none')
    view(2)
    %}
    
    % Centro de masa del perfil transversal de la foto
    [x_c,y_c] = centerOfMass2d(imagen, returnCoMIndexes=true);

    % Marginales normalizados, se integra en la dimensión ortogonal
    marginal_x0 = squeeze(sum(imagen, 1))/max(imagen(:));
    marginal_x0 = marginal_x0 - min(marginal_x0(:));
    marginal_y0 = squeeze(sum(imagen, 2))/max(imagen(:));
    marginal_y0 = marginal_y0 - min(marginal_y0(:));

    % Eliminación de datos en las alas después de que el marginal vale 
    % cero ambos sentidos ya que introduce mucho ruido

    % Se busca el máximo de cada marginal para utilizarlo para dividir el
    % marginal en dos partes, la parte derecha (valores antes del centro) y
    % la parte izquierda (valores después del centro)
    [~,xMax] = max(marginal_x0);
    [~,yMax] = max(marginal_y0);
    xLeftSide = marginal_x0(1:xMax);
    xRightSide = marginal_x0(xMax:end);
    yLeftSide = marginal_y0(1:yMax);
    yRightSide = marginal_y0(yMax:end);

    % Se calcula cuántas cuentas tienen cada valor
    [xLeftCounts, xLeftValues] = groupcounts(xLeftSide');
    [xRightCounts, xRightValues] = groupcounts(xRightSide');
    [yLeftCounts, yLeftValues] = groupcounts(yLeftSide);
    [yRightCounts, yRightValues] = groupcounts(yRightSide);

    % Con la variable repeatedValues, se define a partir de cuantos puntos
    % se define el valor de 0 de la curva. Para estos datos, corresponde en
    % general para 0 y algún otro valor adicional mayor. Puede suceder que
    % el valor adicional se repita más de lo ideal y llega a afectar el
    % valor de la desviación estándar.
    xOffset = max(max(xLeftValues(xLeftCounts > repeatedValues)),...
        max(xRightValues(xRightCounts > repeatedValues)));
    yOffset = max(max(yLeftValues(yLeftCounts > repeatedValues)),...
        max(yRightValues(yRightCounts > repeatedValues)));

    % Se resta el offset a los marginales y los valores negativos se
    % redefinen como 0
    marginal_x0 = marginal_x0 - xOffset;
    marginal_y0 = marginal_y0 - yOffset;
    marginal_x0(marginal_x0 < 0) = 0;
    marginal_y0(marginal_y0 < 0) = 0;
    %}

    % Puede suceder que después de que el marginal sea 0, en algún punto
    % más alejado del centro tome un valor diferente de 0 debido a ruido.
    % Estás líneas eliminan ese ruido lo cual para estos datos sólo se
    % elimina ruido.
    xLeftIndex = find(marginal_x0(1:xMax) == 0, 1, 'last');
    xRightIndex = find(marginal_x0(xMax:end)== 0, 1, 'first') + (xMax-1);
    yLeftIndex = find(marginal_y0(1:yMax) == 0, 1, 'last');
    yRightIndex = find(marginal_y0(yMax:end)== 0, 1, 'first') + (yMax-1);

    marginal_x0(1:xLeftIndex) = 0;
    marginal_x0(xRightIndex:end) = 0;
    marginal_y0(1:yLeftIndex) = 0;
    marginal_y0(yRightIndex:end) = 0;
    %}

    % Vectores de posición en micras los cuales se centran en (0,0)
    X = X_pixel*(1:length(marginal_x0))';
    X = X - x_c*X_pixel;

    Y = Y_pixel*(1:length(marginal_y0))';
    Y = Y - y_c*Y_pixel;
    
    %
    figure(2)
    plot(marginal_x0, '.')
    hold on
    plot(marginal_y0)
    hold off
    %}

    %% Cálculo de segundo momento
    segundo_mom_x = moments1d(X, marginal_x0', n=2, dim=1);
    segundo_mom_y = moments1d(Y', marginal_y0', n=2, dim=2);

    % D4sigma, el diámetro del haz dado por 4 veces la desviación estándar (2
    % veces la desviación estándar corresponde al radio)
    D2sigma_x = 2*sqrt(segundo_mom_x); % micras
    D2sigma_y = 2*sqrt(segundo_mom_y); % micras
    
    % Anchos del haz
    width(1,i) = D2sigma_x/1000; % mm
    width(2,i) = D2sigma_y/1000; % mm
    %}
end
%}

%% Ajuste de la ecuación
tic
% Función de ajuste para hallar M^2 (ecuación (4) de la guía)
model_Sieg = fittype(@(M2, w0, z0, z)...
    sqrt(w0^2 + M2^2*(lambda*1E-6*(z - z0)/(pi*w0)).^2),...
    'independent', 'z', 'coefficients', {'M2', 'w0', 'z0'}); % todo en mm

% Opciones para mejorar la convergencia del ajuste como límites inferiores
% de los coeficientes de ajuste y una condición inicial
options = fitoptions('gauss2', 'lower', [0,0,0], 'Start',...
    [1.00, min(width(:)), Z(floor(length(Z)/2))]);

% Ajuste de curvas en x y en y
[Fit_Sieg_x, gof_x] = fit(Z', width(1,:)', model_Sieg, options)
[Fit_Sieg_y, gof_y] = fit(Z', width(2,:)', model_Sieg, options)

toc

% Vector de de posición en z (recorrido de la platina). En este caso se
% define con una buena resolución para graficar los ajustes
z = Z(1):0.01:Z(end);

% Curvas de los ajustes para graficar
curve_x = model_Sieg(Fit_Sieg_x.M2, Fit_Sieg_x.w0, Fit_Sieg_x.z0, z); % x
curve_y = model_Sieg(Fit_Sieg_y.M2, Fit_Sieg_y.w0, Fit_Sieg_y.z0, z); % y

%% Plots
set(groot,'defaultAxesTickLabelInterpreter','latex');

% Vector de colores
colores = Colores_cb(2);

% Marginal en x
figure(2)
plot(X, marginal_x0, 'DisplayName', 'Marginal en $$x$$', 'MarkerEdgeColor', colores(2))
set(gca,'FontSize',13);
xlabel('$$x$$ ($$\mu$$m)','Interpreter','latex')
ylabel('Intensidad normalizada (u.a.)','Interpreter','latex')
legend('Interpreter','latex')

% Marginal en y
figure(3)
plot(Y, marginal_y0, 'DisplayName', 'Marginal en $$y$$', 'MarkerEdgeColor', colores(2))
set(gca,'FontSize',13);
xlabel('$$y$$ ($$\mu$$m)','Interpreter','latex')
ylabel('Intensidad normalizada (u.a.)','Interpreter','latex')
legend('Interpreter','latex')

% Ecuación de ajuste
figure(4)
scatter(Z, width(1,:), 'x', 'DisplayName', 'Datos en $$x$$', 'MarkerEdgeColor', colores(1))
hold on
plot(z, curve_x, 'DisplayName', 'Ajuste en $$x$$', 'Color', colores(1))
scatter(Z, width(2,:), 'x', 'DisplayName', 'Datos en $$y$$', 'MarkerEdgeColor', colores(2))
plot(z, curve_y, 'DisplayName', 'Ajuste en $$y$$', 'Color', colores(2))
hold off
set(gca,'FontSize',13);
xlabel('$$\textrm{Posici\''on de la platina}$$ (mm)','Interpreter','latex')
ylabel('$$\textrm{Ancho del haz D}2\sigma$$ (mm)','Interpreter','latex')
legend('Interpreter','latex', 'Location','northwest')

