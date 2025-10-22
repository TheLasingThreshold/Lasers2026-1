function moments1d = moments1d(xArray, yArray, parameters, centerOfMass)
% Cálculo del momento de orden n de una función con valores en x dados por
% xArray y valores en y dados por yArray. yArray puede ser una matriz y la
% dimensión en que se quiere calcular el momento de orden n debe ser en la
% dimensión dim y que también debe corresponder a la dimensión del vector
% xArray.

arguments
    xArray                      % Arreglo con valores en x
    yArray                      % Arreglo con valores en y
    parameters.n                % Órden del momento que se quiere calcular
    parameters.dim              % Dimensión tanto de xArray como de la 
    % dimensión en que se calcula el momento de orden n en yArray
    centerOfMass.centerOfMass   % Centro de masa de yArray que se conoce de antemano
end

if isfield(centerOfMass, 'centerOfMass')
    % En caso de que se calculó de antemano el centro de masa en la
    % dimensión dim. Un ejemplo es el caso de coordenadas cilíndricas con
    % simetría radial que tiene que ser 0, pero el vector radial r en el
    % código de la simulación sólo va de 0 hasta rMax. Si no se obliga a
    % que el centro de masa sea 0, se obtiene un valor erróneo.
    centerOfMass = centerOfMass.centerOfMass;
else
    % Si no se introduce este parámetro opcional, se calcula el centro de
    % masa de forma normal
    centerOfMass = centerOfMass1d(xArray, yArray, parameters.dim);
end

% Cálculo del momento de orden n
moments1d = trapz(xArray, (xArray - centerOfMass).^parameters.n .* yArray,...
    parameters.dim) ./ trapz(xArray, yArray, parameters.dim);

%% Función auxiliar
% Centro de masa en la dimensión dim
    function centerOfMassDim = centerOfMass1d(xArray, yArray, dim)
        centerOfMassDim = trapz(xArray, xArray.*yArray, dim)./...
            trapz(xArray, yArray, dim);
    end
end