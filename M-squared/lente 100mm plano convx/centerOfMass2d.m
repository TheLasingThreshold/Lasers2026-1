function [xCoM, yCoM] = centerOfMass2d(A, arrays, outputConditions)
% Calculo del momento de orden 'order' dado por la referencia libro Diels
% ecuacion (1.47)
arguments
    % Arreglo de dos dimensiones que se quiere calcular su centro de masa
    A (:,:) {mustBeNumericOrLogical}                % Arreglo con dimensiones (N,M)

    % Arreglos en x y en y del arreglo A
    arrays.xArray (:,1) {mustBeNumericOrLogical}    % con dimensiones (N,1)
    arrays.yArray (1,:) {mustBeNumericOrLogical}    % con dimensiones (1,M)

    % Indica si se quiere obtener los índices del vector xArray y yArray 
    % que corresponden al centro de masa en vez de los valores reales
    outputConditions.returnCoMIndexes {mustBeNumericOrLogical} = false 
end

% Condicional de si se introduce algún arreglo auxiliar
hasAdditionalArrays = all(isfield(arrays, {'xArray', 'yArray'}));

% Si no se introducen los respectivos arreglos en x y en y, se trabaja en
% términos de los índices del arreglo A
if hasAdditionalArrays
    xArray = arrays.xArray;
    yArray = arrays.yArray;
else
    xArray = (1:size(A,1))';    % (size(A,1),1)
    yArray = 1:size(A,2);       % (1,size(A,2))
end

% Las diferenciales se contienen en el resultado de trapz al dar el vector
% sobre la dimensión de integración, como se realizó en este caso

% Masa del arreglo
mass = trapz(yArray, trapz(xArray, A, 1), 2);

% Cálculo de los momentos
xCoM = trapz(yArray, trapz(xArray, xArray .* A, 1), 2)/mass;
yCoM = trapz(yArray, trapz(xArray, yArray .* A, 1), 2)/mass;

if outputConditions.returnCoMIndexes
    % En caso de que se introduzcan los arreglos en x y en y, se busca el
    % índice de cada uno que se acerca más al centro de masa en esa
    % dimensión. Esto es importante para el caso en que xArray o yArray no
    % tienen una resolución dx,dy constante, como por ejemplo r en
    % coordenadas cilíndricas.
    [~, xCoM] = min(abs(xArray - xCoM));
    [~, yCoM] = min(abs(yArray - yCoM));
else
    % Se realiza este cálculo para pasar los valores de notación científica
    % a flotante en el caso en que se trabaje con indices.
    xCoM = round(xCoM);
    yCoM = round(yCoM);
end

end
