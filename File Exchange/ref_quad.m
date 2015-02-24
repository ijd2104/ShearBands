% Generate two-point Gaussian quadrature over the reference element 

function [ref_quad_pos, quad_weights] = ref_quad()

% Quadrature nodes
ref_quad_pos = [0.5 - 1 / (2 * sqrt(3)), 0.5 + 1 / (2 * sqrt(3))];

% Quadrature weights
quad_weights = [0.5, 0.5];

end