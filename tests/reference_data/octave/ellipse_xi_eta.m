function [xi, eta] = ellipse_xi_eta(center, semi_maj, semi_min, rotation, S)
% ELLIPSE_XI_ETA Compute boundary points and tangent vectors for an ellipse.
%   Replicates bellipse.xi_eta from +bounded_map/bellipse.m
%
%   [xi, eta] = ellipse_xi_eta(center, semi_maj, semi_min, rotation, S)
%
%   Inputs:
%     center   - complex center of the ellipse
%     semi_maj - semi-major axis length
%     semi_min - semi-minor axis length
%     rotation - rotation angle
%     S        - parameter values (column vector)
%
%   Outputs:
%     xi  - boundary points (complex)
%     eta - tangent vectors (complex)

    cosS = cos(S); sinS = sin(S);
    cosR = cos(rotation); sinR = sin(rotation);
    a_cosR = semi_maj * cosR; b_cosR = semi_min * cosR;
    a_sinR = semi_maj * sinR; b_sinR = semi_min * sinR;

    xi = a_cosR*cosS - b_sinR*sinS ...
         + 1i*(b_cosR*sinS + a_sinR*cosS) ...
         + center;

    if nargout > 1
        eta = -a_cosR*sinS - b_sinR*cosS ...
              + 1i*(b_cosR*cosS - a_sinR*sinS);
    end
end
