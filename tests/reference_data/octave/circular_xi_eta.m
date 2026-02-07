function [xi, eta] = circular_xi_eta(center, radius, rotation, S)
% CIRCULAR_XI_ETA Compute boundary points and tangent vectors for a circle.
%   Replicates bcircle.xi_eta from +bounded_map/bcircle.m
%
%   [xi, eta] = circular_xi_eta(center, radius, rotation, S)
%
%   Inputs:
%     center   - complex center of the circle
%     radius   - radius of the circle
%     rotation - rotation angle (usually 0)
%     S        - parameter values (column vector)
%
%   Outputs:
%     xi  - boundary points (complex)
%     eta - tangent vectors (complex)

    reiS = radius * exp(1i*S + rotation);
    xi = center + reiS;
    if nargout > 1
        eta = 1i*reiS;
    end
end
