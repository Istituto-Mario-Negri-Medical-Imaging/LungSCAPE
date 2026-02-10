function FA = computeFractionalAnisotropy(eigenvalues)
%COMPUTEFRACTIONALANISOTROPY Calculate fractional anisotropy from eigenvalues
%
%   FA = computeFractionalAnisotropy(eigenvalues)
%
%   Computes the fractional anisotropy (FA) metric from a set of eigenvalues.
%   FA is a measure of anisotropy ranging from 0 (isotropic) to 1 (highly anisotropic).
%   This is commonly used to characterize the shape of structures in 3D images.
%
%   Input:
%       eigenvalues - 3-element vector of eigenvalues [lambda1, lambda2, lambda3]
%                     where lambda1 >= lambda2 >= lambda3
%
%   Output:
%       FA - Fractional anisotropy value (scalar between 0 and 1)
%
%   Formula:
%       FA = (1/sqrt(2)) * sqrt((λ1-λ2)² + (λ2-λ3)² + (λ1-λ3)²) / sqrt(λ1² + λ2² + λ3²)
%
%   Example:
%       eigenvalues = [100, 50, 10];
%       FA = computeFractionalAnisotropy(eigenvalues);
%       % FA ≈ 0.73 (elongated structure)
%
%   See also: regionprops3, filterByAnisotropy

% Validate input
if length(eigenvalues) ~= 3
    error('computeFractionalAnisotropy:InvalidInput', ...
          'eigenvalues must be a 3-element vector');
end

% Handle edge case of zero eigenvalues
denominator = sqrt(eigenvalues(1)^2 + eigenvalues(2)^2 + eigenvalues(3)^2);
if denominator == 0
    FA = 0;
    return;
end

% Compute FA
numerator = sqrt((eigenvalues(1) - eigenvalues(2))^2 + ...
                 (eigenvalues(2) - eigenvalues(3))^2 + ...
                 (eigenvalues(1) - eigenvalues(3))^2);

FA = (1/sqrt(2)) * (numerator / denominator);

end
