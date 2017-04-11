function F = FiveGaussModel(x,freq)

% x(1) = gaussian amplitude
% x(2) = 1/(2*sigma^2)
% x(3) = centre freq of peak
% x(4-6) = gaussian 2 
% x(7-9) = gaussian 3
% x(10-12) = gaussian 4
% x(13-15) = gaussian 5
% x(16) = offset
% x(17) = slope
% x(18) = quadratic


%F = x(1)*sqrt(-x(2)/pi)*exp(x(2)*(freq-x(3)).*(freq-x(3)))+x(4)*(freq-x(3))+x(5);
F = x(1)*exp(x(2)*(freq-x(3)).*(freq-x(3))) + ...
    x(4)*exp(x(5)*(freq-x(6)).*(freq-x(6))) + ...
    x(7)*exp(x(8)*(freq-x(9)).*(freq-x(9))) +...
    x(10)*exp(x(11)*(freq-x(12)).*(freq-x(12))) +...
    x(13)*exp(x(14)*(freq-x(15)).*(freq-x(15))) + ...
    x(16) + x(17)*(freq - x(3)) + x(18)*(freq-x(3)).^2 ;
