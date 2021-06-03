
function out=von_mises(theta,kappa) 
%calc orientation distribution

Q = integral(@(x)fthetann(x,kappa),0,2*pi);
out=fthetann(theta,kappa)/Q;

end


function  out=fthetann(theta,kappa)

%calc non-normalised orientation distribution

out=exp(kappa*cos(theta-pi/2))
end

