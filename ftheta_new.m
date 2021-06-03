
function out=ftheta_new(theta,beta,lambdaRb,kappa) 
%calc orientation distribution

Q = integral(@(x)fthetann(x,beta,lambdaRb,kappa),0,2*pi);
out=fthetann(theta,beta,lambdaRb,kappa)/Q;

end


function  out=fthetann(theta,beta,lambdaRb,kappa)

%calc non-normalised orientation distribution

temp=lambdaRb*exp(-beta*sin(theta));
temp=1./temp;
out=temp.*exp(kappa*cos(theta-pi/2));
end

