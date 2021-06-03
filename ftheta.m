
function out=ftheta(theta,beta,lambdaRb) 
%calc orientation distribution

Q = integral(@(x)fthetann(x,beta,lambdaRb),0,2*pi);
out=fthetann(theta,beta,lambdaRb)/Q;

end


function  out=fthetann(theta,beta,lambdaRb)

%calc non-normalised orientation distribution

temp=lambdaRb*exp(-beta*sin(theta));
out=1./temp;

end

