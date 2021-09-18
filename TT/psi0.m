function y=psi0(x,p0,x0,alpha)
y=exp(-alpha/2*(x-x0)^2+1i*p0*(x-x0))*(alpha/pi)^0.25;
end

