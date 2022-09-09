% sigma = 0.1;
% x = linspace(-1,5,1000);
% h = mean(diff(x));
% f = sin(x) + randn(size(x))*sigma;
% plot(f)
% 
% sig = estimate_sigma(f,h);
% 
% disp(abs(sig-sigma)/sigma)

function sig = estimate_sigma(f,h)

    Ih = (f(4:end-1)-f(2:end-3))/h/2;
    I2h = (f(5:end)-f(1:end-4))/h/4;    
    sig = sqrt(8/5)*h*rms(Ih-I2h);

    if sig<0.001
        I4h = (f(9:end)-f(1:end-8))/h/8;
        Ih_1 = 4/3*(Ih - 1/4*I2h);
        I2h_1 = 4/3*(I2h(3:end-2) - 1/4*I4h);
        sig = sqrt(576/714)*h*rms(Ih_1(3:end-2)-I2h_1);
    end 
end