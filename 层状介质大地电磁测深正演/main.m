clear all;
p = [100, 600];
h = 1800;
mu=(4e-7)*pi;
T=logspace(-3,4,40);
i=sqrt(-1);
k=zeros(size(p,2),size(T,2));
for N=1:size(p,2)
    k(N,:)=sqrt(-i*2*pi*mu./(T.*p(N)));
end
m=size (p,2);
y=-(i*mu*2*pi)./(T.*k(m,:));
%底层阻抗
for n=m-1:-1:1
    A=-(i*2*pi*mu)./(T.*k(n,:));
    B=exp(-2*k(n,:)*h(n));
    y=A.*(A.*(1-B)+y.*(1+B))./(A.*(1+B)+y.*(1-B)); %表面波阻抗
end
pc=(T./(mu*2*pi)).*(abs(y).^2); %视电阻率
ph=-atan(imag(y)./real(y)).*180/pi; %相位

figure();
subplot(1,2,1);
semilogx(T,pc,'k','LineWidth',1);
xlabel({'T/s';'(a) 视电阻率曲线'});
ylabel('$\rho_{\it{a}}$ / $\Omega\cdot m$','Interpreter','latex');

subplot(1,2,2);
semilogx(T,ph,'k','LineWidth',1);
xlabel({'T/s';'(b) 相位曲线'});
ylabel('度/(°)');
