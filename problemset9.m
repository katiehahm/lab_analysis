figure
t=linspace(10800,40000)
T = 950+428.*exp((-1.435*10^-4).*(t-10800))
plot(t,T)
xlabel('Time (Seconds)')
ylabel('Temperature (K)')
