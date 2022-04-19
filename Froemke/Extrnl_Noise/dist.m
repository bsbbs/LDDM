%% plot distributions of R1 and R2 over iSTDP

%% loading data
filename = sprintf('LDDM_timeCourse_a%1.2f_b%1.2f_sgm%1.1fsinpt%0.3f_sims%i',a0,b0,sgm,sgmInput,sims);
output = fullfile(Simdir,[filename, '.mat']);
%%
pd = fitdist(Vcourse(:,1),'kernel','Kernel','normal');
x = -200:20:1200;
y = pdf(pd,x);
plot(x,y,'r-');