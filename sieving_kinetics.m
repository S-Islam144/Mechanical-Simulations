clc;
clear all;
close all;

%sieving_kinetics solution coded by "Kh. Saidul Islam"

N = 10^5;
d_mean = 100;
d_std = 0.3*d_mean;
d = d_std*randn(1,N)+d_mean;

figure(1);
histogram(d,25);
hold on;
title('Particle Size Distribution');
xlabel('Particle size in micrometer');

Ns = 10^4;
d1 = 100;
d2 = 102;
ds = d1+(d2-d1)*rand(1,Ns);

figure(2);
histogram(ds,25);
hold off;
title('Sieve Opening Width Distribution');
xlabel('Sieve opening size in micrometer');

rho = 2000e-18;
P = 20;
d_input(1,:) = d(d>0);
d_fine = zeros(P,length(d_input(1,:)));
d_coarse = zeros(P,length(d_input(1,:)));
mass_fine = zeros(P,1);
mass_coarse = zeros(P,1);

for i=1:P
    for j=1:length(d_input(i,:))
        if d_input(i,j)~=0
            index=round(1+(Ns-1)*rand);
            if ds(index)>d_input(i,j)
                d_fine(i,j)=d_input(i,j);
                d_fine_product(1,j)=d_input(i,j);
            else
                d_coarse(i,j)=d_input(i,j);
            end
        end
    end
    mass_fine(i,1)=sum(rho*pi*(d_fine(i,:)).^3)/6;
    mass_coarse(i,1)=sum(rho*pi*(d_coarse(i,:)).^3)/6;
    d_input(i+1,1:size(d_coarse(i,:),2))=d_coarse(i,:);
end

figure(3);
plot(1:length(mass_fine),cumsum(mass_fine));
hold on;
plot(1:length(mass_coarse),mass_coarse);
hold off;
title('Calculated masses of fine and coarse particles');
xlabel('Number of Passes');
ylabel('Mass of fine and coarse Particles in micrograms');
legend('Mass of fine particles','Mass of coarse particles');


m_dia_cls=30;  
d_min=min(d);
d_max=max(d);
step_size=(d_max-d_min)/m_dia_cls;
d_size=d_min:step_size:d_max;    

d_size=d_min:step_size:d_max;
g=1;
for h=1:1:30
    d_size(g,h)=(d_size(g,h+1)+d_size(g,h))/2;
    d_mean_size(g,h)=d_size(g,h);
end


fine_mass= zeros(1,m_dia_cls);
for i=1:length(d_fine(end,:))
     if d_fine(end,i)~=0
     new_d_fine=(d_fine(end,i)-d_min+eps)/step_size;
     fine_mass(ceil(new_d_fine))=fine_mass(ceil(new_d_fine))+rho*pi*(d_fine(end,i).^3)/6;
     end
end
fine_mass_fraction= fine_mass./sum(fine_mass);
cumulative_fine_mass=cumsum(fine_mass_fraction);


coarse_mass= zeros(1,m_dia_cls);
for i=1:length(d_coarse(end,:))
     if d_coarse(end,i)~=0
     new_d_coarse=(d_coarse(end,i)-d_min+eps)/step_size;
     coarse_mass(ceil(new_d_coarse))=coarse_mass(ceil(new_d_coarse))+(rho*pi*(d_coarse(end,i).^3))/6;
     end
end
coarse_mass_fraction= coarse_mass./sum(coarse_mass);
cumulative_coarse_mass=cumsum(coarse_mass_fraction);


d_actual= nonzeros(d)';
feed_mass= zeros(1,m_dia_cls);
for i=1:length(d_actual)
    
     if (d_actual(1,i))~=0
     new_d_feed=(d_actual(1,i)-d_min+eps)/step_size;
     feed_mass(ceil(new_d_feed))=feed_mass(ceil(new_d_feed))+ rho*pi*(d_actual(1,i).^3)/6;
     end
end
feed_mass_fraction= feed_mass./sum(feed_mass);
cumulative_feed_mass=cumsum(feed_mass_fraction);


figure(4)
plot(d_mean_size,fine_mass_fraction);
ylim([0 1]);
xlim([d_min d_max]);
hold on;
plot(d_mean_size,coarse_mass_fraction);
hold on;
plot(d_mean_size,feed_mass_fraction);
hold off;
title('Particle Mass Distribution of Fine,Coarse and Feed Particles');
xlabel('Particle Diameter (micrometers)');
ylabel('Mass Fraction (%)');
legend('fine mass fraction','coarse mass fraction','feed mass fraction');


figure(5);
plot(d_mean_size,cumulative_fine_mass);
hold on;
plot(d_mean_size,cumulative_coarse_mass);
hold on;
plot(d_mean_size,cumulative_feed_mass);
ylim([0 1]);
xlim([d_min d_max]);
hold off;
title('Particle Cumulative Mass Distribution of Fine, Coarse and Feed Particles');
xlabel('Particle Diameter (micrometers)');
ylabel('Cumulative Mass Fraction (%)');
legend('fine mass fraction','coarse mass fraction','feed mass fraction');


m_dia_cls=30;  
d_minm= 0;
d_maxm=max(d);
step_size= (d_maxm-d_minm)/m_dia_cls;
d_size_separate_function=d_minm:step_size:d_maxm;
d_size_separate_function=d_minm:step_size:d_maxm;
k=1;
for p=1:1:30
    d_size_separate_function(k,p)=(d_size_separate_function(k,p+1)+d_size_separate_function(k,p))/2;
    d_mean_size_separate_fuction(k,p)=d_size_separate_function(k,p);
end

fine_mass_sep= zeros(1,m_dia_cls);
for i=1:length(d_fine_product(1,:))
    
     if d_fine_product(1,i)~=0
     new_d_fine_sep=(d_fine_product(1,i)-d_minm+eps)/step_size;
     fine_mass_sep(ceil( new_d_fine_sep))=fine_mass_sep(ceil( new_d_fine_sep))+ rho*pi*( d_fine_product(1,i).^3)/6;
     end
end

coarse_mass_sep= zeros(1,m_dia_cls);
for i=1:length(d_coarse(end,:))
    
     if d_coarse(end,i)~=0
     new_d_coarse_sep=(d_coarse(end,i)-d_minm+eps)/step_size;
     coarse_mass_sep(ceil(new_d_coarse_sep))=coarse_mass_sep(ceil(new_d_coarse_sep))+ (rho*pi*(d_coarse(end,i).^3))/6;
     end
end

Separate_Fun=coarse_mass_sep./(coarse_mass_sep+fine_mass_sep+eps); 

figure(6);
plot(d_mean_size_separate_fuction,Separate_Fun);
xlim([d_minm d_maxm]);
ylim([0 1]);
title('Separation Function Plot');
xlabel('30 particle classes');



l=1;
count=0;
while Separate_Fun(l) <= 0.50
       count=count+1;
       l=l+1;
end
d_50=interp1([Separate_Fun(count) Separate_Fun(count+1)],[d_size_separate_function(count) d_size_separate_function(count+1)],0.50,'linear');

m=1;
count=0;
while Separate_Fun(m) <= 0.25
      count=count+1;
      m =m+1;
end
d_25=interp1([Separate_Fun(count) Separate_Fun(count+1)],[d_size_separate_function(count) d_size_separate_function(count+1)],0.25,'linear');

n=1;
count=0;
while Separate_Fun(n) <= 0.75
    count=count+1;
    n=n+1;
end
d_75=interp1([Separate_Fun(count) Separate_Fun(count+1)],[d_size_separate_function(count) d_size_separate_function(count+1)],0.75,'linear');

Selectivity=d_25/d_75;         
disp(Selectivity);
