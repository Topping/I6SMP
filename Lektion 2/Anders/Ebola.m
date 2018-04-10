function Ebola
%% Ebola example

Infected=rand(1,10);

for n=1:5
if Infected(n)<0.2
    Infected(n)=1;
else
    
   Infected(n)=0;
end
end

for n=6:10
    if Infected(n)<0.01
    Infected(n)=1;
    else
   Infected(n)=0;
end
end
day100=0;
for day=1:15

Infected_new = 0;
for n=1:length(Infected)
if Infected(n) == 1
    if rand(1,1)<0.5
    Infected_new = Infected_new+1; %Does the patient infect another or not?
    end
end
if sum(Infected) > 100 && day100 == 0
    day100=day;
end

end
Infected=[Infected ones(1, Infected_new)];


end

day100
sum(Infected)