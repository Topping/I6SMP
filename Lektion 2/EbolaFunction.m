function days = EbolaFunction()
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

days = 0;
while sum(Infected) < 100 && days < 600
   Infected_new = 0;
   for n=1:length(Infected)
        if Infected(n) == 1
            if rand(1,1)<0.01
                Infected_new = Infected_new+1;
            end
        end
   end
   Infected=[Infected ones(1, Infected_new)];
   days = days + 1;
end
days;