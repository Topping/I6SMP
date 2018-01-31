permutations = 10000;
%%
% ES SOM F�RSTE KORT
% Hypotese: Chancen for at tr�kke et es er 4/52
% Ved X permutationer, vil antallet af es'er trukket f�rst v�re X * (4/52)
hypotese = 4 / 52;
expected = permutations * hypotese;
%Ace as first card
ace_count = 0;
for n=0:permutations
    cards = shufflecards();
    draw = cards(1);
    if(contains(draw.name,"Ace"))
        ace_count = ace_count + 1;
    end
end
a = sprintf("DRAWING ACE AS FIRST CARD");
b = sprintf("Expected %0.3f - Got %0.3f. %0.3f", expected, ace_count, expected/ace_count);
disp(a)
disp(b)

%%
% ES SOM F�RSTE KORT - BLACKJACK P� 2. KORT
% Hypotese, chancen for at tr�kke et es er 4/52
% Hypotese, chancen for at tr�kke blackjack p� 2. kort er 16/51
% Chancen for at f� blackjack efter at have trukket et es = (4/52) *
% (16/51)
P_es = 4/52;
P_ti = 16/51;
P_bj = P_es * P_ti;
expected_bj = P_ti * permutations;
bj_count = 0;
for n=0:permutations
    cards = shufflecards();
    draw = cards(1);
    if(contains(cards(1).name,"Ace") && cards(2).value == 10)
        bj_count = bj_count + 1;
    elseif(cards(1).value == 10)
        bj_count = bj_count + 1;
    end    
end
disp(expected_bj)
disp(bj_count)
disp(expected_bj/bj_count)

%%
P_1 = 4/52; % Chance for at tr�kke es f�rst
P_2 = 16/51; % Chance for at tr�kke 10 efter at have trukket es f�rst
P_3 = 16/52; % Chance for at tr�kke 10 f�rst
P_4 = 4/51; % Chance for at tr�kke es efter at have trukket 10 f�rst
P_bj = P_1 * P_2;
P_bj2 = P_3 * P_4;
expected_bj = (P_bj + P_bj2) * permutations;
bj_count = 0;
for n=0:permutations
    cards = shufflecards();
    if( (contains(cards(1).name,"Ace") && cards(2).value == 10) || (cards(1).value == 10 && contains(cards(2).name,"Ace")))
        bj_count = bj_count + 1;
    end    
end
disp("Expected " + expected_bj)
disp("Got " + bj_count)
disp(expected_bj/bj_count)

%%
% Strategi : Es = 1
% Tr�k 3 kort, check om samlet v�rdi er over 21
bust_count = 0;
for n=0:permutations
   cards = shufflecards();
   sum = 0;
   for i=1:3
      if(contains(cards(i).name, "Ace"))
         sum = sum + 1; 
      else
         sum = sum + cards(i).value; 
      end
   end
   if(sum > 21)
      bust_count = bust_count + 1; 
   end
end

disp("Bust chance: " + bust_count/permutations)

%%
bust_count = 0;
for n=0:permutations
   cards = shufflecards();
   sum = 0;
   for i=1:52
      if(contains(cards(i).name, "Ace") && (sum + 11 > 21))
         sum = sum + 1; 
      elseif(contains(cards(i).name, "Ace") && (sum + 11 <= 21))
         sum = sum + 11;
      else
         sum = sum + cards(i).value;
      end
      
      if(sum > 16)
         break 
      end
   end
   if(sum > 21)
      bust_count = bust_count + 1; 
   end
end

disp("Bust chance: " + bust_count/permutations)

%%
lose_count = 0;
for n=0:permutations
   cards = shufflecards();
   sum = 0;
   cards_taken = 1;
   for i=1:52
      if(contains(cards(i).name, "Ace") && (sum + 11 > 21))
         sum = sum + 1; 
      elseif(contains(cards(i).name, "Ace") && (sum + 11 <= 21))
         sum = sum + 11;
      else
         sum = sum + cards(i).value;
      end
      
      if(sum > 20)
          cards_taken = i;
         break 
      end
   end
   if(sum > 21)
      lose_count = lose_count + 1; 
      continue
   end
   dealers_sum = 0;
   for i=(cards_taken+1):52
      if(contains(cards(i).name, "Ace") && (dealers_sum + 11 > 21))
         dealers_sum = dealers_sum + 1; 
      elseif(contains(cards(i).name, "Ace") && (dealers_sum + 11 <= 21))
         dealers_sum = dealers_sum + 11;
      else
         dealers_sum = dealers_sum + cards(i).value;
      end
      
      if(dealers_sum >= sum)
         break 
      end
   end
   if(dealers_sum <= 21)
      lose_count = lose_count + 1; 
   end
end

disp("Win chance: " + (permutations - lose_count)/permutations)