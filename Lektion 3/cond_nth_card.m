% based on the previous card values, input argumemt
function Pr_xn=cond_nth_card(drawncards)
% calculate the probability of the nth card
% based on the previous card values, input argumemt

NumberCards=length(drawncards);

%% Marginal pmf
numbers=1:10;

card_values=[1:10,10,10,10,1:10,10,10,10,1:10,10,10,10,1:10,10,10,10]';
%plot af pmf, we define that ace has a value of 1
[a b]=hist(card_values);
for n=1:NumberCards
    index=drawncards(n);
    a(index)=a(index)-1;
end

stem(numbers,a./sum(a))
Pr_xn=a./sum(a);
end