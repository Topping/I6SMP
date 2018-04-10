function risk=bust_risk(cards_on_tabel, drawncards)
% Calculate bust risk, cards_on_tabel is a subset of drawncards
max_val = 21-sum(cards_on_tabel);
Px = cond_nth_card(drawncards);

% Calculate cdf
Xcdf=cumsum(Px);

% Risk of busting
if max_val < 10
risk = 1-Xcdf(max_val);
else
    risk = 0;
end
end