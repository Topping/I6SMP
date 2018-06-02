% %% Stokastiske Variable
% %Opgave 1 - Tæthadsfunktion(PDF) eller fordelingfunktion (CDF)
% 
%     %Brug differentiation F(x) i wolframalpha for at få PDF(Tæt) f(x)
% 
%     %Brug intergral f(x) i wolframalpha for at få CDF(fordel) F(X)
%     %for x gående mod -inf er F(x)=0 og gående mod inf F(x)=1
%     %For det totale interval fra -inf til inf er F(x)=1
% 
%     %Skal gøres og opskrives for hele intervallet
% 
% 
% %Opgave 2 For hvilken værdi af k er f(x) en gyldig tæthedsfunktion?
% 
%     %Summen af f(x) skal være 1
%     % tag integralet af det derefter
%     % solve k+k..=1 i wolframalpha så bliver k=1/?
% 
% %Opgave 3 vis at f(x) og F(x) er en gyldige
% 
%     %For f(x) i hele intervalet skal det være lig med 1,
%     %altså intergral -inf til inf f(x)=1, 
%     %gør det for det område hvor k er med
% 
%     %fx F(1)=1 solve k i denne ligning
% 
% %Opgave 4 Skitsér tæthedsfunktionen.
% 
%     x = [3 4 5];
% 
%     format rat
%     k = 1/3;
% 
%     fX = [k k k];
% 
%     figure(1);
%     bar(x,fX);
% 
% %% Opgave 5 Sandsynlighed for interval
%     %eksempel 1 Pr(-1/2<=X<=1/2) og Pr(X>0)
%     %Kan ikke køre
%     k = 5/2;
% 
%     F(X) = 1/5*k(x^5+1)
% 
%     Pr(-1/2 <= x <= 1/2) = F(1/2)-F(-1/2) = ...
%         1/5*5/2*(1/2^5+1) - 1/5*5/2*((-1/2)^5+1) = 1/32
% 
%     Pr(x > 0) = 1  
%     F(0) = 1/5*5/2*(0^5+1) = 1/2
% 
%     %eksempel 2 Pr(X <= 1/2) og Pr(X > -1/4)
%     k = 3/2
%     
%     F(x) = 1/3*k(x^3+1)
%     
%     Pr(X <= 1/2) = F(1/2) = 1/3*3/2(1/2^3+1) = 1/2*(1/8 + 1) = 9/16
%     
%     Pr(X > -1/4) = F(-1/4) = 1/2*(-1/4^3 + 1)
%     %eller
%     Pr(X > -1/4) = 1
%     Pr(X <= -1/4) = 1 %sand værdi for den er over -1/4 er 1, 
%     %det kan så trækkes fra sand i punktet? Dette skal gøres de den ikke
%     %sættes lig med ligesom den anden 
%     %F(x <= 1/2)  F(-1/4) = 1 - 1/2*(-1/4^3 + 1) = 1 - 1/128 - 1/2 = (128+1-64)/128 = 68/128
%     
%     
% %% Opgave 6 - Beregn middelværdien og variansen af X
% 
%     %E[X] = integral -1 to 1 x*f(x) 
%     %Brug wolframalpha "integral -1 to 1 x*f(x)"
%     
%     %mid = E[X^2] = integral -1 to 1 x^2*f(x)
%     %Brug wolframalpha "integral -1 to 1 x^2*f(x)"
%     
%     %Var(x) = E[X^2] - (E[X])^2
%     
% %% Opgave 7 - 2 stokastiske variable X og Y værdier Kovarians og korrelation
% 
%     %Eksempel E15 diskret
%     clear
%     clc
% 
%     XMatrix = [1 2 3];
%     YMatrix = [5;6;7;8];
%     Matrix = [0 1/12 0;2/12 0 2/12;2/12 1/12 2/12;0 2/12 0];
% 
%     X = sum(Matrix);
%     Y = sum(Matrix,2);
%     format rat
%     X
%     Y
% 
%     % E[X]
%     Ex = sum(XMatrix .* X) %middelværdi y
% 
%     % E[Y]
%     Ey = sum(YMatrix .* Y) %middelværdi y
% 
%     % E[X*Y]
%     Exy = sum(XMatrix .* X)*sum(YMatrix .* Y)
% 
%     % E[X^2]
%     Ex2 = sum(XMatrix.^2 .* X) 
% 
%     % E[Y^2]
%     Ey2 = sum(YMatrix.^2 .* Y)
% 
%     Vx = Ex2 - Ex^2; % variansen af x
%     Vy = Ey2 - Ey^2; % variansen af y
%     CovXY = Exy-Ex*Ey;
%     RHO = CovXY/sqrt(Vx*Vy); % korrelationskoefficienten RHO
% 
%     %hvis Rho = 0 så er X og Y ukorrelerede
% 
%     % Hvis f(x)*f(y) = f(x,y) så er de identisk ellers ikke identisk 
% 
% 
% %% Stokastiske processer
% %% Opgave 1 skitse af fordelinger/realistationer plot
clear
clc
clear

% En kontinuer stokastisk process er givet ved:
% X(t) = w + 5
% Hvor w er normalfordelt efter w~N(5,1).

% 1) Skitsér fem realisationer af processen X(t) mellem t ? [0; 7]. 
% Brug enGauss-generator, 
% det kan evt. være matlabs indbyggede generator, randn().
% Angiv desuden hvordan de fem realisationer er fremkommet.

w = randn(1,5)+5;
t=0:7;

x1=w(1)+5;
x2=w(2)+5;
x3=w(3)+5;
x4=w(4)+5;
x5=w(5)+5;

figure(1)
plot(t,ones(1,length(t))*x1)
grid
hold on;
plot(t,ones(1,length(t))*x2)
plot(t,ones(1,length(t))*x3)
plot(t,ones(1,length(t))*x4)
plot(t,ones(1,length(t))*x5)

%     %% Eksempel F16 re kontinuert normalfordelt for X(t) = w(n)
%     t = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8, 0.9];
%     X = randn(1,length(t))+t;
%     figure(2)
%     plot(t,X,'--x')
%     grid on
%     xlabel('tid')
%     ylabel('X')
% 
%     %% eksempel E15 diskret Uniform fordeling hvor X(n) = w(n), 10 samples fra 0 til 10
%     for n = 1:10
% 
%         X(n) = rand*10;
% 
%     end
%     figure(3)
%     plot(X,'--x','MarkerSize',10)
% 
%     %% eksempel 2 F13 re Uniform A og B og X(t) kontinuert
%     
%     % ReeksamenF2013 opg. 2.1
% N = 3;
% 
% % A~U(-1,1)
% A = rand(1,N)*2 - 1;
% % A = [ -0.5 0 0.5 ];
% 
% % B~U(-2,2)
% B = rand(1,N)*4 - 2;
% % B = [ -1.5 0.1 1.4 ];
% 
% t = 0:0.1:5;
% 
% figure
% hold on
% for i = 1:N
%     X = A(i)*t+B(i);
%     plot(t,X)
% end
% hold off
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% % Eksempel med normalfordeling
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% N = 3;
% 
% % A~N(0,1)
% mu = 0;
% sigma = 1; % Standardafvigelsen (ikke variansen!)
% A = randn(1,N);
% 
% % B~N(2,4)
% mu = 2;
% sigma = sqrt(4); % Standardafvigelsen (ikke variansen!)
% B = randn(1,N)*sigma + mu;
% 
% t = 0:0.1:5;
% 
% figure
% hold on
% for i = 1:N
%     X = A(i)*t+B(i);
%     plot(t,X)
% end
% hold off
% 
%     
%     %% eksempel E15 re diskret Normalfordeling i.i.d Gaussisk fordelt 
%     %hvor X(n) = w(n), 10 samples 
%     for n = 1:10
% 
%         X(n) = randn+4;
% 
%     end
%     figure(4)
%     plot(X,'--x','MarkerSize',10)
% 
%     
%     %% eksempel Sinus funktion 
%     
%     N = 3;
% 
%     % A~U(-1,1)
%     Theta = rand(1,N)*pi*2 - pi;
%     % A = [ -0.5 0 0.5 ];
% 
%     % B~U(-2,2)
%     Y = randn(1,N);
%     % B = [ -1.5 0.1 1.4 ];
% 
%     t = 0:0.1:10;
% 
%     figure
%     hold on
%     for i = 1:N
%         X = 2*sin(1*t+Theta(i))+Y(i);
%         plot(t,X)
%     end
%     hold off
%     
% %% Opgave 2 - Ensemble middelværdi og varians
% 
%     %Eksempel X(t) = w + 5, kontinuert normalfordelt w  N(5,1)
%     % Ensemble middelværdi
%     E[X(t)] = E[w] + E[5] %middelværdi fra for w og for offset 5 
% 
%     %Ensemble varians
%     Var[X(t)] = var[w] + var[5] = 1 + 0
% 
% 
% %Opgave 3 - En realisation udvalgt middelværdi og varians
% 
%     x(t) = 5.7 + 4; %udvalgt værdi
% 
%     %Middelværdi
%     E_xt = mean(5.7 + 4)
% 
%     %Varians
%     Var_xt = var(5.7 + 4)
% 
% 
% %% Opgave 4 - WSS og ergodisk??
% 
%     % X(t) er stationær WSS i den brede forstand, når ensemble middelværien
%     % E[X(t)] og ensemble variansen Var[X(t)] ikke er tidsafhængige. Altså at
%     % de ikke ændre værdi over tid. Eksempel F16 stationær - eksempel F16 re
%     % ikke stationær.
% 
% 
%     % X(t) er ergodisk når en realisation siger noget om ensample middelværdi og
%     % variansen, altså de er tæt på den samme værdi. Eksempel F16 og F16 re er
%     % ikke ergodisk, men eksempel E15 er ergodisk!
% 
% 
% 
% 
% %% Opgave 5 - Opstil ligningen til bestemmelse af Autocorrelationen og udregning
% 
%     % RX(t1,t2) = E[X(t1)*X(t2)]
% 
%     
%  % Eksempel 1   X(t) = A*t + B    for  t >= 0
%     % RX(t1,t2) = E[X(t1)*X(t2)] = E[(A*t1 + B)*(A*t2 + B)] ...
%     % = t1*t2*E[A^2] + t1*E[A]*E[B] + t2*E[A]*E[B] + E[B^2] ...
%     % = t1*t2*Var(A) + t1*0*0 + t2*0*0 + Var(B) ?= t1*t2*1/3*4/3
% 
%     % Var(X) = E[X^2] - E[X]^2. Da A og B er uniform fordelt og symmetrisk
%     % omkring nul, så er E[A] = E[B] = 0. Derfor har vi Var(A) = E[A^2] 
%     % og Var(B) = E[B^2] 
%     
%     % Processen er ikke stationær da autokorrelationsfunktionen ikke er en
%     % funktion af t1 - t2. Derfor er den heller ikke ergodisk da den ikke
%     % er stationær.
% 
%     
%  % Eksempel 2 E15 - diskret i.i.d uniform AutoCorr Rxx(tau) for tau = 0,1,2,3 og X(n) = w(n)  
%     
%     % Var(X) = E[X^2] - E[X]^2  =>  E[X^2] = Var(X) + E[X]^2
%     
%     % Rxx(0) = E[w(n)*w(n)] = E[X^2] = Var(w(n)) + E[w(n)]^2 = 
%     % 25/3 + 5^2 = 100/3
%     
%     % w(n) og w(n) er angivet som uafhængige (i.i.d)
%     % Rxx(1) = E[w(n)*w(n+1)] = E[w(n)]*E[w(n+1)] = 5*5 = 5^2 = 25
%     
%     % Rxx(2) = E[w(n)*w(n+2)] = E[w(n)]*E[w(n+2)] = 5*5 = 5^2 = 25
%     
%     % Rxx(3) = E[w(n)*w(n+3)] = E[w(n)]*E[w(n+3)] = 5*5 = 5^2 = 25
%     
% %% sansynlighedsrening
% 
% %% Opgave 1 - sandsynlighed for et hændelse
% 
%  % Eksempel 1 divid - Hvis der blev født 29.785 drenge og 28.131 piger i 2012, hvad
%     % er sandsynligheden for hændelse A = at en gravid fødte en pige i 2012?
% 
%     Pr_A = 28.131/(28.131+29.785) %sandsynlighed for at det bliver en pige 0.4857
% 
%  % Eksempel 2 total - Hvad er den totale sandsynlighed for at et af børnene i studiet 
%     % begik alvorlig kriminalitet indenfor de næste 10 år? 
% 
%     %A = flyttet mere end en gang
%     %B = Begået kriminalitet
% 
%     %Pr(B|A) = sandsynlighed for at begå kriminalitet når de er flyttet mere 
%       %end en gang.
%     %Pr(B|A_) = sandsynlighed for at begå kriminalitet når de er flyttet en
%       %gang eller 0
%     %Pr(A) = sandsynlighed for at børnene i studiet tilhørte gruppe der var
%       %flyttet flere gange.
%     %Pr(A_) = 1- Pr(A) = sandynlighed for at børnene i studiet tilhørte gruppen 
%       %der IKKE var flyttet flere gange.
%     %Pr(B,A) = samlet sandsynlighed for at børnene i studiet der har flyttet 
%       %flere gange begår kriminalitet.
%     %Pr(B,A_) = samlet sandsynlighed for at børnene i studiet der IKKE har flyttet 
%       %flere gange begår kriminalitet.
% 
%     %Total sandsynlighed for begået kriminalitet af børnene Pr_B
%     %Pr_B = Pr(B,A) + Pr(B,A_) = Pr(A)*Pr(B|A) + Pr(A_)*Pr(B|A_) = ...
%       %0.31*0.06 + 0.69*0.03 = 0.0393
% 
% %% Opgave 2 - Bayers regl
% 
%     %Bayers bruges når vender sandsynligheden omvendt, altså hændelsen er sket 
%     %som man regnede på før. Findes typisk udfra samlet sandsynlighed Pr(B) og
%     %sandsynlighed for en bestemt gruppe har gjort noget Pr(B|A)*Pr(A).
% 
%  %Eksempel fra opgave 1 med børne kriminalitet
%     %Hvis et barn fra studiet har begået alvorlig kriminalitet indenfor de 10 år 
%     %studiet rakte sig over, hvad er sandsynligheden for at barnet tilhørte
%     %gruppen, der havde flyttet mere end én gang?
% 
% 
%     Pr(A|B) = Pr(B|A)*Pr(A)/Pr(B)
%     Pr_AB = 0.06*0.31/0.0393 % = 0.4733, dvs der er 47% sandsynlighed for at barnet 
%                              % tilhørte gruppen der havde flyttet mere end en gang
% 
% %% Opgave 3 - Simultan sandsynlighed
% 
%  % Eksempel 1 - Hvis hændelse B er at spar es er blandt de syv kort. 
%     % Hvad er den simultane sandsynlighed for hændelserne A og B? 
% 
%     % A = at hjerter konge er blandt syv kort
%     % B = at spar es er blandt de syv kort
% 
%     %Hændelse B er den samme som Pr(A) = 1/52
%     Pr(B) = 1/52;
% 
%     Pr(A|B) = 6*1/51; %sandsynlighed for at spar es er blevet trukket og hjerter konge trækkes 
% 
%     %Simultan sandsynlighed
%     Pr(A,B) = Pr(A|B)*Pr(B)
%     Pr_AB = 6/51*7/52 %Simultan sandsynlighed 0.0158 spar es og hjerter konge bliver trukket
% 
%     %DA Pr(A)*Pr(B)=0.0181 og Pr(A,B) = 0.0155. SÅ ER Pr(A)*Pr(B) IKKE DET
%     %SAMME SOM Pr(A,B) OG DERVED IKKE UAFHÆNGIGE!!!!!!!
% 
% %% Opgave 4 - Kombinationer
% 
%  %Eksempel - Hvor mange forskellige kombinationer af 7 kort kan der trækkes 
%     %fra et spil kort af 52 kort?  
% 
%     % Antal kombinationer af 7 kort. Uden orden og uden tilbagelægning:
%     n = 52; %antal kort
%     k = 7;  %antal kort der skal trækkes
%     Combi1 = factorial(n)/(factorial(n-k)*factorial(k))
% 
%     % Orden og med tilbagelægning
%     % Så skulle kortene trækkes fra en ende af (orden) og lægges tilbage
%     % i bunken hver gang (tilbagelægning)
%     Combi2 = n^k
% 
%     % Orden og uden tilbagelægning
%     %Så skulle kortene trækkes fra en ende af (orden) og kortene lægges IKKE
%     %tilbage i bunken (Uden tilbagelægning)
%     Comb3 = factorial(n)/factorial(n-k)
% 
%     % Uden orden og med tilbagelægning
%     %kortene trækkes fra en random fra bunken (uden orden) og kortene lægges
%     %tilbage i bunken (tilbagelægning)
%     Comb4 = factorial(n+k-1)/factorial(n-1)*factorial(k)
% 
% %% Statistik
% 
% %% Opgave 1 - hypotese test
% 
% %Nul hypotesen er den vi tester for, hvilket typisk går ud fra at fx
% %middelværdierne er ens for begge datasæt.
% 
% %H0: mid1 = mid2 
% 
% %Den alternative er omvendt af nul hypotesen er at der er forskel på den pågældende
% %middelværdi for de to datasæt.
% 
% %H1: mid1 =| mid2 %
% 
% 
% %% Opgave 2 - Hvilken test skal bruges? (KIG I Testkataloger.pdf)
% 
% %Binomial spredning er en Z-test
% 
% %Poisson er en Z-test
% 
% %Når varianserne er ukendte og data er normalfordelt, 
% %bruges en T-test på forskellen af middelværdi.
% %Når det er kendte varianser bruges en z-test.
% 
% %Parret eller uparret? - Når de to data har haft de samme forhold er de
% %parret ellers er de uparret. 
% %Altså hvis man har målt data på nogle mennesker før og efter intagen af medicin så er det parret.
% 
% 
% %% Opgave 3 - Estimer middelværdier og standard afvigelser (KIG I Testkataloger.pdf)
% 
% %For Ukendt varians - Comparing two means (Unknown variance)
% %Estimation af forskellen på middelværdierne mid1 og mid2:
% deltaHat = abs(mid1-mid2)
% 
% %x1 for data 1, x2 for data 2
% %n1 = antal af data i sæt 1
% n1 = length(x1)
% %n2 = antal af data i sæt 2
% n2 = length(x2)
% 
% %estimater varians for data 1
% s_1 = sqrt(1/(n1-1)*sum(x1-mean(x1))^2
% 
% %estimator varians for data 2
% s_2 = = sqrt(1/(n2-1)*sum(x2-mean(x2))^2
% 
% 
% %Esitmator for samlet varians (pooled variance)
% S2 = (1/(n1+n2-2)) * ...
%     (((n1-1)*s_1)+((n2-1)*S_2));
% 
% S = sqrt(S2)
% 
% %% Opgave 3 eksempel to
% 
% %For Ukendt varians - Comparing two means (Unknown variance)
% %Estimation af forskellen på middelværdierne mid1 og mid2:
% deltaHat = abs(1.68-1.78)
% 
% %x1 for data 1, x2 for data 2
% %n1 = antal af data i sæt 1
% n1 = 19
% %n2 = antal af data i sæt 2
% n2 = 35
% 
% %estimater varians for data 1
% %s_1 = sqrt(1/(n1-1)*sum(x1-mean(x1))^2
% S_1 = 0.1
% S_2 = 0.2
% %estimator varians for data 2
% %s_2 = = sqrt(1/(n2-1)*sum(x2-mean(x2))^2
% 
% 
% %Esitmator for samlet varians (pooled variance)
% S2 = (1/(n1+n2-2)) * ...
%     (((n1-1)*S_1)+((n2-1)*S_2));
% 
% S = sqrt(S2)
% 
% %% Opgave 4 - Ukendt varians Test og signifikansniveau på 0.05 (KIG I Testkataloger.pdf)
% 
%  %Eksempel Comparing two means (Unknown variance)
% 
%     %uparret t-test
%     t = (mid1-mid2) / (S*sqrt((1/n1)+(1/n2))) %ex. = 3.75
% 
%     P = 2*(1-tcdf(abs(t),n1+n2-2)) %ex. = 0.0015
%     
%     %Hvis P er mindre et signifikansniveau på 0.05 afvises nul hypotesen,
%     %Hvis P er større afvises null hypotesen IKKE. %ex. afvises
% 
% %% Opgave 5 - 95% konfidensinterval for hældning (KIG I Testkataloger.pdf)
% d_minus = mid1-mid2-t0*S*sqrt(1/n1+1/n2)
% 
% d_plus = mid1-mid2+t0*S*sqrt(1/n1+1/n2)
% 
% %% Eksempel F16
% d_minus = 1.68-1.78-2.0066*0.4067*sqrt(1/19+1/35)
% 
% d_plus = 1.68-1.78+2.0066*0.4067*sqrt(1/19+1/35)
% 
% %% Opgave 6 - Lineær regression og 95% konfidensinterval for hældning (KIG I Testkataloger.pdf)
%     clear
%     Antal = [5562 4357 3471 3078 2309 1285 969 602 238 268]; %y-akse
%     Aar = [1901 1911 1921 1931 1941 1951 1961 1971 1981 1991]; %x-akse
% 
%     y = Antal;
%     x = Aar;
% 
%     %Antal data
%     n = length(x);
% 
%     %alpha (hældning/slobe) og beta (skærring på y-aksen/interception)
%     %De kan være sat omvendt i testkataloget
%     Alpha = sum((x - mean(x)).*(y - mean(y)))/sum((x - mean(x)).^2);
%     Beta = mean(y) - Alpha*mean(x);
% 
%     %Lineær regression formel
%     LinReg = Alpha * x + Beta;
% 
%     figure(1)
%     scatter(x,y)
%     hold on
%     grid on
%     plot(x,LinReg)
%     title('Lineær regression')
%     xlabel('x')
%     ylabel('y')
% 
%     %Residualtegning giver spredningen af data fra linien
%     %Hvis data er spredt tilfældigt omkring nul og ikke afhængig af X
%     %Så er betingelserne for lineær regression opfyldt!
%     Res = y - LinReg; %Lineær regression trækket fra data
%     figure(2)
%     plot(x, Res,'.', [min(x) max(x)], [0 0])
%     title('Residualtegning')
%     xlabel('x')
%     ylabel('Residual')
% 
%     %QQ plot kan også laves til er supplere om det er normaltfordelt
%     figure(3)
%     qqplot(Res)
%     title('Q-Q plot')
%     xlabel('x')
%     ylabel('y')
% 
%     %Empirical Variance
%     sr2 = 1/(n-2)*sum((y-(Beta+Alpha*x)).^2);
%     sr = sqrt(sr2);
% 
%     t0 = tinv(0.975,n-2);
% 
%     % 95% konfidensinteval 
%     % 95% sandsynlighed for at hældningen ligger imellem alpha_minus og alpha_plus
%     Alpha_minus = Alpha - t0*sr*sqrt(1/sum((x-mean(x)).^2))
%     Alpha_plus = Alpha + t0*sr*sqrt(1/sum((x-mean(x)).^2))
%     
% %% Opgave 7 - PoissonFordeling (KIG I Testkataloger.pdf)
% 
% %Poisson beskriver en spontan hændelse. Det giver sandsynligheden for en k
% %spantan hændelser pr. en tidsenhed(t).
% 
% %Eksempel - En central databaseserver modtager i gennemsnit 25 forespørgsler 
% %per sekund fra dens klienter. 
% %Serveren modtager i gennemsnittet 0,25 forespørgsler per 10 millisekunder.
% 
% %1) Sandsynlighed for at serveren ikke modtager nogen forespørgsler inden
% %for 10 millisekunders interval?
% 
% %X = antal forespørgsler inden for 10 millisekunder = 0.25
% 
% Pr_X0 = poisspdf(0,0.25) % 0.7788 dvs. 77,9% sandsynlighed for at ingen modtagelse inden for 10 ms
% 
% %2) databaseserveren modtager 89500 forespørgsler på en time 
% %   = 24,8611 forespørgsler per sekund. Stemmer dette overens med 25 per s?
% 
% %x er de observerede data 89500
% x = 89500
% 
% %tidsperioden for observeret er t = 60 * 60 = 3600 sekunder
% t = 60*60
% 
% %Nul hypotese opsættes til lambda = 25
% lambda = 25
% 
% %Z-test udføres
% Z = (x-lambda*t)/sqrt(lambda*t) %1.6667
% 
% %P-værdi
% %P = 2*(1-Phi(|Z|)) = 2*(1-normcdf(1.6667)) 
% 
% P = 2*(1-normcdf(abs(Z))) %0.0956
% 
% %Da p-værdien er over signifikansniveauet på 0.05 kan nulhypotesen ikke
% %afvises. Altså stemmer observationen overens med antagelsen at serveren i
% %gennemsnit modtager 25 forespørgsler per sekund.
% 
% %% Opgave 8 - Binomialfordeling (KIG I Testkataloger.pdf)
% 
% %1) Eksempel F13 - dartspillere
% %Sandsynlighed for at spiller A rammer sit mål mindre
% %end 5 gange ud af 10 uafhængige kast.
% 
% k = 4; %k = rammer sit mål mindre end < 5 gange så derfor fra 4
% n = 10; %n = det er ud af 10 uafhængige dart kast. (Antalsparamteren)
% P_A = 3/4; %P_A = sandsynlighed for spiller A rammer sit mål. (Sandsynligheds parameter)
% 
% % Pr(X < 5) = Pr(X <= 4) = sum(binomial(k,n,p)) = sum( n k )*P^k(1-P)^n-k
% 
% binocdf(k,n,P_A) %Binomialfordeling
% 
% %% 2) eksempel F13 - dartspillere
% %Beregn et estimat af sandsynligheden P_B for, at spiller B rammer sit mål i et givet kast.
% 
% %Spiller B ramte sit mål 41 ud af 60 uafhængige kast. 
% x = 41; %Observerede antal successer
% n = 60; %Antalsparamteren
% 
% %Paramter estimat:
% %P_hat_B = x/n = 41/60 = 0.6833 
% 
% P_B = x/n %Sandsynligheden er 0.6833 for at han rammer
% 
% 
% %3) Eksempel videre - Kan man ud fra ovenstående oplysninger konkludere, 
% %at spiller A er bedre end spiller B.
% 
% %Nulhypotesen for testen er at de er ens altså at de er lige gode
% 
% %H: P_A = P_B = 0.75
% 
% 
% %z-test
% Z = (x-n*P_A)/sqrt(n*P_A*(1-P_A)) % Hvor Z~N(0,1)
% 
% %Approksimeret p-værdi
% %P = 2*(1-Phi(abs(Z)))= 2*(1-Phi(1.1926)) %0.233
% 
% P1 = 2*(1-normcdf(abs(Z),0,1)) %0.233
% 
% %Da P er over et signifikansniveau på 0.05 er der ikke forskel på P_A og P_B,
% %, derved kan vi ikke sige at spiller A er bedre end B. 
% 
% 
% %Dette kan også regnes ved direkte brug af binomialfordeling 
% %(Mere præcis men ikke så sjov at skrive):
% 
% P2 = 2*binocdf(41,60,0.75)
% 
% 
% 
