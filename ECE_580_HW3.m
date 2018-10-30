%calculate Zr in secondary
CTR = 800/5;
VTR = (500000/sqrt(3))/67;
ZTR = VTR/CTR;

Zl1 = 0.073+0.8*1i;
Zl0 = 0.1+2.6*1i;

Zrs=(0.8*Zl1*100)/ZTR;

%% draw MHO
CxPS=0.2169/2;
CyPS=2.3767/2;
% plot(CxPS,CyPS,'g*')
% hold on


ZrS=(sqrt((0.2169^2)+(2.3767^2)))/2;

ang=0:0.01:2*pi; 
xp=ZrS*cos(ang);
yp=ZrS*sin(ang);
plot(CxPS+xp,CyPS+yp,'b');
grid on
hold on

%%
%traditional fault calculation
E=500000/sqrt(3);
Zs1=2+30*1i;
Zs0=2+30*1i;
L=20;


I1=E/(Zs1+Zl1*L+Zs1+Zl1*L+Zs0+Zl0*L);
I2=I1;
I0=I1;
V1=E-Zs1*I1;
V2=-Zs1*I2;
V0=-Zs0*I0;
kk=180/pi;
a=1*exp(1i*120/kk);

Ia=(I0+I1+I2)/CTR;
Ib=(I0+(a^2)*I1+a*I2)/CTR;
Ic=(I0+a*I1+(a^2)*I2)/CTR;
Va=(V0+V1+V2)/VTR;
Vb=(V0+(a^2)*V1+a*V2)/VTR;
Vc=(V0+a*V1+(a^2)*V2)/VTR;
Ires = Ia+Ib+Ic;
k0 = (1/3)*(Zl0-Zl1)/Zl1;

Zab = (Va-Vb)/(Ia-Ib);
Zbc = (Vb-Vc)/(Ib-Ic);
Zca = (Vc-Va)/(Ic-Ia);
Zag = (Va)/(Ia + k0*Ires);
Zbg = (Vb)/(Ib + k0*Ires);
Zcg = (Vc)/(Ic + k0*Ires);

%plot all Z    bc is infinity
plot(-1.7368,2.2338,'r*',1.9938,1.9996,'b*',0.0542,0.5942,'g*',2.2783,-1.7765,'c*',-2.6777,-1.0848,'y*')
hold on


%%
%Problem 2

Zl1 = 0.073+0.8*1i;
Zl0 = 0.1+2.6*1i;

va=30.8557+9.2392*1i;
vb=-15.5917-64.2484*1i;
vc=-49.0146+46.5991*1i;
ia=10.7296-28.8735*1i;
ib=0.2039-0.0428*1i;
ic=-0.0889-0.1311*1i;

ires=ia+ib+ic;
k0 = (1/3)*(Zl0-Zl1)/Zl1;

zab = (va-vb)/(ia-ib);
zbc = (vb-vc)/(ib-ic);
zca = (vc-va)/(ic-ia);
zag = (va)/(ia + k0*ires);
zbg = (vb)/(ib + k0*ires);
zcg = (vc)/(ic + k0*ires);

Va_all=[0.00000000000000 + 0.00000000000000i,7.82365554306144 + 0.00000000000000i,13.7956046812324 - 2.47366232683272i,16.5893632192514 - 5.26742086485177i,16.9080884235869 - 6.03689157582865i,16.9080884235869 - 3.62478242594663i,18.9324356936198 + 1.26242420831975i,24.1383638873155 + 6.46835240201546i,31.8193794295164 + 9.64993321239360i,39.8191189558014 + 9.64993321239360i,45.7945170827380 + 7.17484226763765i,48.5883367398272 + 4.38102261054845i,48.9070649176578 + 3.61154472091948i,48.9070649176578 + 6.02365079196723i,50.9314074509670 + 10.9108459907710i,56.1373295484390 + 16.1167680882430i,63.8183434099435 + 19.2983482024540i,61.0614200508903 + 19.2983482024539i,56.8452987215255 + 21.0447228376873i,55.3434232138831 + 22.5465983453297i,55.6546003686646 + 21.7953502379557i,55.6546003686646 + 21.2045133863226i,54.8150050639512 + 19.1775510147787i,51.1791888355988 + 15.5417347864265i,47.5385557428896 + 14.0337351838019i,44.1659263335182 + 14.0337351838019i,40.0145288641551 + 15.7533003184135i,38.5003714059097 + 17.2674577766588i,38.6833566953571 + 16.8256922091602i,38.6833566953570 + 16.3784629574666i,37.7944550701462 + 14.2324645982671i,34.3553667870502 + 10.7933763151711i,30.8179411581955 + 9.32812664381331i,30.2178319385425 + 9.32812664381323i,30.4712584305084 + 9.22315395377635i,30.5686833723708 + 9.12572901191397i,30.3704008222482 + 9.60442543360189i,30.3704008222482 + 9.74770835083558i,30.2676519931630 + 9.49965073414037i,30.6263470093087 + 9.85834575028599i,30.7843144534441 + 9.92377800806035i,30.4877401776908 + 9.92377800806033i,30.8382278396965 + 9.77860126501316i,31.0592831610373 + 9.55754594367233i,30.9550487396331 + 9.80919009749253i,30.9550487396330 + 9.69039638741011i,30.8316320745828 + 9.39244220082308i,31.0867534237778 + 9.64756355001800i,31.1119910181617 + 9.65801730389339i,30.8105812822848 + 9.65801730389332i,31.0656707234990 + 9.55235579772425i,31.2147165729836 + 9.40330994823969i,31.1322609074407 + 9.60237553428785i,31.1322609074407 + 9.37685496125856i,31.0054756999550 + 9.07076839383823i,31.1648397128646 + 9.23013240674790i,31.0695667570656 + 9.19066905632863i,30.7422223401466 + 9.19066905632866i,30.8785691276017 + 9.13419236777873i,30.8974039258020 + 9.11535756957841i,30.7802856219083 + 9.39810616724076i,30.7802856219083 + 9.25882303196098i,30.7063485825656 + 9.08032322881809i,30.8657095204192 + 9.23968416667174i,30.7720195349375 + 9.20087650402645i,30.5188415292211 + 9.20087650402649i,30.6436743401185 + 9.14916906072367i,30.6169159871439 + 9.17592741369835i,30.4933771179767 + 9.47417662712196i,30.4933771179767 + 9.41476082056891i,30.4776875250595 + 9.37688279256024i,30.6741006024426 + 9.57329586994331i,30.6479578119885 + 9.56246717157875i,30.5265985567152 + 9.56246717157879i,30.7208969393141 + 9.48198614635920i,30.7341613746601 + 9.46872171101324i,30.6560471779818 + 9.65730606404802i,30.6560471779818 + 9.54978076235405i,30.6480052052032 + 9.53036572260401i,30.8196137305716 + 9.70197424797228i,30.7977495015682 + 9.69291778778830i,30.7388478555470 + 9.69291778778834i,30.9525374242535 + 9.60440467029253i,30.9891062073091 + 9.56783588723691i,30.9554526618233 + 9.64908273317059i,30.9554526618233 + 9.46292876603659i,30.9266442365409 + 9.39337907500924i,31.0220301337522 + 9.48876497222028i,30.9339053857521 + 9.45226250641798i,30.8474424416758 + 9.45226250641801i,30.9917759035282 + 9.39247762901472i,30.9853873671579 + 9.39886616538492i,30.9545533965516 + 9.47330595540464i,30.9545533965515 + 9.30165129215590i,30.9283354237072 + 9.23835550653734i,30.9811874038491 + 9.29120748667897i,30.8556987799507 + 9.23922839673674i];
Vb_all=[0.00000000000000 + 0.00000000000000i,-1.80761777901446 - 1.11022302462516e-16i,-0.560872758451288 - 0.516418696338389i,2.55070108238349 - 3.62799253717317i,5.14666168847453 - 9.89519583978443i,5.14666168847454 - 18.0292216903130i,1.99099525868329 - 25.6476743834405i,-3.03148840789735 - 30.6701580500211i,-7.53838767305390 - 32.5369768498983i,-9.44932647017420 - 32.5369768498983i,-8.20460504243327 - 33.0525573466450i,-5.09306292930009 - 36.1640994597782i,-2.49709955703566 - 42.4313094405228i,-2.49709955703563 - 50.5653335944566i,-5.65276026561858 - 58.1837724753654i,-10.6752380992023 - 63.2062503089491i,-15.1821362574403 - 65.0730686503256i,-14.5484240855365 - 65.0730686503256i,-14.4187044570899 - 65.1268002797342i,-14.4737293466543 - 65.0717753901698i,-14.5066643016319 - 64.9922633751867i,-14.5066643016319 - 64.4825555648478i,-14.3266719141820 - 64.0480155019423i,-14.3968145846855 - 64.1181581724459i,-14.3885755608354 - 64.1147454570263i,-14.1488253783548 - 64.1147454570264i,-13.9310728361577 - 64.2049415132456i,-13.5627718417032 - 64.5732425077002i,-13.3674688795690 - 65.0447455676563i,-13.3674688795689 - 65.0124306164497i,-13.3192973277426 - 64.8961342027100i,-13.3968075768606 - 64.9736444518282i,-13.3459113966974 - 64.9525625637315i,-14.1208038874655 - 64.9525625637315i,-14.4795939872623 - 64.8039468383505i,-14.3100535501293 - 64.9734872754834i,-14.0955059044034 - 65.4914511115700i,-14.0955059044035 - 66.2237970541041i,-14.3520939446384 - 66.8432553807822i,-14.4774222283080 - 66.9685836644517i,-14.4866601807627 - 66.9724101496470i,-14.8105670047774 - 66.9724101496470i,-15.4009798062379 - 66.7278531598834i,-15.9418228556833 - 66.1870101104381i,-16.1546227494203 - 65.6732657209065i,-16.1546227494202 - 65.5985837334688i,-16.2387788427177 - 65.8017545152641i,-16.3191513217221 - 65.8821269942683i,-16.3301669138604 - 65.8866898019295i,-16.2663368567396 - 65.8866898019295i,-16.2287261024583 - 65.9022686864439i,-16.4501312644223 - 65.6808635244798i,-16.7029993496200 - 65.0703859637046i,-16.7029993496199 - 64.4841307839884i,-16.5623487262165 - 64.1445701414113i,-16.3708856708765 - 63.9531070860714i,-16.2235033969616 - 63.8920593493626i,-16.0004029687095 - 63.8920593493626i,-15.6003097455402 - 64.0577833886129i,-15.3221270988912 - 64.3359660352618i,-15.2751215226849 - 64.4494475348461i,-15.2751215226848 - 64.3824418433985i,-15.1965074887827 - 64.1926507765587i,-15.0095452932459 - 64.0056885810220i,-14.8950609500764 - 63.9582676134016i,-14.8149463013852 - 63.9582676134016i,-14.5860851808847 - 64.0530649934128i,-14.3188293020958 - 64.3203208722017i,-14.1568413441876 - 64.7113943971245i,-14.1568413441878 - 65.1562060987094i,-14.2696098700790 - 65.4284534033251i,-14.3457691790339 - 65.5046127122802i,-14.4611718291871 - 65.5524140551072i,-14.6349743366282 - 65.5524140551071i,-14.7566236929659 - 65.5020252418581i,-14.8235927439598 - 65.4350561908641i,-14.8182716154570 - 65.4479025314631i,-14.8182716154569 - 65.6664559799060i,-14.9235333017552 - 65.9205801705658i,-15.0736465535647 - 66.0706934223751i,-15.2717241799037 - 66.1527398616075i,-15.4877657003747 - 66.1527398616075i,-15.6874923886280 - 66.0700103585651i,-15.9045970166368 - 65.8529057305562i,-16.0338971281946 - 65.5407476476180i,-16.0338971281945 - 65.3205070706971i,-15.9836736835230 - 65.1992569494224i,-15.9276986387833 - 65.1432819046822i,-15.9259419103903 - 65.1425542439565i,-15.9274222381205 - 65.1425542439565i,-15.9089681047308 - 65.1501981962883i,-15.9613764496639 - 65.0977898513554i,-16.0297559265683 - 64.9327071908248i,-16.0297559265682 - 64.7245239821541i,-15.9343338038533 - 64.4941545993458i,-15.7571540026218 - 64.3169747981139i,-15.5916791512910 - 64.2484328704610i];
Vc_all=[0.00000000000000 + 0.00000000000000i,-6.01603778573793 - 4.44089209850063e-16i,-13.2347319490293 + 2.99008102505879i,-19.1400643206867 + 8.89541339671622i,-22.0547501377001 + 15.9320874262064i,-22.0547501377001 + 21.6540040874678i,-20.9234309680246 + 24.3852501702711i,-21.1068755024578 + 24.2018056358378i,-24.2809917832768 + 22.8870436237734i,-30.3697924547093 + 22.8870436237734i,-37.5899119460958 + 25.8777150390600i,-43.4952737181388 + 31.7830768111030i,-46.4099652627898 + 38.8197646683336i,-46.4099652627898 + 44.5416826566572i,-45.2786470865080 + 47.2729263411959i,-45.4620913216143 + 47.0894821060896i,-48.6362070083914 + 45.7747203400853i,-47.9719381279216 + 45.7747203400853i,-47.8416303401423 + 45.7207450871043i,-47.8966533044031 + 45.7757680513650i,-47.9295879859634 + 45.8552794062604i,-47.9295879859634 + 46.3649769837648i,-47.7496025344275 + 46.7995003018926i,-47.8197457960732 + 46.7293570402469i,-47.8115010476763 + 46.7327721268514i,-47.5717495538532 + 46.7327721268514i,-47.3539976826387 + 46.6425763485621i,-46.9856883451414 + 46.2742670110649i,-46.7903763920847 + 45.8027422451019i,-46.7903763920848 + 45.8350569414664i,-46.7421980828763 + 45.9513696689697i,-46.8197017240024 + 45.8738660278436i,-46.7688005463965 + 45.8949499859487i,-47.5436859994331 + 45.8949499859488i,-47.9024658428021 + 46.0435614629782i,-47.7329227639848 + 45.8740183841610i,-47.5183833325906 + 45.3560743792251i,-47.5183833325905 + 44.6237361259774i,-47.7749718276688 + 44.0042767012105i,-47.9003018106591 + 43.8789467182203i,-47.9095495499074 + 43.8751161792023i,-48.2334641650611 + 43.8751161792023i,-48.8238911756114 + 44.1196790545636i,-49.3647421573329 + 44.6605300362852i,-49.5775420518344 + 45.1742744276622i,-49.5775420518345 + 45.2489454981784i,-49.6617046366569 + 45.0457590444560i,-49.7420761742819 + 44.9653875068309i,-49.7530880547135 + 44.9608262366090i,-49.6892594574565 + 44.9608262366091i,-49.6516423213649 + 44.9452447086623i,-49.8730417644457 + 45.1666441517429i,-50.1259011800399 + 45.7771007822446i,-50.1259011800400 + 46.3633599132637i,-49.9852451580018 + 46.7029335892976i,-49.7937816202572 + 46.8943971270422i,-49.6463915738800 + 46.9554480832104i,-49.4232861283457 + 46.9554480832103i,-49.0231781004795 + 46.7897179116539i,-48.7449946033490 + 46.5115344145233i,-48.6979974860144 + 46.3980733364614i,-48.6979974860144 + 46.4650861475409i,-48.6193831899264 + 46.6548778473530i,-48.4324366622113 + 46.8418243750681i,-48.3179530454980 + 46.8892450417804i,-48.2378563131894 + 46.8892450417803i,-48.0090072443704 + 46.7944526537390i,-47.7417611121924 + 46.5272065215609i,-47.5797726614267 + 46.1361318067743i,-47.5797726614266 + 45.6913043219579i,-47.6925428372029 + 45.4190530341677i,-47.7687034744037 + 45.3428923969670i,-47.8840902802211 + 45.2950976170786i,-48.0578946022245 + 45.2950976170786i,-48.1795475324387 + 45.3454879106756i,-48.2465087179576 + 45.4124490961944i,-48.2411784017486 + 45.3995805745113i,-48.2411784017487 + 45.1810310445012i,-48.3464380108578 + 44.9269118686200i,-48.4965340055727 + 44.7768158739051i,-48.6946290500294 + 44.6947622198522i,-48.9106463898099 + 44.6947622198522i,-49.1103730879323 + 44.7774917269824i,-49.3274702162860 + 44.9945888553359i,-49.4567796298990 + 45.3067693954238i,-49.4567796298990 + 45.5270268569657i,-49.4065608272253 + 45.6482657714669i,-49.3505761613613 + 45.7042504373308i,-49.3488390325974 + 45.7049699796244i,-49.3503050453200 + 45.7049699796244i,-49.3318583082383 + 45.6973290909431i,-49.3842777142666 + 45.7497484969719i,-49.4526573351445 + 45.9148315050850i,-49.4526573351446 + 46.1229980463258i,-49.3572368226137 + 46.3533635418064i,-49.1800690051075 + 46.5305313593126i,-49.0145722055566 + 46.5990823782159i];
Ia_all=[0.00000000000000 + 0.00000000000000i,-0.00641574954564015 - 4.33680868994202e-19i,-0.0193227456571296 + 0.00534625283887574i,-0.0330710703010380 + 0.0190945774827842i,-0.0414735245218497 + 0.0393798964198866i,-0.0414735245218496 + 0.0605075343212814i,-0.0349364750091461 + 0.0762893679127547i,-0.0275571047533852 + 0.0836687381685156i,-0.0255235362721595 + 0.0845110698134538i,-0.0318923784192621 + 0.0845110698134538i,-0.0447982429087738 + 0.0898568539191588i,-0.0585464614768438 + 0.103605072487229i,-0.0669489098318480 + 0.123890377263020i,-0.0669489098318481 + 0.145018049227287i,-0.0604118266774166 + 0.160799964037076i,-0.0530323308367773 + 0.168179459877715i,-0.0509985426577385 + 0.169021882524467i,1.03253750516412 + 0.169021882524468i,3.22348429678567 - 0.738497993002994i,5.56623674451050 - 3.08125044072782i,6.95494066407782 - 6.43387827746791i,6.95494066407783 - 9.81421598898842i,5.94302771059266 - 12.2571899652333i,4.95358265119679 - 13.2466350246292i,5.04176647491411 - 13.2101080888636i,6.57955577729558 - 13.2101080888636i,9.19181379560531 - 14.2921407884653i,11.8415183370175 - 16.9418453298775i,13.3927052647450 - 20.6867418485730i,13.3927052647450 - 24.4878551150746i,12.2275607993221 - 27.3007626856225i,10.9626458511775 - 28.5656776337671i,10.6925931319315 - 28.6775371326346i,10.7750595120662 - 28.6775371326346i,10.8658136819128 - 28.7151287406270i,10.9201139727181 - 28.7694290314322i,10.9506567024376 - 28.8431657037527i,10.9506567024375 - 28.9333750130063i,10.9239932974936 - 28.9977461668411i,10.8766589545865 - 29.0450805097481i,10.8032837319462 - 29.0754735221079i,10.7382468818224 - 29.0754735221079i,10.6871342067775 - 29.0543019588951i,10.6383970963645 - 29.0055648484821i,10.6159748163260 - 28.9514326759139i,10.6159748163261 - 28.9052072839556i,10.6386472835369 - 28.8504711061229i,10.6783158289722 - 28.8108025606875i,10.7175178245332 - 28.7945645624540i,10.7698577661071 - 28.7945645624541i,10.8207766661024 - 28.8156558614132i,10.8494962001759 - 28.8443753954866i,10.8680255520820 - 28.8891092081605i,10.8680255520821 - 28.9447396022105i,10.8514813493017 - 28.9846808409415i,10.8189037908840 - 29.0172583993592i,10.7673995356933 - 29.0385921603791i,10.7223899759639 - 29.0385921603791i,10.6822930721459 - 29.0219834790085i,10.6437247413574 - 28.9834151482200i,10.6259897361057 - 28.9405990580124i,10.6259897361057 - 28.8998925280833i,10.6452016247734 - 28.8535109259029i,10.6769992478284 - 28.8217133028479i,10.7100573657417 - 28.8080201820617i,10.7538692392587 - 28.8080201820617i,10.7926202558700 - 28.8240713786978i,10.8145302841287 - 28.8459814069565i,10.8284061233704 - 28.8794806462430i,10.8284061233703 - 28.9166333943622i,10.8184021768216 - 28.9407850577974i,10.7980230153174 - 28.9611642193017i,10.7674637559564 - 28.9738222789851i,10.7440913277581 - 28.9738222789850i,10.7223250429339 - 28.9648063886084i,10.7012692117501 - 28.9437505574244i,10.6927249689385 - 28.9231229305489i,10.6927249689386 - 28.9026094197607i,10.7034700987152 - 28.8766683817246i,10.7202534678534 - 28.8598850125863i,10.7386442403442 - 28.8522673051981i,10.7660853716868 - 28.8522673051981i,10.7897368257911 - 28.8620640582579i,10.8041560162670 - 28.8764832487338i,10.8144248652421 - 28.9012744432000i,10.8144248652422 - 28.9287551535105i,10.8063101276494 - 28.9483458630622i,10.7883145698320 - 28.9663414208795i,10.7627438816769 - 28.9769331467125i,10.7418604481407 - 28.9769331467125i,10.7208611601848 - 28.9682349568408i,10.7024064758066 - 28.9497802724628i,10.6948747583297 - 28.9315970979821i,10.6948747583298 - 28.9126327948716i,10.7036242736454 - 28.8915095963326i,10.7159440527676 - 28.8791898172103i,10.7295715721041 - 28.8735451138796i];
Ib_all=[0.00000000000000 + 0.00000000000000i,0.0211015757815669 + 8.67361737988404e-19i,0.0412143751774269 - 0.00833099428705460i,0.0544792036612935 - 0.0215958227709212i,0.0594100928240186 - 0.0335000422621306i,0.0594100928240186 - 0.0385490361952385i,0.0607708068885076 - 0.0352639818462373i,0.0689867465343685 - 0.0270480422003765i,0.0855367136246613 - 0.0201928213747488i,0.107017548466545 - 0.0201928213747488i,0.127137668892407 - 0.0285268481317206i,0.140402584092656 - 0.0417917633319699i,0.145333478555152 - 0.0536959956179566i,0.145333478555152 - 0.0587450187942908i,0.146694163763796 - 0.0554600341094619i,0.154910025933212 - 0.0472441719400459i,0.171459878633816 - 0.0403889984961787i,0.182427993620287 - 0.0403889984961787i,0.198095067650125 - 0.0468785130420406i,0.219390582832897 - 0.0681740282248123i,0.227423988205362 - 0.0875683844270588i,0.227423988205362 - 0.0985511674564492i,0.222287108612203 - 0.110952691838530i,0.217683432394817 - 0.115556368055917i,0.216876239247637 - 0.115890718404934i,0.209881597815374 - 0.115890718404934i,0.208487025561674 - 0.115313067663742i,0.216860155734897 - 0.123686197836965i,0.221310026777037 - 0.134429136857711i,0.221310026777037 - 0.148042522669128i,0.213742083490207 - 0.166313153991463i,0.204825308752422 - 0.175229928729248i,0.197467231500324 - 0.178277744120056i,0.171226747119477 - 0.178277744120057i,0.142341073114494 - 0.166312906188905i,0.115477928046196 - 0.139449761120606i,0.103445564832477 - 0.110401066662648i,0.103445564832477 - 0.0936136466695635i,0.107068224088036 - 0.0848677735629366i,0.110909564667238 - 0.0810264329837351i,0.112920385562110 - 0.0801935236975756i,0.114668819617398 - 0.0801935236975757i,0.111424530210281 - 0.0788496950248841i,0.100015779727571 - 0.0674409445421750i,0.0908643955605032 - 0.0453475487715534i,0.0908643955605032 - 0.0190906188114960i,0.101206097353082 + 0.00587645791416660i,0.117441986625703 + 0.0221123471867877i,0.135451149679672 + 0.0295719867707297i,0.158176037226851 + 0.0295719867707293i,0.180546092610869 + 0.0203060064396317i,0.194071246872589 + 0.00678085217791200i,0.199390730599015 - 0.00606151757924932i,0.199390730599015 - 0.0125572917660231i,0.200488530101048 - 0.00990696931944912i,0.204843858865853 - 0.00555164055464422i,0.210531271137220 - 0.00319583725703645i,0.221034826472719 - 0.00319583725703612i,0.235220349485957 - 0.00907167327847545i,0.248416199715123 - 0.0222675235076408i,0.257421194911700 - 0.0440075050403216i,0.257421194911701 - 0.0684584607278146i,0.249845527655848 - 0.0867477393609192i,0.237913316990560 - 0.0986799500262069i,0.223416292183735 - 0.104684814315254i,0.210867215378297 - 0.104684814315254i,0.201108152031291 - 0.100642477920866i,0.195365775969390 - 0.0949001018589643i,0.194998349825387 - 0.0940130566789423i,0.194998349825387 - 0.0996104726231898i,0.191738339201814 - 0.107480834484099i,0.183448349288281 - 0.115770824397632i,0.170507153708135 - 0.121131243120252i,0.156556315106896 - 0.121131243120252i,0.141591169641251 - 0.114932476905495i,0.128019258176450 - 0.101360565440694i,0.121070237414086 - 0.0845841452709840i,0.121070237414086 - 0.0686512297232371i,0.126370212425646 - 0.0558559581700913i,0.133063998955324 - 0.0491621716404132i,0.138412017289871 - 0.0469469499144237i,0.142995578957348 - 0.0469469499144235i,0.144283203526789 - 0.0474803014743303i,0.141990132221902 - 0.0451872301694435i,0.139565371465025 - 0.0393333398646824i,0.139565371465025 - 0.0294014373894782i,0.144684444345235 - 0.0170429022153002i,0.154262436972270 - 0.00746490958826605i,0.166498734064405 - 0.00239646937947682i,0.180869482859142 - 0.00239646937947661i,0.193718120175906 - 0.00771854921409165i,0.202646479492060 - 0.0166469085302467i,0.207220041903127 - 0.0276884649314046i,0.207220041903127 - 0.0367932360612885i,0.205372183828264 - 0.0412543600869629i,0.203728735232426 - 0.0428978086828028i,0.203863292918901 - 0.0428420730641429i];
Ic_all=[0.00000000000000 + 0.00000000000000i,-0.0146858261915187 - 8.67361737988404e-19i,-0.0218916294436777 + 0.00298474143483638i,-0.0214081331947695 + 0.00250124518592817i,-0.0179365681219932 - 0.00587985429542902i,-0.0179365681219932 - 0.0219584982827639i,-0.0258343316961640 - 0.0410253862159433i,-0.0414296416884232 - 0.0566206962082025i,-0.0600131772896028 - 0.0643182486910544i,-0.0751251698632281 - 0.0643182486910544i,-0.0823394258836199 - 0.0613300060049764i,-0.0818561223869979 - 0.0618133095015983i,-0.0783845684962644 - 0.0701943819871163i,-0.0783845684962644 - 0.0862730306514861i,-0.0862823368391996 - 0.105339930097481i,-0.101877694871900 - 0.120935288130182i,-0.120461335700569 - 0.128632884199687i,-0.110298675689743 - 0.128632884199687i,-0.0946474334893596 - 0.135115840987072i,-0.0733522068459087 - 0.156411067630522i,-0.0653188231742539 - 0.175805371442378i,-0.0653188231742539 - 0.186788086855827i,-0.0704556705086473 - 0.199189533358359i,-0.0750593213031156 - 0.203793184152828i,-0.0758665486437144 - 0.204127548665222i,-0.0828612178352498 - 0.204127548665222i,-0.0842557999285440 - 0.203549893848337i,-0.0758826837515049 - 0.211923010025376i,-0.0714328287688399 - 0.222665910275119i,-0.0714328287688396 - 0.236279251660166i,-0.0790007699372923 - 0.254549877868286i,-0.0879175336377567 - 0.263466641568751i,-0.0952755866001886 - 0.266514446898450i,-0.121516034551015 - 0.266514446898450i,-0.150401662218435 - 0.254549628160945i,-0.177264831832190 - 0.227686458547190i,-0.189297217470667 - 0.198637709951076i,-0.189297217470667 - 0.181850286860036i,-0.185674544845112 - 0.173104381475385i,-0.181833157747436 - 0.169262994377708i,-0.179822258813049 - 0.168430052766523i,-0.178073803400423 - 0.168430052766523i,-0.181318131712300 - 0.167086207978953i,-0.192726894845104 - 0.155677444846149i,-0.201878282562048 - 0.133584040505369i,-0.201878282562048 - 0.107327067449864i,-0.191536552052810 - 0.0823599213960561i,-0.175300638403265 - 0.0661240077465098i,-0.157291505702433 - 0.0586643807352480i,-0.134566630267998 - 0.0586643807352476i,-0.112196581887257 - 0.0679303581654930i,-0.0986714126662810 - 0.0814555273864695i,-0.0933519274323205 - 0.0942979007831407i,-0.0933519274323203 - 0.100793612822969i,-0.0922541356191703 - 0.0981433089389992i,-0.0878988096906028 - 0.0937879830104319i,-0.0822113719745292 - 0.0914321691732822i,-0.0717077696251537 - 0.0914321691732822i,-0.0575221967647180 - 0.0973080258421064i,-0.0443263493546167 - 0.110503873252208i,-0.0353213617610920 - 0.132243836429496i,-0.0353213617610917 - 0.156694812862636i,-0.0428970189078043 - 0.174984067090118i,-0.0548292095714665 - 0.186916257753780i,-0.0693261566089831 - 0.192921089829723i,-0.0818751999128564 - 0.192921089829723i,-0.0916342974781885 - 0.188878739261640i,-0.0973766669856800 - 0.183136369754149i,-0.0977441035424347 - 0.182249299435519i,-0.0977441035424351 - 0.187846732850757i,-0.101004075311285 - 0.195717000908069i,-0.109294051638943 - 0.204006977235727i,-0.122235278889164 - 0.209367409076520i,-0.136186128996641 - 0.209367409076520i,-0.151151254063535 - 0.203168651311203i,-0.164723153416516 - 0.189596751958222i,-0.171672180775989 - 0.172820315861681i,-0.171672180775989 - 0.156887353939948i,-0.166372229564567 - 0.144092139845419i,-0.159678429825999 - 0.137398340106849i,-0.154330398243791 - 0.135183112893499i,-0.149746806559756 - 0.135183112893500i,-0.148459096970754 - 0.135716499669662i,-0.150752198396104 - 0.133423398244311i,-0.153176941813507 - 0.127569549800745i,-0.153176941813506 - 0.117637595633432i,-0.148057879655070 - 0.105279086343904i,-0.138479838568846 - 0.0957010452576783i,-0.126243475361849 - 0.0906325776632175i,-0.111872721892482 - 0.0906325776632176i,-0.0990241278469595 - 0.0959546395742999i,-0.0900957521101449 - 0.104883015311114i,-0.0855221829521363 - 0.115924588000830i,-0.0855221829521361 - 0.125029396057036i,-0.0873700065451421 - 0.129490436836145i,-0.0890134637796577 - 0.131133894070660i,-0.0888789433509627 - 0.131078173884678i];

Ires_all=[0.00000000000000 + 0.00000000000000i,-4.46005211929279e-12 + 4.33680868994202e-19i,-1.48227576046711e-11 + 4.29237635019009e-12i,5.20801943237892e-11 - 6.26105712414615e-11i,6.25179942537191e-11 - 8.78096918094373e-11i,6.25180254787416e-11 - 1.17393907095309e-10i,5.38375455327866e-11 - 1.38350379763619e-10i,-6.93034102328127e-11 - 2.61491328590324e-10i,-1.49183415620868e-10 - 2.94578694770564e-10i,-8.24977169910923e-11 - 2.94578708648352e-10i,-2.09308195775471e-10 - 2.42052100585255e-10i,-1.02474834173805e-10 - 3.48885469125815e-10i,-1.08502928863885e-10 - 3.34332325713227e-10i,-1.08502914986097e-10 - 2.21302254299616e-10i,-1.00063055064048e-10 - 2.00926608684426e-10i,-1.55208026986209e-10 - 2.56071511217648e-10i,-1.54448009936914e-10 - 2.55756749112379e-10i,1.10466682253898 - 2.55755860933959e-10i,3.32693192991271 - 0.920492346918452i,5.71227511863014 - 3.30583553563587i,7.11704582693756 - 6.69725203165612i,7.11704582693757 - 10.0995552410888i,6.09485914679891 - 12.5673321875567i,5.09620676066330 - 13.5659845736923i,5.18277616417664 - 13.5301263526706i,6.70657615666199 - 13.5301263526706i,9.31604502106428 - 14.6110037468963i,11.9824958087992 - 17.2774545346312i,13.5425824625558 - 21.0438368926076i,13.5425824625558 - 24.8721768860777i,12.3623021128931 - 27.7216257136358i,11.0795536264759 - 29.0043742000530i,10.7947847771855 - 29.1223293195705i,10.8247702257628 - 29.1223293195705i,10.8577530950026 - 29.1359912713356i,10.8583270723270 - 29.1365652486600i,10.8648050537963 - 29.1522044793798i,10.8648050537962 - 29.2088389466231i,10.8453869804062 - 29.2557183227565i,10.8057353646744 - 29.2953699384884i,10.7363818612537 - 29.3240971002033i,10.6748418998385 - 29.3240971002033i,10.6172406065657 - 29.3002378633195i,10.5456859825607 - 29.2286832393145i,10.5049609308064 - 29.1303642670408i,10.5049609308064 - 29.0316249727557i,10.5483168299939 - 28.9269545729286i,10.6204571778347 - 28.8548142250877i,10.6956774685133 - 28.8236569605230i,10.7934671720856 - 28.8236569605230i,10.8891261747575 - 28.8632802167928i,10.9448960315281 - 28.9190500735633i,10.9740643519999 - 28.9894686284381i,10.9740643519999 - 29.0580905079611i,10.9597157406729 - 29.0927311200277i,10.9358488371343 - 29.1165980235664i,10.8957194321772 - 29.1332201673496i,10.8717170301932 - 29.1332201673496i,10.8599912221001 - 29.1283631786077i,10.8478145886767 - 29.1161865451843i,10.8480895659328 - 29.1168503990053i,10.8480895659329 - 29.1250458003164i,10.8521501304966 - 29.1152427302757i,10.8600833527678 - 29.1073095080045i,10.8641474995752 - 29.1056260832774i,10.8828612537328 - 29.1056260832774i,10.9020941101580 - 29.1135925932519i,10.9125193933837 - 29.1240178764775i,10.9256603701072 - 29.1557430007061i,10.9256603701072 - 29.2040905984047i,10.9091364411358 - 29.2439828918312i,10.8721773134669 - 29.2809420195001i,10.8157356314888 - 29.3043209296586i,10.7644615148504 - 29.3043209296586i,10.7127649599160 - 29.2829075154768i,10.6645653183925 - 29.2347078739532i,10.6421230277359 - 29.1805273914795i,10.6421230277359 - 29.1281480039651i,10.6634680834530 - 29.0766164809633i,10.6936390384337 - 29.0464455259826i,10.7227258603649 - 29.0343973698524i,10.7593341445345 - 29.0343973698523i,10.7855609324147 - 29.0452608610898i,10.7953939499889 - 29.0550938786639i,10.8008132947462 - 29.0681773342766i,10.8008132947462 - 29.0757941880316i,10.8029366920844 - 29.0706678533797i,10.8040971677464 - 29.0695073777175i,10.8029991395317 - 29.0699621958959i,10.8108572077523 - 29.0699621958959i,10.8155551506507 - 29.0719081475595i,10.8149572009127 - 29.0713101978217i,10.8165726147596 - 29.0752101518396i,10.8165726147596 - 29.0744554273425i,10.8216264485763 - 29.0622543932003i,10.8306593221360 - 29.0532215196405i,10.8445559198515 - 29.0474653603959i];

Zag_all = (Va_all)./(Ia_all + k0.*Ires_all);

i=1;
x=real(Zag_all(1,i:97));
y=imag(Zag_all(1,i:97));
plot(x,y,'r*')
hold on

%% All data six MHO

t=[0,0.00104166666666667,0.00208333333333333,0.00312500000000000,0.00416666666666667,0.00520833333333333,0.00625000000000000,0.00729166666666667,0.00833333333333333,0.00937500000000000,0.0104166666666667,0.0114583333333333,0.0125000000000000,0.0135416666666667,0.0145833333333333,0.0156250000000000,0.0166666666666667,0.0177083333333333,0.0187500000000000,0.0197916666666667,0.0208333333333333,0.0218750000000000,0.0229166666666667,0.0239583333333333,0.0250000000000000,0.0260416666666667,0.0270833333333333,0.0281250000000000,0.0291666666666667,0.0302083333333333,0.0312500000000000,0.0322916666666667,0.0333333333333333,0.0343750000000000,0.0354166666666667,0.0364583333333333,0.0375000000000000,0.0385416666666667,0.0395833333333333,0.0406250000000000,0.0416666666666667,0.0427083333333333,0.0437500000000000,0.0447916666666667,0.0458333333333333,0.0468750000000000,0.0479166666666667,0.0489583333333333,0.0500000000000000,0.0510416666666667,0.0520833333333333,0.0531250000000000,0.0541666666666667,0.0552083333333333,0.0562500000000000,0.0572916666666667,0.0583333333333333,0.0593750000000000,0.0604166666666667,0.0614583333333333,0.0625000000000000,0.0635416666666667,0.0645833333333333,0.0656250000000000,0.0666666666666667,0.0677083333333333,0.0687500000000000,0.0697916666666667,0.0708333333333333,0.0718750000000000,0.0729166666666667,0.0739583333333333,0.0750000000000000,0.0760416666666667,0.0770833333333333,0.0781250000000000,0.0791666666666667,0.0802083333333333,0.0812500000000000,0.0822916666666667,0.0833333333333333,0.0843750000000000,0.0854166666666667,0.0864583333333333,0.0875000000000000,0.0885416666666667,0.0895833333333333,0.0906250000000000,0.0916666666666667,0.0927083333333333,0.0937500000000000,0.0947916666666667,0.0958333333333333,0.0968750000000000,0.0979166666666667,0.0989583333333333,0.100000000000000];

MHOAB=zeros(1,97);
MHOBC=zeros(1,97);
MHOCA=zeros(1,97);
MHOAG=zeros(1,97);
MHOBG=zeros(1,97);
MHOCG=zeros(1,97);

Sab1=((Ia_all-Ib_all)*Zrs)-(Va_all-Vb_all);
Sab2=Va_all-Vb_all;
 
Sbc1=((Ib_all-Ic_all)*Zrs)-(Vb_all-Vc_all);
Sbc2=Vb_all-Vc_all;
 
Sca1=((Ic_all-Ia_all)*Zrs)-(Vc_all-Va_all);
Sca2=Vc_all-Va_all;
 
Sag1=((Ia_all+k0*Ires_all)*Zrs)-Va_all;
Sag2=Va_all;
 
Sbg1=((Ib_all+k0*Ires_all)*Zrs)-Vb_all;
Sbg2=Vb_all;
 
Scg1=((Ic_all+k0*Ires_all)*Zrs)-Vc_all;
Scg2=Vc_all;
    
i=1;
while (i <= 97)
    
    if real(Sab1(1,i)*conj(Sab2(1,i)))>0
        MHOAB(1,i)=1;
       
    elseif real(Sbc1(1,i)*conj(Sbc2(1,i)))>0
        MHOBC(1,i)=1;
        
    elseif real(Sca1(1,i)*conj(Sca2(1,i)))>0
        MHOCA(1,i)=1;
       
    elseif real(Sag1(1,i)*conj(Sag2(1,i)))>0
        MHOAG(i,1)=1;
    
    elseif real(Sbg1(1,i)*conj(Sbg2(1,i)))>0
        MHOBG(1,i)=1;
        
    elseif real(Scg1(1,i)*conj(Scg2(1,i)))>0
        MHOCG(1,i)=1;
    end
    i=i+1;
end

subplot(6,1,1),
plot(t,MHOAB);
subplot(6,1,2),
plot(t,MHOBC);
subplot(6,1,3),
plot(t,MHOCA);
subplot(6,1,4),
plot(t,MHOAG);
subplot(6,1,5),
plot(t,MHOBG);
subplot(6,1,6),
plot(t,MHOCG);

%% AG element  last 

CTR = 800/5;
VTR = (500000/sqrt(3))/67;
ZTR = VTR/CTR;
Zl1 = 0.073+0.8*1i;
Zr=(0.8*Zl1*100)/ZTR;


K0=(1/3)*(Zl0-Zl1)/Zl1;
V=Va_all(1,97);
I=Ia_all(1,97)+K0*Ires_all(1,97);
Vpol=(1)*(Vb_all(1,97)-Vc_all(1,97));

Vp=V-Vpol;
Zp=Vp/I;
% Zr=V/(I);

Zc=(Zr+Zp)/2;
r=((abs(Zr-Zp))/2);



% draw Cross MHO
CxPS=real(Zc)+0.9;
CyPS=imag(Zc)-1.71;



% ZrS=(sqrt((0.2169^2)+(2.3767^2)))/2;
ZrS=r+1.15;

ang=0:0.01:2*pi; 
xp=ZrS*cos(ang);
yp=ZrS*sin(ang);
plot(CxPS+xp,CyPS+yp,'r');
grid on
hold on



