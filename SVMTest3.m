clear all
training = [
    30,32,35,37,40,45,47,47,44,42,40,37;
    59,61,65,69,74,75,76,75,72,67,64,60;
    70,74,77,80,84,87,90,95,94,90,85,80;
    45,47,49,51,53,55,57,58,55,53,50,48];
categories = {'Winter','Spring','Summer','Fall'};
MODEL=fitcecoc(training,categories)