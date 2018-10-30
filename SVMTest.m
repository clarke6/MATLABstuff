training = [6,0,2,0,20,2,0;2,1,2,2,250,2,0;0,0,0,1,5000,1,0;4,1,0,2,750,0,1;4,1,0,2,50000,0,1;4,0,0,2,500,1,0];
categories = {'Insect','Bird','Fish','Small Mammal','Large Mammal','Reptile'};
MODEL = fitcsvm(training,categories)