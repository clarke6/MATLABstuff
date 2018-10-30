clear all
training = [
    30,32,35,37,40,45,47,47,44,42,40,37;
    69,71,75,79,84,88,87,85,82,77,74,70;
    70,74,77,80,84,87,90,95,94,90,85,80];
categories = {'Winter','Summer','Summer'};
MODEL=fitcsvm(training,categories);
new_day = [65,69,70,72,72,74,71,69,67,64,61,60];
disp(new_day)
predict(MODEL,new_day)
k = 1;
while k <= length(new_day),
    new_day(k) = new_day(k) - 30;
    k = k + 1;
end
disp(new_day)
predict(MODEL,new_day)