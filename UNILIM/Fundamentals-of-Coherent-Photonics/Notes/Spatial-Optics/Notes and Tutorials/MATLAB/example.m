%% Lab Two Code %%
clear all

set(0, 'defaulttextinterpreter','Latex');
%% Task 4 Checking
T4_tf = tf(1,[1 5 6]);

% Open output file to write variables for Latex
[L2_T4_Out]=fopen('Outputs/MATLAB_output_example.txt','w'); 
tf_string = evalc('T4_tf');
fprintf(L2_T4_Out,'%s',tf_string);
fclose(L2_T4_Out);

poles_T4_tf = pole(T4_tf);

% Open output file to write variables for Latex
[L2_T4_Out]=fopen('Outputs/MATLAB_output_example.txt','a'); 
tf_string = evalc('T4_tf');
fprintf(L2_T4_Out,'Poles =');
fclose(L2_T4_Out);

writematrix(poles_T4_tf,'Outputs/MATLAB_output_example.txt','WriteMode','append');

%% Task 5 Checking

T5_tf = tf(1,[2 2]);

% Open output file to write variables for Latex
[L2_T5_Out]=fopen('Outputs/MATLAB_output_example_2.txt','w'); 
tf_string = evalc('T5_tf');
fprintf(L2_T5_Out,'%s',tf_string);
fclose(L2_T5_Out);

poles_T5_tf = pole(T5_tf);

% Open output file to write variables for Latex
[L2_T5_Out]=fopen('Outputs/MATLAB_output_example_2.txt','a'); 
tf_string = evalc('T5_tf');
fprintf(L2_T5_Out,'Poles =');
fclose(L2_T5_Out);

writematrix(poles_T5_tf,'Outputs/MATLAB_output_example_2.txt','WriteMode','append');