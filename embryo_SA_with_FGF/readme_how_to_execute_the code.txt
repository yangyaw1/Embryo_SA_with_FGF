1: To execute the code, please copy the following into the commend line:
gcc -o exe main.c -lOpenCL -std=c99 -lm
./exe

2: You need to choose the platform to execute the code in line 474
(It is show as platform_id[0] in the code, if you have error when executing,
you can try platform_id[1], platform_id[2], ...)

3: You can set the parameter related to the results of the paper at line 20-line 44

4: Once you executed the code, you will get space files (in "PQR" folder) and gene files (in "signal" folder).


