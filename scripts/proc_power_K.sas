proc power; 
  twosamplefreq test=pchi 
  refproportion=0.10 /*index p, since the exposure is protective, reverse this for calculation here but reverse back once get the RR*/
  power = .80
  relativerisk=1.18 to 1.23 by 0.01
  npergroup=.;
run;

*using the smallest group number from HRS dataset. 1294+150=1444;
*RR = 1.33;
*1/1.33 = 0.76;
/*Computed N per Group
Index	Relative Risk	Actual Power	N per Group
1	1.30	0.800	1774
2	1.31	0.800	1668
3	1.32	0.800	1571
4	1.33	0.800	1483
5	1.34	0.800	1402
6	1.35	0.800	1328
7	1.36	0.800	1260
8	1.37	0.800	1197
9	1.38	0.800	1139
10	1.39	0.800	1085
11	1.40	0.800	1035*/

*using the towards larger side group number; *3093; 
*RR = 1.21;
*1/1.21 = 0.83;
/*
Computed N per Group
Index	Relative Risk	Actual Power	N per Group
1	1.18	0.800	4705
2	1.19	0.800	4239
3	1.20	0.800	3841
4	1.21	0.800	3498
5	1.22	0.800	3200
6	1.23	0.800	2939
*/