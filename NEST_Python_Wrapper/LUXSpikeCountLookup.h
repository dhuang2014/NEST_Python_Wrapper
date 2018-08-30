#ifndef __LUXSpikeCountLookup_H__
#define __LUXSpikeCountLookup_H__

// Probability coefficients for partitioning number of PHE into top and bottom arrays.
// Comes from Evan's Tritium analysis  
double topArrayProbabilityCoeffs[3] = {1.10677683e-07, -7.23082313e-04, 3.36092810e-01};

// Tomasz's spike calibration lookup tables
//#True_S1(phd) Spikes_S1(phd) SDEV_Spikes_S1(phd)

double true_S1[100] = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
                       22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,
                       43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,
                       63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,
                       84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101 };


double spikes_S1[100] = {
1.99679,2.96876,3.94683,4.90814,5.86752,6.80822,7.75766,8.68038,9.60107,10.5086,11.4352,
12.3059,13.2187,14.1029,14.9768,15.8248,16.7015,17.5502,18.4067,19.233,20.0594,20.881,
21.7004,22.4926,23.3636,24.1025,24.8751,25.6897,26.4774,27.1799,27.8757,28.7142,29.3949,30.2456,30.9857,
31.6475,32.3673,33.0385,33.7713,34.5531,35.2869,35.9793,36.7272,37.1953,38.0447,38.7265,39.3292,
39.9918,40.6495,41.3445,41.9674,42.7357,43.3772,43.772,44.7766,45.1056,45.6693,46.3948,46.9735,47.5282,
48.1869,48.7751,49.293,49.923,50.5167,51.1609,51.7998,52.3593,52.8243,53.6332,53.9315,54.44,55.1031,
55.7808,56.1518,56.662,57.2615,57.6953,58.3935,58.6762,59.1871,59.7136,60.2945,60.8945,61.2837,
61.698,62.2539,62.6438,63.6058,64.0361,64.029,64.9631,65.1398,65.9543,66.3628,66.8205,67.3376,67.7828,
67.7412,68.7244};

double sdev_spikes_S1[100] = {
0.070045,0.178146,0.232193,0.303243,0.358738,0.422288,0.479704,0.546591,0.609808,0.672361,0.715478,
0.784271,0.838726,0.882284,0.945144,0.997971,1.04986,1.10885,1.17028,1.22058,1.27285,1.32461,1.35563,
1.39081,1.44916,1.49548,1.55712,1.57177,1.65802,1.68188,1.78498,1.82452,1.83148,1.97509,1.96888,1.85443,
1.96372,2.0752,2.04838,2.15156,2.2021,2.21821,2.13843,2.31059,2.36143,2.27017,2.34899,2.4843,2.51542,
2.53412,2.61608,2.53476,2.61154,2.63662,2.51964,2.77855,2.85784,2.83935,2.83767,2.95355,2.80571,2.96081,
2.95098,2.91599,3.01217,2.96976,3.0199,2.98849,3.06377,3.12445,3.12768,3.09304,3.1166,3.17492,3.30526,
3.19687,3.23392,3.39693,3.51929,3.33133,3.43802,3.35899,3.47708,3.53378,3.46362,3.58777,3.55388,3.5137,
3.73149,3.59951,3.70084,3.63223,3.72674,3.68116,3.58975,3.85882,3.50127,4.14561,3.80175,3.82282};


double num_spikes_lookup[80] = {
 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,
13,14,15,16,17,18,19,20,21,22,23,24,
25,26,27,28,29,30,31,32,33,34,35,36,
37,38,39,40,41,42,43,44,45,46,47,48,
49,50,51,52,53,54,55,56,57,58,59,60,
61,62,63,64,65,66,67,68,69,70,71,72,
73,74,75,76,77,78,79,80};

double SC_offset[80] = {
9.4898e-01,8.8959e-01,9.9934e-01,1.0013e+00,9.9362e-01,9.0993e-01,
8.9202e-01,9.2374e-01,8.2519e-01,1.3029e+00,1.1644e+00,1.0565e+00,
1.0789e+00,1.1383e+00,1.2057e+00,7.5864e-01,8.0802e-01,4.1297e-01,
1.2876e+00,-2.7253e-02,1.4167e+00,5.5438e-01,2.2096e-01,2.6715e-01,
2.0470e+00,1.2796e+00,1.1017e+00,3.6675e-01,-8.3168e-02,1.0918e+00,
2.5046e-01,-4.8226e-01,7.6286e-02,2.6616e+00,1.2315e-02,1.3216e+00,
1.2644e+00,1.8070e+00,7.7750e-01,1.5503e+00,7.2700e-01,7.8651e-01,
1.3119e+00,8.6602e-01,-8.7957e-01,2.2035e-01,-5.2797e-01,1.6043e+00,
1.8140e+00,2.7407e+00,-2.2064e-01,-4.1029e-01,1.0536e+00,1.0569e+00,
8.2461e-01,2.1863e+00,2.3877e+00,-1.8133e-01,5.7445e-01,-1.2000e+00,
2.3802e+00,8.4822e-01,5.2737e-01,1.2374e+00,7.8560e-01,1.1424e+00,
1.7438e+00,4.5310e-01,1.2468e-01,1.6962e+00,3.5355e-01,2.0620e+00,
5.8929e-01,8.6273e-01,2.4041e-01,1.1105e+00,-1.7123e-01,1.7225e+00,1.2714e+00,-5.0883e-01};

double SC_lin[80] = {
1.8887e-01,1.8858e-01,3.0562e-02,1.1891e-02,1.1372e-02,4.7620e-02,5.5361e-02,2.8823e-02,5.8942e-02,-8.7477e-02,-3.8534e-02,-1.1232e-02,
-1.2889e-02,-2.5960e-02,-3.3841e-02,4.1433e-02,3.7722e-02,7.9976e-02,-4.7081e-02,1.3231e-01,-5.8732e-02,5.1039e-02,8.0301e-02,7.8182e-02,
-1.1205e-01,-2.4931e-02,-9.4881e-03,5.2296e-02,8.2768e-02,-1.5298e-02,5.2305e-02,1.0859e-01,6.2261e-02,-1.2391e-01,6.4751e-02,-2.3787e-02,
-2.2201e-02,-5.9814e-02,6.6715e-03,-3.7362e-02,1.0740e-02,7.5974e-03,-2.2918e-02,3.6668e-03,8.6120e-02,3.2829e-02,7.1771e-02,-3.3276e-02,
-3.9546e-02,-8.0963e-02,4.9240e-02,5.1887e-02,-6.6502e-03,-6.8269e-03,-7.5931e-04,-5.0691e-02,-5.7287e-02,3.8109e-02,9.5614e-03,6.7948e-02,
-5.0377e-02,-2.6256e-03,1.1003e-02,-1.2888e-02,1.9955e-03,-7.0537e-03,-2.4831e-02,1.3619e-02,1.9359e-02,-2.2920e-02,1.3439e-02,-3.3675e-02,
7.2703e-03,9.5540e-04,1.4802e-02,-7.6185e-03,2.4034e-02,-2.2186e-02,-1.0171e-02,3.1258e-02};

double SC_quad[80] = {
-8.3162e-02,-8.0222e-02,-1.2496e-02,-3.1370e-03,-2.0615e-03,-7.4027e-03,-8.0890e-03,-3.0603e-03,-6.2534e-03,8.3106e-03,3.0373e-03,7.7053e-04,
6.5205e-04,1.5969e-03,1.8320e-03,-2.4071e-03,-2.3424e-03,-3.6214e-03,2.4541e-03,-5.6141e-03,2.6180e-03,-1.9412e-03,-2.7277e-03,-2.7921e-03,
3.9476e-03,7.4313e-04,3.0627e-04,-1.4059e-03,-2.0170e-03,6.3919e-04,-1.1673e-03,-2.6028e-03,-1.3541e-03,3.0701e-03,-1.3899e-03,6.1153e-04,
6.0899e-04,1.4606e-03,2.3622e-05,8.3762e-04,-8.0877e-05,-2.7392e-05,5.3757e-04,1.6768e-05,-1.2602e-03,-4.1245e-04,-1.0847e-03,6.1458e-04,
6.7328e-04,1.2819e-03,-6.1516e-04,-5.8858e-04,1.8416e-04,1.9081e-04,1.4499e-04,7.4230e-04,8.0995e-04,-3.6410e-04,-1.3020e-05,-6.4771e-04,
6.4131e-04,1.4654e-04,-3.3370e-05,2.2911e-04,6.6761e-05,1.3327e-04,3.0475e-04,-7.4230e-05,-9.3400e-05,2.8012e-04,-4.8187e-05,3.7744e-04,
3.1476e-06,4.8946e-05,-5.1896e-05,1.3577e-04,-1.2482e-04,2.4817e-04,1.4184e-04,-1.7761e-04};

double SC_cubic[80] = {
2.9553e-02,1.2880e-02,2.2564e-03,6.7160e-04,3.5331e-04,5.3805e-04,4.8418e-04,1.8342e-04,2.8057e-04,-2.0495e-04,-3.6887e-05,1.6266e-05,
1.9660e-05,-7.2855e-06,-1.1872e-05,6.5040e-05,6.1960e-05,6.7565e-05,-2.9236e-05,8.8826e-05,-2.8006e-05,3.2922e-05,3.7934e-05,4.0166e-05,
-3.9658e-05,-1.8127e-06,1.6243e-06,1.6883e-05,1.9722e-05,-3.8674e-06,1.1874e-05,2.3737e-05,1.2521e-05,-2.2351e-05,1.2440e-05,-2.8343e-06,
-3.0883e-06,-9.5352e-06,6.9222e-07,-4.2526e-06,1.3480e-06,9.5643e-07,-2.4451e-06,8.2003e-07,7.2102e-06,2.7234e-06,6.5119e-06,-2.5842e-06,
-2.8135e-06,-5.7888e-06,3.3292e-06,2.9183e-06,-4.9355e-07,-5.8230e-07,-5.3101e-07,-2.8835e-06,-3.1213e-06,1.6429e-06,2.1131e-07,2.4677e-06,
-2.1770e-06,-4.8990e-07,2.4545e-07,-7.1359e-07,-1.4094e-07,-2.7333e-07,-8.2912e-07,3.9309e-07,3.3199e-07,-7.6214e-07,2.1628e-07,-1.0548e-06,
6.1077e-08,-4.9284e-08,1.8268e-07,-3.3613e-07,3.7203e-07,-6.3064e-07,-3.2549e-07,4.8505e-07};

#endif
