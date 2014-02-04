C
C Subroutine to return the Hydrogenic gaunt factors in the COMMON block
C HYDCROSS.
C
C No other action is taken by the routine
C
C Altered 09-May-00 - TMP_HYD inserted to avoid changing to a BLOCK data
C                       subroutine.
C Altered 27-Nov-91 - Changed from function to subroutine.
C Altered 26-Mar-88 - Extended to n=6 although g(5s,5p,6s,6p) set to 1.
C
	SUBROUTINE HYDANGGAUNT
	COMMON /HYDCROSS/ HYD(0:24,1:18)
	REAL*8 HYD
	REAL*8 TMP_HYD(0:24,1:18)
C
C The first 50 elements correspond to 1s and 2s, and are never used
C hence they are set to Zero.
C
	DATA TMP_HYD/50*0.00,
C 2p
	9    -0.0672,-0.0672,-0.0672,-0.0672,-0.0673,-0.0673,-0.0675
	9 , -0.0678,-0.0683,-0.0693,-0.0713,-0.0754,-0.0842,-0.1021
	9 , -0.1366,-0.1969,-0.2919,-0.4271,-0.6037,-0.8190,-1.0681
	9 , -1.3451,-1.6449,-1.9623,-2.2935,
C 3s
	9     0.0261, 0.0263, 0.0268, 0.0278, 0.0297, 0.0328, 0.0377, 0.0464
	9 , 0.0618, 0.0873, 0.1284, 0.1909, 0.2774, 0.3866, 0.5098, 0.6331
	9 , 0.7423, 0.8265, 0.8802, 0.9027, 0.8963, 0.8649, 0.8126, 0.7438
	9 , 0.6617,
C 3p
	9     0.0406, 0.0408, 0.0412, 0.0418, 0.0428, 0.0442, 0.0472, 0.0524
	9 , 0.0607, 0.0745, 0.0954, 0.1242, 0.1568, 0.1845, 0.1930, 0.1671
	9 , 0.0966,-0.0230,-0.1902,-0.4001, -0.6461,-0.9215,-1.2201,-1.5371
	9 ,-1.8681,
C 3d
	9    -0.1179,-0.1180,-0.1181,-0.1186,-0.1194,-0.1207,-0.1231,-0.1272
	9 ,-0.1347,-0.1479,-0.1711,-0.2114,-0.2795,-0.3883,-0.5510,-0.7764
	9 ,-1.0674,-1.4208,-1.8295,-2.2850,-2.7791,-3.3030,-3.8512,-4.4179
	9 ,-4.9987,
C 4s
	9    0.0722, 0.0726, 0.0734, 0.0752, 0.0784, 0.0836, 0.0922, 0.1072
	9 , 0.1329, 0.1743, 0.2374, 0.3273, 0.4431, 0.5790, 0.7229, 0.8606
	9 , 0.9789, 1.0686, 1.1259, 1.1502, 1.1447, 1.1139, 1.0622, 0.9934
	9 , 0.9115,
C 4p
	9    0.0962, 0.0965, 0.0971, 0.0985, 0.1008, 0.1047, 0.1115, 0.1231
	9 , 0.1424, 0.1729, 0.2169, 0.2749, 0.3389, 0.3952, 0.4263, 0.4161
	9 , 0.3555, 0.2419, 0.0783,-0.1296,-0.3745,-0.6492,-0.9476,-1.2643
	9 ,-1.5952
C 4d
	1 ,0.0535,0.0537,0.0539,0.0542,0.0548,0.0559,0.0577
	1 ,0.0604,0.0641,0.0685,0.0698,0.0603,0.0255,-0.0530,-0.1914
	1 ,-0.3998,-0.6799,-1.0266,-1.4314,-1.8846,-2.3773,-2.9008
	1 ,-3.4484,-4.0148,-4.5953
C 4f
	1 ,-0.2200,-0.2222,-0.2229,-0.2233,-0.2250,-0.2290,-0.2366
	1 ,-0.2498,-0.2729,-0.3128,-0.3804,-0.4908,-0.6615,-0.9098
	1 ,-1.2458,-1.6720,-2.1827,-2.7686,-3.4167,-4.1164,-4.8570
	1 ,-5.6292,-6.4264,-7.2427,-8.0728
C 5s
	1 ,25*0.00
C 5p
	1 ,25*0.00
C 5d
	1 ,0.1316,0.1330,0.1336,0.1343,0.1362,0.1401,0.1468,0.1576
	1 ,0.1738,0.1956,0.2196,0.2353,0.2242,0.1647,0.0400,-0.1594
	1 ,-0.4340,-0.7778,-1.1804,-1.6328,-2.1247,-2.6478,-3.1955
	1 ,-3.7616,-4.3421
C 5fg
	1 ,-0.1538,-0.1557,-0.1563,-0.1568,-0.1584,-0.1624,-0.1699
	1 ,-0.1828,-0.2055,-0.2452,-0.3129,-0.4240,-0.5962,-0.8461
	1 ,-1.1834,-1.6105,-2.1220,-2.7084,-3.3568,-4.0567,-4.7978
	1 ,-5.5715,-6.3499,-7.2789,-8.8584
C 6s
	1 ,25*0.00
C 6p
	1 ,25*0.00
C 6d
	1 ,0.1801,0.1826,0.1839,0.1855,0.1891,0.1961,0.2080,0.2272
	1 ,0.2562,0.2952,0.3390,0.3741,0.3795,0.3326,0.2163,0.0219
	1 ,-0.2497,-0.5916,-0.9934,-1.4450,-1.9367,-2.4596,-3.0070
	1 ,-3.5734,-4.1537
C 6fgh
	1 ,-0.1179,-0.1201,-0.1207,-0.1211,-0.1227,-0.1266,-0.1336
	1 ,-0.1469,-0.1698,-0.2078,-0.2820,-0.3663,-0.5442,-0.9124
	1 ,-0.9628,-2.2915,-2.0916,-2.6727,-3.3259,-4.0265
	1 ,-4.7676,-5.5414,-6.3198,-7.2490,-8.8284/
C
	HYD=TMP_HYD
C
	END