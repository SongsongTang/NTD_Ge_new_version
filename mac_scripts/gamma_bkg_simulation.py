import os


directory = "executable"


# gamma backgrond simulation
nfiles = 2              #require 2 background files in total
run = 1500                 #number of runs one file is divided into
events_per_run = 10000000

for j in range(0,nfiles):
    for i in range(0, run):
        filename = "gamma_allaluminum_frame_5mmAl_2cmPb_"
        command = f'''cat >{directory}/gamma_run{j}_{i}.mac <<EOF
#Sets some gps default verbose
#
/mydet/generator fGParticleSource
/control/verbose 2
/run/verbose 1
#/tracking/storeTrajectory 1
/analysis/setFileName {filename}{j}_{i}
#
# For Co-60

#the following gamma flux data is from 104 lab's experiment result

/gps/particle gamma
/gps/ene/min 0 keV
/gps/ene/max 3100 keV
/gps/ene/type Arb
/gps/hist/type arb

/gps/hist/point 0.011839844	0
/gps/hist/point 0.023679688	0
/gps/hist/point 0.035519531	0
/gps/hist/point 0.047359375	0
/gps/hist/point 0.059199219	0
/gps/hist/point 0.071039063	0
/gps/hist/point 0.082878906	0
/gps/hist/point 0.09471875	1241679.822
/gps/hist/point 0.106558594	1133108.72
/gps/hist/point 0.118398438	1208353.396
/gps/hist/point 0.130238281	1250545.853
/gps/hist/point 0.142078125	1213497.358
/gps/hist/point 0.153917969	1118879.07
/gps/hist/point 0.165757813	1042089.051
/gps/hist/point 0.177597656	962093.1294
/gps/hist/point 0.1894375	1010931.906
/gps/hist/point 0.201277344	859789.5865
/gps/hist/point 0.213117188	833867.3605
/gps/hist/point 0.224957031	697457.0268
/gps/hist/point 0.236796875	842516.4924
/gps/hist/point 0.248636719	987712.8802
/gps/hist/point 0.260476563	451502.0678
/gps/hist/point 0.272316406	551417.7406
/gps/hist/point 0.28415625	493551.1777
/gps/hist/point 0.295996094	588389.5417
/gps/hist/point 0.307835938	455538.256
/gps/hist/point 0.319675781	338619.997
/gps/hist/point 0.331515625	371487.7249
/gps/hist/point 0.343355469	336704.8147
/gps/hist/point 0.355195313	670089.4554
/gps/hist/point 0.367035156	261188.8087
/gps/hist/point 0.378875	230565.1875
/gps/hist/point 0.390714844	257808.3375
/gps/hist/point 0.402554688	269518.2476
/gps/hist/point 0.414394531	255030.3696
/gps/hist/point 0.426234375	150588.7756
/gps/hist/point 0.438074219	310295.0882
/gps/hist/point 0.449914063	130176.7136
/gps/hist/point 0.461753906	250498.2598
/gps/hist/point 0.47359375	247819.3038
/gps/hist/point 0.485433594	107109.4438
/gps/hist/point 0.497273438	133685.1218
/gps/hist/point 0.509113281	283904.731
/gps/hist/point 0.520953125	264508.0578
/gps/hist/point 0.532792969	131939.2394
/gps/hist/point 0.544632813	128241.3972
/gps/hist/point 0.556472656	156307.4285
/gps/hist/point 0.5683125	56770.56874
/gps/hist/point 0.580152344	279768.0304
/gps/hist/point 0.591992188	406872.8533
/gps/hist/point 0.603832031	172885.0225
/gps/hist/point 0.615671875	652947.9284
/gps/hist/point 0.627511719	56025.63659
/gps/hist/point 0.639351563	24150.64528
/gps/hist/point 0.651191406	60423.86362
/gps/hist/point 0.66303125	129655.9997
/gps/hist/point 0.674871094	173900.2704
/gps/hist/point 0.686710938	36829.70747
/gps/hist/point 0.698550781	91918.24018
/gps/hist/point 0.710390625	82472.15055
/gps/hist/point 0.722230469	95338.11125
/gps/hist/point 0.734070313	106527.1584
/gps/hist/point 0.745910156	116623.6501
/gps/hist/point 0.75775	155286.1472
/gps/hist/point 0.769589844	198803.4829
/gps/hist/point 0.781429688	181916.6328
/gps/hist/point 0.793269531	46377.17138
/gps/hist/point 0.805109375	98930.99375
/gps/hist/point 0.816949219	81773.99937
/gps/hist/point 0.828789063	56658.99872
/gps/hist/point 0.840628906	153878.7929
/gps/hist/point 0.85246875	121295.0744
/gps/hist/point 0.864308594	102254.1374
/gps/hist/point 0.876148438	13299.49259
/gps/hist/point 0.887988281	53203.928
/gps/hist/point 0.899828125	118738.6873
/gps/hist/point 0.911667969	589070.1961
/gps/hist/point 0.923507813	52716.05044
/gps/hist/point 0.935347656	92195.59173
/gps/hist/point 0.9471875	34337.13774
/gps/hist/point 0.959027344	143492.526
/gps/hist/point 0.970867188	353016.134
/gps/hist/point 0.982707031	65558.94131
/gps/hist/point 0.994546875	35216.03562
/gps/hist/point 1.006386719	35872.02664
/gps/hist/point 1.018226563	30037.60012
/gps/hist/point 1.030066406	129982.7752
/gps/hist/point 1.04190625	77534.71403
/gps/hist/point 1.053746094	40627.55301
/gps/hist/point 1.065585938	92694.89043
/gps/hist/point 1.077425781	73456.73171
/gps/hist/point 1.089265625	13560.76619
/gps/hist/point 1.101105469	94942.26549
/gps/hist/point 1.112945313	86178.14058
/gps/hist/point 1.124785156	146965.0622
/gps/hist/point 1.136625	45155.8418
/gps/hist/point 1.148464844	41564.29026
/gps/hist/point 1.160304688	38807.42923
/gps/hist/point 1.172144531	74562.69037
/gps/hist/point 1.183984375	18230.40451
/gps/hist/point 1.195824219	140017.3666
/gps/hist/point 1.207664063	47215.9451
/gps/hist/point 1.219503906	29269.02825
/gps/hist/point 1.23134375	81297.30273
/gps/hist/point 1.243183594	110919.8976
/gps/hist/point 1.255023438	64436.84557
/gps/hist/point 1.266863281	92504.99424
/gps/hist/point 1.278703125	142639.0434
/gps/hist/point 1.290542969	39383.10793
/gps/hist/point 1.302382813	18764.02281
/gps/hist/point 1.314222656	1488.186678
/gps/hist/point 1.3260625	6757.612551
/gps/hist/point 1.337902344	10071.71567
/gps/hist/point 1.349742188	9403.403336
/gps/hist/point 1.361582031	54066.04274
/gps/hist/point 1.373421875	100936.9618
/gps/hist/point 1.385261719	22961.52518
/gps/hist/point 1.397101563	46916.19147
/gps/hist/point 1.408941406	133717.6526
/gps/hist/point 1.42078125	118450.6366
/gps/hist/point 1.432621094	88036.69266
/gps/hist/point 1.444460938	40980.5905
/gps/hist/point 1.456300781	1924569.508
/gps/hist/point 1.468140625	160759.5757
/gps/hist/point 1.479980469	2690.523408
/gps/hist/point 1.491820313	3738.427572
/gps/hist/point 1.503660156	19856.94286
/gps/hist/point 1.5155	1104.083645
/gps/hist/point 1.527339844	217.092321
/gps/hist/point 1.539179688	888.064636
/gps/hist/point 1.551019531	23270.33229
/gps/hist/point 1.562859375	41215.0193
/gps/hist/point 1.574699219	47582.79362
/gps/hist/point 1.586539063	118403.0698
/gps/hist/point 1.598378906	56237.28094
/gps/hist/point 1.61021875	131382.2644
/gps/hist/point 1.622058594	52144.89588
/gps/hist/point 1.633898438	8766.52304
/gps/hist/point 1.645738281	116.226135
/gps/hist/point 1.657578125	56.793693
/gps/hist/point 1.669417969	413.737095
/gps/hist/point 1.681257813	6241.445173
/gps/hist/point 1.693097656	23694.71495
/gps/hist/point 1.7049375	90813.40864
/gps/hist/point 1.716777344	35115.24219
/gps/hist/point 1.728617188	15002.47205
/gps/hist/point 1.740457031	97590.84611
/gps/hist/point 1.752296875	243947.8443
/gps/hist/point 1.764136719	66041.25462
/gps/hist/point 1.775976563	24473.06703
/gps/hist/point 1.787816406	26385.29391
/gps/hist/point 1.79965625	7758.873276
/gps/hist/point 1.811496094	10211.60349
/gps/hist/point 1.823335938	17283.84297
/gps/hist/point 1.835175781	18124.09071
/gps/hist/point 1.847015625	4778.968881
/gps/hist/point 1.858855469	19333.71157
/gps/hist/point 1.870695313	15198.66212
/gps/hist/point 1.882535156	4752.781234
/gps/hist/point 1.894375	7146.42286
/gps/hist/point 1.906214844	2803.173027
/gps/hist/point 1.918054688	1653.286752
/gps/hist/point 1.929894531	1389.230594
/gps/hist/point 1.941734375	733.581724
/gps/hist/point 1.953574219	929.079176
/gps/hist/point 1.965414063	3253.653662
/gps/hist/point 1.977253906	1109.639457
/gps/hist/point 1.98909375	395.462882
/gps/hist/point 2.000933594	78.798385
/gps/hist/point 2.012773438	321.271617
/gps/hist/point 2.024613281	4585.644339
/gps/hist/point 2.036453125	2011.691019
/gps/hist/point 2.048292969	3658.774923
/gps/hist/point 2.060132813	12076.68144
/gps/hist/point 2.071972656	5708.257503
/gps/hist/point 2.0838125	10133.45367
/gps/hist/point 2.095652344	49864.50627
/gps/hist/point 2.107492188	40916.37799
/gps/hist/point 2.119332031	8525.357687
/gps/hist/point 2.131171875	7289.067582
/gps/hist/point 2.143011719	7067.597377
/gps/hist/point 2.154851563	1702.668639
/gps/hist/point 2.166691406	560.18554
/gps/hist/point 2.17853125	117.23941
/gps/hist/point 2.190371094	80.267022
/gps/hist/point 2.202210938	44.191212
/gps/hist/point 2.214050781	51.453533
/gps/hist/point 2.225890625	47.349683
/gps/hist/point 2.237730469	26.56
/gps/hist/point 2.249570313	56.394677
/gps/hist/point 2.261410156	98.70112
/gps/hist/point 2.27325	52.764313
/gps/hist/point 2.285089844	71.046823
/gps/hist/point 2.296929688	775.375353
/gps/hist/point 2.308769531	1155.700428
/gps/hist/point 2.320609375	349.667262
/gps/hist/point 2.332449219	118.764387
/gps/hist/point 2.344289063	263.166139
/gps/hist/point 2.356128906	1959.235261
/gps/hist/point 2.36796875	790.584006
/gps/hist/point 2.379808594	979.63204
/gps/hist/point 2.391648438	2885.901089
/gps/hist/point 2.403488281	329.156396
/gps/hist/point 2.415328125	61.097067
/gps/hist/point 2.427167969	268.771352
/gps/hist/point 2.439007813	261.175393
/gps/hist/point 2.450847656	195.587923
/gps/hist/point 2.4626875	34.07114
/gps/hist/point 2.474527344	709.938699
/gps/hist/point 2.486367188	181.628961
/gps/hist/point 2.498207031	1.144952
/gps/hist/point 2.510046875	0.867718
/gps/hist/point 2.521886719	11.845484
/gps/hist/point 2.533726563	17.484288
/gps/hist/point 2.545566406	2.933168
/gps/hist/point 2.55740625	19.316477
/gps/hist/point 2.569246094	818.301328
/gps/hist/point 2.581085938	8250.414839
/gps/hist/point 2.592925781	60749.48628
/gps/hist/point 2.604765625	528876.4054
/gps/hist/point 2.616605469	121207.9059
/gps/hist/point 2.628445313	7888.040583
/gps/hist/point 2.640285156	875.109832
/gps/hist/point 2.652125	253.925493
/gps/hist/point 2.663964844	0.335662
/gps/hist/point 2.675804688	0.005038
/gps/hist/point 2.687644531	0.006795
/gps/hist/point 2.699484375	0.194361
/gps/hist/point 2.711324219	1.109926
/gps/hist/point 2.723164063	13.95665
/gps/hist/point 2.735003906	7.108394
/gps/hist/point 2.74684375	0.249363
/gps/hist/point 2.758683594	1.624757
/gps/hist/point 2.770523438	93.43647
/gps/hist/point 2.782363281	474.483109
/gps/hist/point 2.794203125	1083.680571
/gps/hist/point 2.806042969	6008.213111
/gps/hist/point 2.817882813	14237.83085
/gps/hist/point 2.829722656	29058.90105
/gps/hist/point 2.8415625	9389.906888
/gps/hist/point 2.853402344	2040.968258
/gps/hist/point 2.865242188	958.626842
/gps/hist/point 2.877082031	9467.736434
/gps/hist/point 2.888921875	75592.94578
/gps/hist/point 2.900761719	5747.187176
/gps/hist/point 2.912601563	7207.977109
/gps/hist/point 2.924441406	8739.422257
/gps/hist/point 2.93628125	51993.18799
/gps/hist/point 2.948121094	6970.414491
/gps/hist/point 2.959960938	2686.367073
/gps/hist/point 2.971800781	35938.46868
/gps/hist/point 2.983640625	102139.7951
/gps/hist/point 2.995480469	34626.6102
/gps/hist/point 3.007320313	12507.94081
/gps/hist/point 3.019160156	43714.41602
/gps/hist/point 3.031	139099.2431


/gps/hist/inter Lin

/gps/pos/type Surface
/gps/pos/shape Sphere
/gps/pos/centre 0 0 0 cm
/gps/pos/radius 30 cm
#/gps/pos/rot1  1. 0 0
#/gps/pos/rot2  0 1. 0
/gps/ang/surfnorm true
/gps/ang/type cos
/gps/ang/mintheta 0. deg
/gps/ang/maxtheta 90. deg
#/gps/direction 0 0 1.
#
#/gps/direction  0. 0. 1.0
#/gps/ene/tupe Brem
#/gps/ene/temp 2e9
#/gps/ene/min  1. keV
#/gps/ene/max  2. MeV
#
# /tracking/verbose 2
# /run/beamOn 1
#
/tracking/verbose 0
#
#/analysis/h1/set 1  150  0. 1500 keV	#e+ e-
#/analysis/h1/set 2  150  0. 1500 keV	#neutrino
#/analysis/h1/set 3  150  0. 1500 keV	#gamma
#/analysis/h1/set 6  100  0. 2500 keV	#EkinTot (Q)
#/analysis/h1/set 7  150  0. 15e3 keV	#P balance
#/analysis/h1/set 8  100  0. 100. year	#time of life
#/analysis/h1/set 9  150  0. 600. keV    #e- spectrum in Gas
#/analysis/h1/set 10 100  0. 100. cm     #e- track length in Gas

# /analysis/h1/set 11 256  0. 3050. keV    #primarty particle spectrum
# /analysis/h1/set 13 100  -15. 15. cm    #x dim start position
# /analysis/h1/set 14 100  -15. 15. cm    #y dim start position

# /analysis/h1/set 15 256 0. 3050. keV    #primarty particle spectrum after selection
# /analysis/h1/set 16 256 0. 3050. keV    #primarty particle spectrum after selection from PCB
# /analysis/h1/set 17 256 0. 3050. keV    #primarty particle spectrum after selection from Gas
# /analysis/h1/set 18 256 0. 3050. keV    #primarty particle spectrum after selection from film
# /analysis/h1/set 19 256 0. 3050. keV    #primarty particle spectrum after selection from Cu board
# /analysis/h1/set 20 256 0. 3050. keV    #primarty particle spectrum after selection from other places
# /analysis/h1/set 21 256 0. 3050. keV    #primarty particle spectrum after selection from compt
# /analysis/h1/set 22 256 0. 3050. keV    #primarty particle spectrum after selection from phot
# /analysis/h1/set 23 256 0. 3050. keV    #primarty particle spectrum after selection from eIoni
# /analysis/h1/set 24 256 0. 3050. keV    #primarty particle spectrum after selection from other processes


#
#/DBDecay/event/printModulo 100000
#  
#/run/beamOn 265500579
#/run/beamOn 132750289
/run/beamOn {events_per_run}
EOF'''
        os.system(command)

        command = f'''cat >{directory}/gamma_simulation{j}_{i}.sh <<EOF
cd /home/rzhang/ustcfs/TPC_simulation/MMD_G4/NTD_Ge_new_version/build
rm Rawroot_{filename}{j}_{i}.root
./NTD_Ge {directory}/gamma_run{j}_{i}.mac
EOF'''
        os.system(command)

        command = f'''cat >gamma_simulation.condor <<EOF
universe = vanilla
Notification         = Never
GetEnv               = True
next_job_start_delay = 3
executable = {directory}/gamma_simulation{j}_{i}.sh
output = {directory}/gamma_output.txt
error = {directory}/gamma_error.txt
log = {directory}/gamma_log.txt
request_memory = 4096
# priority = 1
queue
EOF'''
        os.system(command)
        command = f"chmod 764 {directory}/gamma_simulation{j}_{i}.sh"
        os.system(command)
        command = f"condor_submit gamma_simulation.condor"
        os.system(command)
    #command=f"sed -i 's/allpix_run{i}.sh/allpix_run{i+1}.sh/g' test.condor"
    # os.system(command)
    # run+=1

#command=f"sed -i 's/allpix_run{run}.sh/allpix_run0.sh/g' test.condor"
# os.system(command)
