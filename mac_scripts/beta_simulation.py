import os


directory = "executable"


# beta simulation
nfiles = 2             #require 3 beta files in total
run = 10                 #number of runs one file is divided into
events_per_run = 800000

for j in range(1,nfiles):
    for i in range(0, run):
        filename = "Sr90_2.5iC4H10_trainsample_"
        command = f'''cat >{directory}/beta_run{j}_{i}.mac <<EOF
#Sets some gps default verbose
#
/mydet/generator fGParticleSource
/control/verbose 2
/run/verbose 1
/analysis/setFileName {filename}{j}_{i}
#
# For Co-60

/gps/particle ion
/gps/ion 38 90 0 0
/gps/energy    0
/gps/pos/type Plane
/gps/pos/shape Circle
/gps/pos/centre 0 0 -11.8 cm
/gps/pos/radius 0.5 cm
/gps/ang/surfnorm true
/gps/ang/type iso
/gps/ang/maxtheta 90. deg
#/globalField/setValue 0.1 0 0 tesla
#/gps/pos/rot1  1. 0 0
#/gps/pos/rot2  0 1. 0
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
#/analysis/h1/set 11 100  0. 1000 keV    #primarty particle spectrum
#/analysis/h1/set 12 100  0. 1000 keV    #secondary beta spectrum (for radioactive decay betas)
#/analysis/h1/set 13 100  -15. 15. cm    #x dim start position
#/analysis/h1/set 14 100  -15. 15. cm    #y dim start position
#
#/DBDecay/event/printModulo 100000
#  
/run/beamOn {events_per_run}
EOF'''
        os.system(command)

        command = f'''cat >{directory}/beta_simulation{j}_{i}.sh <<EOF
cd /home/rzhang/ustcfs/TPC_simulation/MMD_G4/NTD_Ge_new_version/build
rm Rawroot_{filename}{j}_{i}.root
./NTD_Ge {directory}/beta_run{j}_{i}.mac
EOF'''
        os.system(command)

        command = f'''cat >beta_simulation.condor <<EOF
universe = vanilla
Notification         = Never
GetEnv               = True
next_job_start_delay = 3
executable = {directory}/beta_simulation{j}_{i}.sh
output = {directory}/beta_output.txt
error = {directory}/beta_error.txt
log = {directory}/beta_log.txt
request_memory = 4096
priority = 1
queue
EOF'''
        os.system(command)
        command = f"chmod 764 {directory}/beta_simulation{j}_{i}.sh"
        os.system(command)
        command = f"condor_submit beta_simulation.condor"
        os.system(command)
    #command=f"sed -i 's/allpix_run{i}.sh/allpix_run{i+1}.sh/g' test.condor"
    # os.system(command)
    # run+=1

#command=f"sed -i 's/allpix_run{run}.sh/allpix_run0.sh/g' test.condor"
# os.system(command)
