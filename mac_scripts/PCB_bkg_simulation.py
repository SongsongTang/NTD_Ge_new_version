import os


directory = "executable"


# beta simulation
nfiles = 2              #require 3 beta files in total
run = 10                 #number of runs one file is divided into
# these event numbers are corresponding to 4e8 gamma bkgg events
events_Th232 = 5698
events_U238 = 3538
events_K40 = 2431

for j in range(0,nfiles):
    for i in range(0, run):
        filename = "PCBbkg_1of2_digi_1cmPb_"
        command = f'''cat >{directory}/PCBbkg_run{j}_{i}.mac <<EOF
#Sets some gps default verbose
#
/mydet/generator fGParticleSource
/control/verbose 2
/run/verbose 1
/analysis/setFileName {filename}{j}_{i}

# ==== from the central PCB =====

/gps/particle ion
/gps/pos/type Volume
/gps/pos/shape Para
/gps/pos/centre 0 0 0.115 cm
/gps/pos/halfx 11.65 cm
/gps/pos/halfy 11.65 cm
/gps/pos/halfz 0.115 cm
/gps/ang/type iso

# Th-232
/gps/ion 90 232 0 0
/gps/energy    0
/tracking/verbose 0
/run/beamOn {events_Th232}

# U-238
/gps/ion 92 238 0 0
/gps/energy    0
/tracking/verbose 0
/run/beamOn {events_U238}

# K-40
/gps/ion 19 40 0 0
/gps/energy    0
/tracking/verbose 0
/run/beamOn {events_K40}

EOF'''
        os.system(command)

        command = f'''cat >{directory}/PCBbkg_simulation{j}_{i}.sh <<EOF
cd /home/rzhang/ustcfs/TPC_simulation/MMD_G4/NTD_Ge_new_version/build
rm Rawroot_{filename}{j}_{i}.root
./NTD_Ge {directory}/PCBbkg_run{j}_{i}.mac
EOF'''
        os.system(command)

        command = f'''cat >PCBbkg_simulation.condor <<EOF
universe = vanilla
Notification         = Never
GetEnv               = True
next_job_start_delay = 3
executable = {directory}/PCBbkg_simulation{j}_{i}.sh
output = {directory}/PCB_output.txt
error = {directory}/PCB_error.txt
log = {directory}/PCB_log.txt
request_memory = 4096
priority = 1
queue
EOF'''
        os.system(command)
        command = f"chmod 764 {directory}/PCBbkg_simulation{j}_{i}.sh"
        os.system(command)
        command = f"condor_submit PCBbkg_simulation.condor"
        os.system(command)
    #command=f"sed -i 's/allpix_run{i}.sh/allpix_run{i+1}.sh/g' test.condor"
    # os.system(command)
    # run+=1

#command=f"sed -i 's/allpix_run{run}.sh/allpix_run0.sh/g' test.condor"
# os.system(command)
