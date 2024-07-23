import os


directory = "executable"


# muon backgrond simulation
nfiles = 1              #require 3 muon bkg files in total
run = 100                 #number of runs one file is divided into
events_per_run = 10000

for j in range(0,nfiles):
    for i in range(0, run):
        filename = "muon_30cm_30cm_coincident_efficiency_95_trainsample_"
        command = f'''cat >{directory}/muon_run{j}_{i}.mac <<EOF
/analysis/setFileName {filename}{j}_{i}
/tracking/verbose 0
/run/beamOn {events_per_run}
EOF'''
        os.system(command)

        command = f'''cat >{directory}/muon_simulation{j}_{i}.sh <<EOF
cd /home/rzhang/ustcfs/TPC_simulation/MMD_G4/NTD_Ge_new_version/build
rm Rawroot_{filename}{j}_{i}.root
./NTD_Ge {directory}/muon_run{j}_{i}.mac
EOF'''
        os.system(command)

        command = f'''cat >muon_simulation.condor <<EOF
universe = vanilla
Notification         = Never
GetEnv               = True
next_job_start_delay = 3
executable = {directory}/muon_simulation{j}_{i}.sh
output = {directory}/muon_output.txt
error = {directory}/muon_error.txt
log = {directory}/muon_log.txt
request_memory = 4096
priority = 1
queue
EOF'''
        os.system(command)
        command = f"chmod 764 {directory}/muon_simulation{j}_{i}.sh"
        os.system(command)
        command = f"condor_submit muon_simulation.condor"
        os.system(command)
    #command=f"sed -i 's/allpix_run{i}.sh/allpix_run{i+1}.sh/g' test.condor"
    # os.system(command)
    # run+=1

#command=f"sed -i 's/allpix_run{run}.sh/allpix_run0.sh/g' test.condor"
# os.system(command)
