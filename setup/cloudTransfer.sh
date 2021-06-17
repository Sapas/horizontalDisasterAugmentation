# Here are a set of commands I keep on calling to use scp (transferring between the cluster and my machine)
# It is important I do not forget them so I am writing them here


# ================================================    CLUSTER   ================================================    
# First of all want to log into the cluster
ssh nandres@spartan.hpc.unimelb.edu.au
# Should be prompted by password
# Next will want to go to my project folder
cd /data/gpfs/projects/punim1489/horizontalDisasterAugmentation/

# First transfer everything as a folder, don't expect to do this in the future
# Go to horizontalDisasterAugmentation, then call
scp -r . nandres@spartan.hpc.unimelb.edu.au:/data/gpfs/projects/punim1489/horizontalDisasterAugmentation/ 

# In order to move code into the cpp folder, go to horizontalDisasterAugmentation/cpp_code/ and call
scp augmentation.cpp nandres@spartan.hpc.unimelb.edu.au:/data/gpfs/projects/punim1489/horizontalDisasterAugmentation/cpp_code/


# Might want to transfer setup file, in that case navigate to horizontalDisasterAugmentation/ and call
scp clusterSetup.sh nandres@spartan.hpc.unimelb.edu.au:/data/gpfs/projects/punim1489/horizontalDisasterAugmentation/

# Will want to transfer all the graphs so that the cluster does not create them and only reads them
# Go to horizontalDisasterAugmentation/data/graphs/ and call
scp -r . nandres@spartan.hpc.unimelb.edu.au:/data/gpfs/projects/punim1489/horizontalDisasterAugmentation/data/graphs/

# If want to transfer the cluster run information, need to go to horizontalDisasterAugmentation/data/runScripts and call
scp mainRun-runScript.txt nandres@spartan.hpc.unimelb.edu.au:/data/gpfs/projects/punim1489/horizontalDisasterAugmentation/data/runScripts/
scp mainRunLarge-runScript.txt nandres@spartan.hpc.unimelb.edu.au:/data/gpfs/projects/punim1489/horizontalDisasterAugmentation/data/runScripts/
scp preliminaryRun-runScript.txt nandres@spartan.hpc.unimelb.edu.au:/data/gpfs/projects/punim1489/horizontalDisasterAugmentation/data/runScripts/
scp smallTest-runScript.txt nandres@spartan.hpc.unimelb.edu.au:/data/gpfs/projects/punim1489/horizontalDisasterAugmentation/data/runScripts/

# Might want to transfer the cluster run information, to do so go to horizontalDisasterAugmentation/setup/ and call
scp clusterRun.sh nandres@spartan.hpc.unimelb.edu.au:/data/gpfs/projects/punim1489/horizontalDisasterAugmentation/setup/
scp longClusterRun.sh nandres@spartan.hpc.unimelb.edu.au:/data/gpfs/projects/punim1489/horizontalDisasterAugmentation/setup/

# In order to transfer the data, go to horizontalDisasterAugmentation/data/ and call
scp -r nandres@spartan.hpc.unimelb.edu.au:/data/gpfs/projects/punim1489/horizontalDisasterAugmentation/data/runResults/ .



# =================================================    CLOUD   =================================================   
# Could technically use github, but want to fix memory leaks, so will be doing it here
# If want to transfer C++ files, go to horizontalDisasterAugmentation/cpp_code/ and call
scp -i ../../access/euclideanSkeletonsKey.pem augmentation.cpp ubuntu@45.113.234.177:/home/ubuntu/horizontalDisasterAugmentation/cpp_code/






