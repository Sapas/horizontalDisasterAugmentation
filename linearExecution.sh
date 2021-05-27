# Run this script to run the augmentation algorithms as specified in the runScripts folder
# For now will comment / uncomment the different options

# NOTE: please make sure you first run the setup script to make sure everything is as should be

# ATTENTION: Only change the values under this statement

runName="smallTest"
lines=48

# runName="preliminaryRun"
# lines=880

# runName="mainRun"
# lines=576

# runName="mainRunLarge"
# lines=864

print_info=0
create_graphs=0

# NO CHANGES PAST THIS POINT

# Go into code folder
cd cpp_code/
# Get time
start=`date +%s`

# Cycle through lines in runScript
for ((n = 1; n <= $lines; n = n + 1)); do
	./augmentation run ${runName} ${n} ${print_info} ${create_graphs}
done

end=`date +%s`
echo "Total execution time: $((end-start)) seconds"
