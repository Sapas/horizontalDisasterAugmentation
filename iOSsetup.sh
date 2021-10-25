# First install CGAL
brew install cgal
# Navigate to folder with the code
cd cpp_code/
# Create CMakeLists
# Note if this step does not work, it is possible cgal was installed in a different folder. For help look at
# https://doc.cgal.org/latest/Manual/usage.html
/Users/nandresthio/Documents/CGAL-5.0/scripts/cgal_create_CMakeLists -c Core -s augmentation
# Create cmake lists and make the project
cmake -DCMAKE_BUILD_TYPE=Release .
make

# Run program to create the desired run files
#./augmentation setup smallTest \{10,20,5\} \{1,5\} \{M,T\} \{10,200\} \{o2,o3,q2,q4\} 500 500
#./augmentation setup preliminaryRun \{10,100,10\} \{1,25\} \{M,T\} \{0.001,5,100,300\} \{a,b,c,d,e,i,j,k,o2,p,s\} 500 500
#./augmentation setup mainRun \{10,100,10\} \{1,100\} \{M,T\} \{0.001,5,100,300\} \{o1,o2,o3,o4,q1,q2,q3,q4,u1,u2,u3,u4\} 500 500
#./augmentation setup polynomialRun \{10,100,10\} \{1,100\} \{M,T\} \{0.001,5,100,300\} \{v1,v2,v3,v4,w1,w2,w3,w4\} 500 500
#./augmentation setup intermediateWeightsRun \{10,100,10\} \{1,100\} \{M,T\} \{0.001,5,100,300\} \{o1.5,o2.5,o3.5,q1.5,q2.5,q3.5,u1.5,u2.5,u3.5,v1.5,v2.5,v3.5,w1.5,w2.5,w3.5\} 500 500
#./augmentation setup mainRunLarge \{100,500,50\} \{1,25\} \{M,T\} \{0.001,5,100,300\} \{o1,o2,o3,o4,q1,q2,q3,q4,u1,u2,u3,u4\} 500 500
#./augmentation setup polynomialComparisonRun \{10,100,10\} \{1,100\} \{M,T\} \{0.001,5,100,300\} \{o3,q2,v3,w2\} 500 500

# Run program to create the graphs
#./augmentation graphCreation preliminaryRun
#./augmentation graphCreation mainRun
#./augmentation graphCreation polynomialRun
#./augmentation graphCreation intermediateWeightsRun
#./augmentation graphCreation mainRunLarge
#./augmentation graphCreation polynomialComparisonRun


./augmentation setup largeInstancesTest \{750,750,250\} \{1,5\} \{M\} \{100\} \{q2\} 500 500
#./augmentation graphCreation largeInstancesTest

runName="largeInstancesTest"
lines=6

print_info=1
create_graphs=0

# NO CHANGES PAST THIS POINT

# Get time
start=`date +%s`

# First create the required graphs
./augmentation graphCreation ${runName}
# Cycle through lines in runScript
for ((n = 1; n <= $lines; n = n + 1)); do
	./augmentation run ${runName} ${n} ${print_info} ${create_graphs}
done

end=`date +%s`
echo "Total execution time: $((end-start)) seconds"

