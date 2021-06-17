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
./augmentation setup preliminaryRun \{10,100,10\} \{1,25\} \{M,T\} \{0.001,5,100,300\} \{a,b,c,d,e,i,j,k,o2,p,s\} 500 500
./augmentation setup mainRun \{10,100,10\} \{1,25\} \{M,T\} \{0.001,5,100,300\} \{o1,o2,o3,o4,q1,q2,q3,q4,u1,u2,u3,u4\} 500 500
#./augmentation setup mainRunLarge \{100,500,50\} \{1,25\} \{M,T\} \{0.001,5,100,300\} \{o1,o2,o3,o4,q1,q2,q3,q4,u1,u2,u3,u4\} 500 500
./augmentation setup polynomialRun \{10,100,10\} \{1,25\} \{M,T\} \{0.001,5,100,300\} \{v1,v2,v3,v4,w1,w2,w3,w4\} 500 500

# Run program to create the graphs
#./augmentation graphCreation smallTest
./augmentation graphCreation preliminaryRun
./augmentation graphCreation mainRun
#./augmentation graphCreation mainRunLarge
./augmentation graphCreation polynomialRun


