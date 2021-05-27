# First install dependencies of CGAL and CGAL itself
# g++
sudo apt install g++
# cmake
sudo apt install cmake
# boost
sudo apt-get install libboost-all-dev
# cgal itself
sudo apt-get install libcgal-dev
# Now navigate to folder containing code
cd cpp_code
# Create CMakeLists
# Note if this step does not work, it is possible cgal was installed in a different folder. For help look at
# https://doc.cgal.org/latest/Manual/usage.html
/usr/bin/cgal_create_CMakeLists -c Core -s augmentation
# Run the cmake script, and then make the program
cmake -DCMAKE_BUILD_TYPE=Release .
make
# Run program to create the desired run files
./augmentation setup smallTest \{10,20,5\} \{1,5\} \{M,T\} \{10,200\} \{o2,o3,q2,q4\} 500 500
./augmentation setup preliminaryRun \{10,100,10\} \{1,25\} \{M,T\} \{0.001,5,100,300\} \{a,b,c,d,e,i,j,k,o2,p,s\} 500 500
./augmentation setup mainRun \{50,100,10\} \{1,25\} \{M,T\} \{0.001,5,100,300\} \{o1,o2,o3,o4,p1,p2,p3,p4,u1,u2,u3,u4\} 500 500
./augmentation setup mainRunLarge \{100,500,50\} \{1,25\} \{M,T\} \{0.001,5,100,300\} \{o1,o2,o3,o4,p1,p2,p3,p4,u1,u2,u3,u4\} 500 500

