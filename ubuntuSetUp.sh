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
