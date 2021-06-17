# cd /data/gpfs/projects/punim1489/horizontalDisasterAugmentation/
# First get modules needed for cgal to work
module load gcc/8.3.0
module load cmake/3.18.4 
module load cgal/4.14.1-python-3.7.4
module load gmp/6.1.2
module load mpfr/4.0.2
# Navigate to folder with the code
cd cpp_code/
# Create CMakeLists
/usr/local/easybuild-2019/easybuild/software/mpi/gcc/8.3.0/openmpi/3.1.4/cgal/4.14.1-python-3.7.4/bin/cgal_create_CMakeLists -c Core -s augmentation
# Run the cmake script, and then make the program
cmake -DCMAKE_BUILD_TYPE=Release -DGMP_LIBRARIES=$EBROOTGMP/lib -DGMP_INCLUDE_DIR=$EBROOTGMP/include -DMPFR_LIBRARIES=$EBROOTMPFR/lib -DMPFR_INCLUDE_DIR=$EBROOTMPFR/include .
make
# For some weird reason, need this last call to make it all work
c++ -O3 -DNDEBUG -rdynamic -lCGAL -lgmp -lboost_thread -lmpfr CMakeFiles/augmentation.dir/augmentation.cpp.o CMakeFiles/augmentation.dir/cgal_interface.cpp.o CMakeFiles/augmentation.dir/input_output.cpp.o CMakeFiles/augmentation.dir/line_geometry.cpp.o CMakeFiles/augmentation.dir/simple_svg_1.0.0.cpp.o CMakeFiles/augmentation.dir/sweep_binary_search_tree.cpp.o -o augmentation
# These are the standard run scripts generated in other platforms, should not be necessary to run them if files are already present
#./augmentation setup smallTest \{10,20,5\} \{1,5\} \{M,T\} \{10,200\} \{o2,o3,q2,q4\} 500 500
#./augmentation setup preliminaryRun \{10,100,10\} \{1,25\} \{M,T\} \{0.001,5,100,300\} \{a,b,c,d,e,i,j,k,o2,p,s\} 500 500
#./augmentation setup mainRun \{50,100,10\} \{1,25\} \{M,T\} \{0.001,5,100,300\} \{o1,o2,o3,o4,q1,q2,q3,q4,u1,u2,u3,u4\} 500 500
#./augmentation setup mainRunLarge \{100,500,50\} \{1,25\} \{M,T\} \{0.001,5,100,300\} \{o1,o2,o3,o4,q1,q2,q3,q4,u1,u2,u3,u4\} 500 500