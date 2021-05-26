// FOR NOW JUST HAVE COMMANDS
// THIS IS THE COMMANDS FOR NICO'S MAC
// Note first step you want path/to/cgal/CGAL-5.0/etc...

/Users/nandresthio/Documents/CGAL-5.0/scripts/cgal_create_CMakeLists -c Core -s augmentation



// Navigate to horizontalDisaster/cpp_code
cmake -DCMAKE_BUILD_TYPE=Release .
make

// RUN SPECS ARE HARDCODED IN augmentation.cpp
// To run, call
./augmentation