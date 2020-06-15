Stabilised Finite Element Method for incompressible Navier-Stokes. This code uses the Petsc library for solving the matrix system. The formulation is presented in https://arxiv.org/abs/2001.04925


### Build and compile

    * modify the CMakeLists.txt with paths to the libraries on your system
    * Create `build` and `bin` folders
    * `mkdir build`
    * `mkdir bin`
    * `cd build`
    * `cmake ..`
    * If you get errors, then fix them and rerun `cmake ..`
    * After configuring successfully, compile and build the executable
    * `make install`
    * If you get compilation errors, then fix them and run `make install` untill successful build and installation of the executables

Once successfully compiled, the executable `cfdsfemmpi` and the executables for the tests are copied to the `bin` folder

### Tests
* Test cases are provided as separate executables. All the inputs are hardcoded in the respective source file. So, the test cases do not require any input arguments. See the screen output to check if the tests are executed successfully or not.

* To run the test cases with 4 processors
    * mpirun -n 4 ./test-ldc-stru-nelem10
    * mpirun -n 4 ./test-ldc-stru-mesh2

### Usage with other input files

* To run the simulations, from `bin` directory
    * `mpirun -n <number-of-processors>  ./cfdsfemmpi  <input-file-prefix>  <control-parameters-file-name>`<petsc-options-file-name>
    * Example,
    * `mpirun -n 10 cfdsfemmpi  LDC2Dquad-stru-FI  control-parameters-Re100-dt0p1.dat petsc_options.dat`
