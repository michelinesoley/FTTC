To get c3:

1. Download BLAS, LAPACK, and SWIG-python as follows:
brew install openblas
brew install lapack
sudo port install swig-python

2. Install c3 as follows: 
git clone https://github.com/goroda/Compressed-Continuous-Computation.git c3
cd c3
mkdir build
cd build
cmake ..
sudo make install
(Note the final line differs from the suggested installation instructions. Using sudo make install will facilitate linking of c3 but centrally change your c3 code. Check that any c3-based python codes function with the newly installed version of c3.)
(Note you don't need to reclone c3 if you already have an old version, just do "git pull" in the c3 directory to get the new code.)

To get c3py working, it is helpful to create a conda environment:
Write:
conda create -n c3env python=3.7 numpy scipy
conda activate c3env

Then enter the folder c3:
python setup.py build
python setup.py install
pip install matplotlib

Then one can open the environment again by first entering:
conda activate c3env

To run:
python fttc.py

Troubleshooting:

If Accelerate is not found but Accelerate.tbd exists in any folder, copy the Accelerate.tbd file into the same folder and rename it Accelerate.

If make yields an error that c3/array.h is missing, -lc3 cannot be found, etc., try the following:
1. Ensure that c3 is downloaded. Then update the version of cmake as follows:
cmake_minimum_required (VERSION 3.0)
2. Change the final argument of target_link_libraries to the following:
${c3lib} 
3. Update the location of the c3 lib files in CMakeLists.txt inside the if APPLE loop as follows:
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/lib)
4. Update the location of the c3 include and lib files outside of the if APPLE loop as follows:
find_library(c3lib c3 PATHS “/usr/local/lib”)
include_directories(/usr/local/include)
