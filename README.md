# FTTC
Simulation of chemical dynamics in high dimensionality with low-rank functional tensor-train Chebyshev (FTTC) quantum dynamics 

# Dependencies

The functional tensor-train Chebyshev (FTTC) program requires the following:

- C compiler
- CMake
- Compressed Continuous Computation (C3) library (which requires BLAS, LAPACK, and optionally SWIG-python)

The tensor-train (TT) Chebyshev program requires the following:

- MATLAB
- H-Tucker Tensor Toolbox
- TT-Toolbox

Comparison of FTTC and TT results with the provided scripts requires:

- Gnuplot

# Instructions

## FTTC

To run the FTTC code in folder ``FunctionTrain``, compile with commands:
```
cmake .
make
```
and run the program as follows:
```
./ft_cheb
```
Troubleshooting details are included in the ``readme.txt`` file in the FT-Chebyshev folder. 

## TT-Chebyhsev

To run the TT-Chebyshev code in folder ``TensorTrain``, run ``tt_CHEBSOFT.m`` in MATLAB. Data files can then be converted to Gnuplot format by running ``GnuplotConversion.m`` in MATLAB.

## Comparison Code

FTTC and TT-Chebyshev results can be compared by moving the data files to folder ``TensorTrainData`` and ``FunctionTrainData`` in folder ``Comparisons`` and using the following commands:

```
gnuplot wavefunctionsnapshots.gpt
gnuplot autocorrelationcompare.gpt
```

# Additional Documentation

Additional documentation is included in the ``CodeDocumentation.pdf`` file in the main folder.
