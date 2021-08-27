# TTandFTChebyshev
Simulation of chemical dynamics in high dimensionality with low-rank tensor-train Chebyshev (TT-Chebyshev) and function-train Chebyshev (FT-Chebyshev) quantum dynamics 

# Dependencies

The TT-Chebyshev program requires the following:

- MATLAB
- H-Tucker Tensor Toolbox
- TT-Toolbox

The FT-Chebyshev program requires the following:

- C compiler
- CMake
- Compressed Continuous Computation (C3) library (which requires BLAS, LAPACK, and optionally SWIG-python)

Comparison of TT- and FT-Chebyshev results with the provided scripts requires:

- Gnuplot

# Instructions

To run the TT-Chebyshev code in folder ``TensorTrain``, run ``tt_CHEBSOFT.m`` in MATLAB. Data files can then be converted to Gnuplot format by running ``GnuplotConversion.m`` in MATLAB.

To run the FT-Chebyshev code in folder ``FunctionTrain``, compile with commands:
``
cmake .
make
``
and run the program as follows:
``
./ft_cheb
``
Troubleshooting details are included in the ``readme.txt`` file in the FT-Chebyshev folder. 

TT- and FT-Chebyshev results can be compared by moving the TT- and FT-Chebyshev data files to folder ``TensorTrainData`` and ``FunctionTrainData`` in folder ``Comparisons`` and using the following commands:

``
gnuplot wavefunctionsnapshots.gpt
gnuplot autocorrelationcompare.gpt
``

Additional documentation is included in the ``CodeDocumentation.pdf`` file in the main folder.
