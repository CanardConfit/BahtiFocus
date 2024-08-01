using System;

namespace CanardConfit.NINA.BahtiFocus.Bahtinov {
    public static class LeastSquares
    {
        public static void BuildGeneralLeastSquares(
            ref double[] y, ref double[] weights, ref double[,] fMatrix, int numDataPoints, int numFunctions,
            ref double[] coefficients)
        {
            // Initialize matrices and variables for LQ decomposition and SVD
            var qMatrix = new double[0, 0];
            var ptMatrix = new double[0, 0];
            var tauLq = new double[0];
            var tauQ = new double[0];
            var tauP = new double[0];
            var singularValuesD = new double[0];
            var superDiagonalE = new double[0];
            var isUpper = false;
            var maxDataFunctions = Math.Max(numDataPoints, numFunctions);
            coefficients = new double[numFunctions];

            // Initialize matrices for least squares computation
            var leastSquaresMatrix = new double[numFunctions + 1, maxDataFunctions + 1];
            var rhsVector = new double[maxDataFunctions + 1];

            // Construct the right-hand side vector
            for (var i = 1; i <= numDataPoints; ++i)
                rhsVector[i] = weights[i - 1] * y[i - 1];
            for (var i = numDataPoints + 1; i <= numFunctions; ++i)
                rhsVector[i] = 0.0;

            // Construct the weighted function matrix
            for (var i = 1; i <= numFunctions; ++i)
            {
                for (var j = 1; j <= numDataPoints; ++j)
                    leastSquaresMatrix[i, j] = fMatrix[j - 1, i - 1];
            }

            // Apply weights to the function matrix
            for (var i = 1; i <= numFunctions; ++i)
            {
                for (var j = 1; j <= numDataPoints; ++j)
                    leastSquaresMatrix[i, j] *= weights[j - 1];
            }

            // Zero out the extra rows in the function matrix
            for (var i = 1; i <= numFunctions; ++i)
            {
                for (var j = numDataPoints + 1; j <= numFunctions; ++j)
                    leastSquaresMatrix[i, j] = 0.0;
            }

            // Perform LQ decomposition on the function matrix
            LQ.LQDecomposition(ref leastSquaresMatrix, numFunctions, maxDataFunctions, ref tauLq);
            LQ.UnPackQFromLQ(ref leastSquaresMatrix, numFunctions, maxDataFunctions, ref tauLq, numFunctions, ref qMatrix);

            // Solve the linear system using the Q matrix
            var rhsTransformed = new double[2, numFunctions + 1];
            for (var i = 1; i <= numFunctions; ++i)
                rhsTransformed[1, i] = 0.0;
            for (var i = 1; i <= numFunctions; ++i)
            {
                var sum = 0.0;
                for (var j = 1; j <= maxDataFunctions; ++j)
                    sum += rhsVector[j] * qMatrix[i, j];
                rhsTransformed[1, i] = sum;
            }

            // Transform the upper triangle of the function matrix
            for (var i = 1; i <= numFunctions - 1; ++i)
            {
                for (var j = i + 1; j <= numFunctions; ++j)
                    leastSquaresMatrix[i, j] = leastSquaresMatrix[j, i];
            }

            for (var i = 2; i <= numFunctions; ++i)
            {
                for (var j = 1; j <= i - 1; ++j)
                    leastSquaresMatrix[i, j] = 0.0;
            }

            // Convert to bidiagonal form and perform SVD
            bidiagonal.tobidiagonal(ref leastSquaresMatrix, numFunctions, numFunctions, ref tauQ, ref tauP);
            bidiagonal.multiplybyqfrombidiagonal(ref leastSquaresMatrix, numFunctions, numFunctions, ref tauQ,
                ref rhsTransformed, 1, numFunctions, true, false);
            bidiagonal.unpackptfrombidiagonal(ref leastSquaresMatrix, numFunctions, numFunctions, ref tauP, numFunctions,
                ref ptMatrix);
            bidiagonal.unpackdiagonalsfrombidiagonal(ref leastSquaresMatrix, numFunctions, numFunctions, ref isUpper,
                ref singularValuesD, ref superDiagonalE);

            if (!bdsvd.bidiagonalsvddecomposition(ref singularValuesD, superDiagonalE, numFunctions, isUpper, false,
                    ref rhsTransformed, 1, ref qMatrix, 0, ref ptMatrix, numFunctions))
            {
                for (var i = 0; i < numFunctions; ++i)
                    coefficients[i] = 0.0;
            }
            else
            {
                if (singularValuesD[1] != 0.0)
                {
                    for (var i = 1; i <= numFunctions; ++i)
                        rhsTransformed[1, i] = singularValuesD[i] <=
                                               5.0000000000000008E-15 * Math.Sqrt(numFunctions) * singularValuesD[1]
                            ? 0.0
                            : rhsTransformed[1, i] / singularValuesD[i];
                }

                for (var i = 1; i <= numFunctions; ++i)
                    rhsVector[i] = 0.0;
                for (var i = 1; i <= numFunctions; ++i)
                {
                    var value = rhsTransformed[1, i];
                    for (var j = 1; j <= numFunctions; ++j)
                        rhsVector[j] += value * ptMatrix[i, j];
                }

                for (var i = 0; i < numFunctions; ++i)
                    coefficients[i] = rhsVector[i + 1];
            }
        }

        public static void BuildLinearLeastSquares(ref double[] x, ref double[] y, int n, ref double slope,
            ref double intercept)
        {
            double numPoints = n;
            var sumXSquared = 0.0;
            var sumX = 0.0;
            var sumY = 0.0;
            var sumXY = 0.0;

            // Calculate sums needed for linear regression
            for (var i = 0; i < n; ++i)
            {
                sumX += x[i];
                sumXSquared += CustomMath.Sqr(x[i]);
                sumY += y[i];
                sumXY += x[i] * y[i];
            }

            // Calculate angle for rotation
            var angle = Math.Atan2(2.0 * sumX, sumXSquared - numPoints) / 2.0;
            var cosAngle = Math.Cos(angle);
            var sinAngle = Math.Sin(angle);

            // Calculate rotated sums
            var rotatedSum1 = CustomMath.Sqr(cosAngle) * numPoints + CustomMath.Sqr(sinAngle) * sumXSquared -
                              2.0 * sinAngle * cosAngle * sumX;
            var rotatedSum2 = CustomMath.Sqr(sinAngle) * numPoints + CustomMath.Sqr(cosAngle) * sumXSquared +
                              2.0 * sinAngle * cosAngle * sumX;
            var maxAbsRotatedSum = Math.Max(Math.Abs(rotatedSum1), Math.Abs(rotatedSum2));

            // Calculate rotated intercepts
            var rotatedIntercept1 = cosAngle * sumY - sinAngle * sumXY;
            var rotatedIntercept2 = sinAngle * sumY + cosAngle * sumXY;

            // Determine final slope and intercept
            var rotatedSlope1 = Math.Abs(rotatedSum1) <= maxAbsRotatedSum * 5E-16 * 1000.0
                ? 0.0
                : rotatedIntercept1 / rotatedSum1;
            var rotatedSlope2 = Math.Abs(rotatedSum2) <= maxAbsRotatedSum * 5E-16 * 1000.0
                ? 0.0
                : rotatedIntercept2 / rotatedSum2;

            slope = cosAngle * rotatedSlope1 + sinAngle * rotatedSlope2;
            intercept = -(sinAngle * rotatedSlope1) + cosAngle * rotatedSlope2;
        }

        public static void BuildSplineLeastSquares(
            ref double[] x, ref double[] y, ref double[] weights, double a, double b, int numDataPoints, int splineDegree,
            ref double[] splineCoefficients)
        {
            // Initialize matrices and variables for LQ decomposition and SVD
            var qMatrix = new double[0, 0];
            var ptMatrix = new double[0, 0];
            var tauLQ = new double[0];
            var tauQ = new double[0];
            var tauP = new double[0];
            var singularValuesD = new double[0];
            var superDiagonalE = new double[0];
            var isUpper = false;
            var numRows = numDataPoints;
            var numCols = splineDegree;
            var splineX = new double[numCols];
            var splineY = new double[numCols];
            var splineMatrix = new double[numCols + 1, Math.Max(numRows, numCols) + 1];
            var rhsVector = new double[Math.Max(numRows, numCols) + 1];

            // Create a set of equally spaced points for the spline
            for (var i = 0; i < numCols; ++i)
                splineX[i] = a + (b - a) * i / (numCols - 1);

            // Build the spline basis functions
            for (var i = 0; i < numCols; ++i)
            {
                // Initialize spline basis to zero
                for (var j = 0; j < numCols; ++j)
                    splineY[j] = 0.0;

                // Set the i-th basis function to 1
                splineY[i] = 1.0;

                // Build the cubic spline for the i-th basis function
                Spline3.BuildCubicSpline(splineX, splineY, numCols, 0, 0.0, 0, 0.0, ref splineCoefficients);

                // Evaluate the spline at each data point and build the matrix
                for (var j = 0; j < numRows; ++j)
                    splineMatrix[i + 1, j + 1] = weights[j] * Spline3.SplineInterpolation(ref splineCoefficients, x[j]);
            }

            // Zero out the extra rows in the spline matrix
            for (var i = 1; i <= numCols; ++i)
            {
                for (var j = numRows + 1; j <= numCols; ++j)
                    splineMatrix[i, j] = 0.0;
            }

            // Construct the right-hand side vector
            for (var i = 0; i < numRows; ++i)
                rhsVector[i + 1] = weights[i] * y[i];
            for (var i = numRows + 1; i <= numCols; ++i)
                rhsVector[i] = 0.0;

            var maxDim = Math.Max(numRows, numCols);

            // Perform LQ decomposition on the spline matrix
            LQ.LQDecomposition(ref splineMatrix, numCols, maxDim, ref tauLQ);
            LQ.UnPackQFromLQ(ref splineMatrix, numCols, maxDim, ref tauLQ, numCols, ref qMatrix);

            // Solve the linear system using the Q matrix
            var rhsTransformed = new double[2, numCols + 1];
            for (var i = 1; i <= numCols; ++i)
                rhsTransformed[1, i] = 0.0;
            for (var i = 1; i <= numCols; ++i)
            {
                var sum = 0.0;
                for (var j = 1; j <= maxDim; ++j)
                    sum += rhsVector[j] * qMatrix[i, j];
                rhsTransformed[1, i] = sum;
            }

            // Transform the upper triangle of the spline matrix
            for (var i = 1; i <= numCols - 1; ++i)
            {
                for (var j = i + 1; j <= numCols; ++j)
                    splineMatrix[i, j] = splineMatrix[j, i];
            }

            for (var i = 2; i <= numCols; ++i)
            {
                for (var j = 1; j <= i - 1; ++j)
                    splineMatrix[i, j] = 0.0;
            }

            // Convert to bidiagonal form and perform SVD
            bidiagonal.tobidiagonal(ref splineMatrix, numCols, numCols, ref tauQ, ref tauP);
            bidiagonal.multiplybyqfrombidiagonal(ref splineMatrix, numCols, numCols, ref tauQ, ref rhsTransformed, 1,
                numCols, true, false);
            bidiagonal.unpackptfrombidiagonal(ref splineMatrix, numCols, numCols, ref tauP, numCols, ref ptMatrix);
            bidiagonal.unpackdiagonalsfrombidiagonal(ref splineMatrix, numCols, numCols, ref isUpper, ref singularValuesD,
                ref superDiagonalE);

            if (!bdsvd.bidiagonalsvddecomposition(ref singularValuesD, superDiagonalE, numCols, isUpper, false,
                    ref rhsTransformed, 1, ref qMatrix, 0, ref ptMatrix, numCols))
            {
                for (var i = 1; i <= numCols; ++i)
                {
                    singularValuesD[i] = 0.0;
                    rhsTransformed[1, i] = 0.0;
                    for (var j = 1; j <= numCols; ++j)
                        ptMatrix[i, j] = i != j ? 0.0 : 1.0;
                }

                rhsTransformed[1, 1] = 1.0;
            }

            // Solve the linear system for the coefficients
            for (var i = 1; i <= numCols; ++i)
                rhsTransformed[1, i] =
                    singularValuesD[i] <= 5.0000000000000008E-15 * Math.Sqrt(numCols) * singularValuesD[1]
                        ? 0.0
                        : rhsTransformed[1, i] / singularValuesD[i];
            for (var i = 1; i <= numCols; ++i)
                rhsVector[i] = 0.0;
            for (var i = 1; i <= numCols; ++i)
            {
                var value = rhsTransformed[1, i];
                for (var j = 1; j <= numCols; ++j)
                    rhsVector[j] += value * ptMatrix[i, j];
            }

            for (var i = 0; i < numCols; ++i)
                splineY[i] = rhsVector[i + 1];

            // Build the final cubic spline with the calculated coefficients
            Spline3.BuildCubicSpline(splineX, splineY, numCols, 0, 0.0, 0, 0.0, ref splineCoefficients);
        }

        public static void BuildPolynomialLeastSquares(ref double[] x, ref double[] y, int numDataPoints,
            int polynomialDegree, ref double[] coefficients)
        {
            var chebyshevCoefficients = new double[0];
            var maxX = x[0];
            var minX = x[0];

            // Find the minimum and maximum values of x
            for (var i = 1; i < numDataPoints; ++i)
            {
                if (x[i] > maxX)
                    maxX = x[i];
                if (x[i] < minX)
                    minX = x[i];
            }

            // Adjust minX and maxX if they are equal
            if (minX == maxX)
            {
                minX -= 0.5;
                maxX += 0.5;
            }

            // Initialize weights to 1.0 for all data points
            var weights = new double[numDataPoints];
            for (var i = 0; i < numDataPoints; ++i)
                weights[i] = 1.0;

            // Build Chebyshev least squares coefficients
            BuildChebyshevLeastSquares(ref x, ref y, ref weights, minX, maxX, numDataPoints, polynomialDegree,
                ref chebyshevCoefficients);

            var polynomialCoefficients = new double[polynomialDegree + 1];
            for (var i = 0; i <= polynomialDegree; ++i)
                polynomialCoefficients[i] = 0.0;

            var temp = 0.0;
            for (var i = 0; i <= polynomialDegree; ++i)
            {
                for (var j = i; j <= polynomialDegree; ++j)
                {
                    var previousTemp = polynomialCoefficients[j];
                    polynomialCoefficients[j] = 0.0;
                    if ((i <= 1) && (j == i))
                    {
                        polynomialCoefficients[j] = 1.0;
                    }
                    else
                    {
                        if (i != 0)
                            polynomialCoefficients[j] = 2.0 * temp;
                        if (j > i + 1)
                            polynomialCoefficients[j] -= polynomialCoefficients[j - 2];
                    }

                    temp = previousTemp;
                }

                temp = polynomialCoefficients[i];
                var sum = 0.0;
                for (var j = i; j <= polynomialDegree; j += 2)
                    sum += polynomialCoefficients[j] * chebyshevCoefficients[j];
                polynomialCoefficients[i] = sum;
            }

            var scale = 2.0 / (chebyshevCoefficients[polynomialDegree + 2] - chebyshevCoefficients[polynomialDegree + 1]);
            var shift = -(2.0 * chebyshevCoefficients[polynomialDegree + 1] / (chebyshevCoefficients[polynomialDegree + 2] -
                                                                               chebyshevCoefficients
                                                                                   [polynomialDegree + 1])) - 1.0;

            coefficients = new double[polynomialDegree + 1];
            var tempArray1 = new double[polynomialDegree + 1];
            var tempArray2 = new double[polynomialDegree + 1];
            coefficients[0] = polynomialCoefficients[0];
            tempArray2[0] = 1.0;
            tempArray1[0] = 1.0;

            for (var i = 1; i <= polynomialDegree; ++i)
            {
                tempArray1[i] = 1.0;
                tempArray2[i] = shift * tempArray2[i - 1];
                coefficients[0] += polynomialCoefficients[i] * tempArray2[i];
            }

            for (var i = 1; i <= polynomialDegree; ++i)
            {
                tempArray1[0] *= scale;
                coefficients[i] = polynomialCoefficients[i] * tempArray1[0];
                for (var j = i + 1; j <= polynomialDegree; ++j)
                {
                    var k = j - i;
                    tempArray1[k] = scale * tempArray1[k] + tempArray1[k - 1];
                    coefficients[i] += polynomialCoefficients[j] * tempArray1[k] * tempArray2[k];
                }
            }
        }


        public static void BuildChebyshevLeastSquares(ref double[] x, ref double[] y, ref double[] weights, double a,
            double b, int numPoints, int degree, ref double[] chebyshevCoefficients)
        {
            // Initialize matrices and variables for LQ decomposition and SVD
            var qMatrix = new double[0, 0];
            var ptMatrix = new double[0, 0];
            var tauLq = new double[0];
            var tauQ = new double[0];
            var tauP = new double[0];
            var singularValuesD = new double[0];
            var superDiagonalE = new double[0];
            var isUpper = false;
            var numRows = numPoints;
            var numCols = degree + 1;
            var chebyshevMatrix = new double[numCols + 1, Math.Max(numRows, numCols) + 1];
            var rhsVector = new double[Math.Max(numRows, numCols) + 1];

            // Construct the Chebyshev matrix
            for (var col = 1; col <= numCols; ++col)
            {
                for (var row = 1; row <= numRows; ++row)
                {
                    var normalizedX = 2.0 * (x[row - 1] - a) / (b - a) - 1.0;
                    if (col == 1)
                        chebyshevMatrix[col, row] = 1.0;
                    if (col == 2)
                        chebyshevMatrix[col, row] = normalizedX;
                    if (col > 2)
                        chebyshevMatrix[col, row] = 2.0 * normalizedX * chebyshevMatrix[col - 1, row] -
                                                    chebyshevMatrix[col - 2, row];
                }
            }

            // Apply weights to the Chebyshev matrix
            for (var col = 1; col <= numCols; ++col)
            {
                for (var row = 1; row <= numRows; ++row)
                    chebyshevMatrix[col, row] = weights[row - 1] * chebyshevMatrix[col, row];
            }

            // Zero out the extra rows in the Chebyshev matrix
            for (var col = 1; col <= numCols; ++col)
            {
                for (var row = numRows + 1; row <= numCols; ++row)
                    chebyshevMatrix[col, row] = 0.0;
            }

            // Construct the right-hand side vector
            for (var i = 0; i < numRows; ++i)
                rhsVector[i + 1] = weights[i] * y[i];
            for (var i = numRows + 1; i <= numCols; ++i)
                rhsVector[i] = 0.0;

            // Perform LQ decomposition on the Chebyshev matrix
            var maxDim = Math.Max(numRows, numCols);
            LQ.LQDecomposition(ref chebyshevMatrix, numCols, maxDim, ref tauLq);
            LQ.UnPackQFromLQ(ref chebyshevMatrix, numCols, maxDim, ref tauLq, numCols, ref qMatrix);

            // Solve the linear system using the Q matrix
            var rhsTransformed = new double[2, numCols + 1];
            for (var col = 1; col <= numCols; ++col)
                rhsTransformed[1, col] = 0.0;
            for (var col = 1; col <= numCols; ++col)
            {
                var sum = 0.0;
                for (var row = 1; row <= maxDim; ++row)
                    sum += rhsVector[row] * qMatrix[col, row];
                rhsTransformed[1, col] = sum;
            }

            // Transform the upper triangle of the Chebyshev matrix
            for (var row = 1; row <= numCols - 1; ++row)
            {
                for (var col = row + 1; col <= numCols; ++col)
                    chebyshevMatrix[row, col] = chebyshevMatrix[col, row];
            }

            for (var row = 2; row <= numCols; ++row)
            {
                for (var col = 1; col <= row - 1; ++col)
                    chebyshevMatrix[row, col] = 0.0;
            }

            // Convert to bidiagonal form and perform SVD
            bidiagonal.tobidiagonal(ref chebyshevMatrix, numCols, numCols, ref tauQ, ref tauP);
            bidiagonal.multiplybyqfrombidiagonal(ref chebyshevMatrix, numCols, numCols, ref tauQ, ref rhsTransformed, 1,
                numCols, true, false);
            bidiagonal.unpackptfrombidiagonal(ref chebyshevMatrix, numCols, numCols, ref tauP, numCols, ref ptMatrix);
            bidiagonal.unpackdiagonalsfrombidiagonal(ref chebyshevMatrix, numCols, numCols, ref isUpper,
                ref singularValuesD, ref superDiagonalE);

            if (!bdsvd.bidiagonalsvddecomposition(ref singularValuesD, superDiagonalE, numCols, isUpper, false,
                    ref rhsTransformed, 1, ref qMatrix, 0, ref ptMatrix, numCols))
            {
                for (var i = 1; i <= numCols; ++i)
                {
                    singularValuesD[i] = 0.0;
                    rhsTransformed[1, i] = 0.0;
                    for (var j = 1; j <= numCols; ++j)
                        ptMatrix[i, j] = i != j ? 0.0 : 1.0;
                }

                rhsTransformed[1, 1] = 1.0;
            }

            // Solve the linear system for the coefficients
            for (var i = 1; i <= numCols; ++i)
                rhsTransformed[1, i] =
                    singularValuesD[i] <= 5.0000000000000008E-15 * Math.Sqrt(numCols) * singularValuesD[1]
                        ? 0.0
                        : rhsTransformed[1, i] / singularValuesD[i];
            for (var i = 1; i <= numCols; ++i)
                rhsVector[i] = 0.0;
            for (var i = 1; i <= numCols; ++i)
            {
                var value = rhsTransformed[1, i];
                for (var j = 1; j <= numCols; ++j)
                    rhsVector[j] += value * ptMatrix[i, j];
            }

            // Store the resulting coefficients and boundaries
            chebyshevCoefficients = new double[numCols + 1 + 1];
            for (var i = 1; i <= numCols; ++i)
                chebyshevCoefficients[i - 1] = rhsVector[i];
            chebyshevCoefficients[numCols] = a;
            chebyshevCoefficients[numCols + 1] = b;
        }


        public static bool BuildChebyshevLeastSquaresConstrained(ref double[] x, ref double[] y, ref double[] weights,
            double a,
            double b, int numPoints, ref double[] constrainedX, ref double[] constrainedY, ref int[] constraintTypes,
            int numConstraints,
            int degree, ref double[] chebyshevCoefficients)
        {
            // Initialize matrices and variables for SVD and least squares computation
            var uMatrix = new double[0, 0];
            var vtMatrix = new double[0, 0];
            var singularValues = new double[0];
            var successFlag = true;
            var chebyshevMatrix = new double[Math.Max(numPoints, degree + 1) + 1, degree + 1 + 1];
            var rhsVector = new double[Math.Max(numPoints, degree + 1) + 1];

            // Construct the Chebyshev matrix
            for (var i = 1; i <= numPoints; ++i)
            {
                for (var j = 1; j <= degree + 1; ++j)
                {
                    var normalizedX = 2.0 * (x[i - 1] - a) / (b - a) - 1.0;
                    if (j == 1)
                        chebyshevMatrix[i, j] = 1.0;
                    if (j == 2)
                        chebyshevMatrix[i, j] = normalizedX;
                    if (j > 2)
                        chebyshevMatrix[i, j] = 2.0 * normalizedX * chebyshevMatrix[i, j - 1] - chebyshevMatrix[i, j - 2];
                }
            }

            // Apply weights to the Chebyshev matrix
            for (var i = 1; i <= numPoints; ++i)
            {
                for (var j = 1; j <= degree + 1; ++j)
                    chebyshevMatrix[i, j] = weights[i - 1] * chebyshevMatrix[i, j];
            }

            // Zero out the extra rows in the Chebyshev matrix
            for (var i = numPoints + 1; i <= degree + 1; ++i)
            {
                for (var j = 1; j <= degree + 1; ++j)
                    chebyshevMatrix[i, j] = 0.0;
            }

            // Construct the right-hand side vector
            for (var i = 0; i < numPoints; ++i)
                rhsVector[i + 1] = weights[i] * y[i];
            for (var i = numPoints + 1; i <= degree + 1; ++i)
                rhsVector[i] = 0.0;

            numPoints = Math.Max(numPoints, degree + 1);

            // Initialize constraint matrices and vectors
            var constraintMatrix = new double[degree + 1 + 1, degree + 1 + 1];
            var constrainedRhs = new double[degree + 1 + 1];

            if (numConstraints == 0)
            {
                // If there are no constraints, set the constraint matrix to identity
                for (var i = 1; i <= degree + 1; ++i)
                {
                    for (var j = 1; j <= degree + 1; ++j)
                        constraintMatrix[i, j] = 0.0;
                    constrainedRhs[i] = 0.0;
                }

                for (var i = 1; i <= degree + 1; ++i)
                    constraintMatrix[i, i] = 1.0;

                numConstraints = degree + 1;
            }
            else
            {
                // Construct the constraint matrix based on the constraints
                var constraintFunctionMatrix = new double[numConstraints + 1, degree + 1 + 1];
                var constraintValues = new double[numConstraints + 1];
                var chebyshevValues = new double[degree + 1];
                var chebyshevDerivatives = new double[degree + 1];
                var chebyshevSecondDerivatives = new double[degree + 1];

                for (var i = 0; i < numConstraints; ++i)
                {
                    var normalizedX = 2.0 * (constrainedX[i] - a) / (b - a) - 1.0;
                    for (var j = 0; j <= degree; ++j)
                    {
                        if (j == 0)
                        {
                            chebyshevValues[j] = 1.0;
                            chebyshevDerivatives[j] = 1.0;
                            chebyshevSecondDerivatives[j] = 0.0;
                        }

                        if (j == 1)
                        {
                            chebyshevValues[j] = normalizedX;
                            chebyshevDerivatives[j] = 2.0 * normalizedX;
                            chebyshevSecondDerivatives[j] = 1.0;
                        }

                        if (j > 1)
                        {
                            chebyshevValues[j] = 2.0 * normalizedX * chebyshevValues[j - 1] - chebyshevValues[j - 2];
                            chebyshevDerivatives[j] = 2.0 * normalizedX * chebyshevDerivatives[j - 1] -
                                                      chebyshevDerivatives[j - 2];
                            chebyshevSecondDerivatives[j] = j * chebyshevDerivatives[j - 1];
                        }

                        if (constraintTypes[i] == 0)
                            constraintFunctionMatrix[i + 1, j + 1] = chebyshevValues[j];
                        if (constraintTypes[i] == 1)
                            constraintFunctionMatrix[i + 1, j + 1] = chebyshevSecondDerivatives[j];
                    }

                    constraintValues[i + 1] = constrainedY[i];
                }

                // Perform SVD decomposition on the constraint matrix
                if (!SVD.SVDDecomposition(constraintFunctionMatrix, numConstraints, degree + 1, 2, 2, 2, ref singularValues,
                        ref uMatrix, ref vtMatrix) ||
                    singularValues[1] == 0.0 || singularValues[numConstraints] <=
                    5.0000000000000008E-15 * Math.Sqrt(numConstraints) * singularValues[1])
                    return false;

                // Construct the unconstrained solution
                constraintMatrix = new double[degree + 1 + 1, degree + 1 - numConstraints + 1];
                constrainedRhs = new double[degree + 1 + 1];

                for (var i = 1; i <= degree + 1 - numConstraints; ++i)
                {
                    for (var j = 1; j <= degree + 1; ++j)
                        constraintMatrix[j, i] = vtMatrix[numConstraints + i, j];
                }

                var singularValuesInverse = new double[numConstraints + 1];
                for (var i = 1; i <= numConstraints; ++i)
                {
                    var sum = 0.0;
                    for (var j = 1; j <= numConstraints; ++j)
                        sum += uMatrix[j, i] * constraintValues[j];
                    singularValuesInverse[i] = sum / singularValues[i];
                }

                for (var i = 1; i <= degree + 1; ++i)
                    constrainedRhs[i] = 0.0;
                for (var i = 1; i <= numConstraints; ++i)
                {
                    var value = singularValuesInverse[i];
                    for (var j = 1; j <= degree + 1; ++j)
                        constrainedRhs[j] += value * vtMatrix[i, j];
                }

                for (var i = 1; i <= numPoints; ++i)
                {
                    var sum = 0.0;
                    for (var j = 1; j <= degree + 1; ++j)
                        sum += chebyshevMatrix[i, j] * constrainedRhs[j];
                    rhsVector[i] -= sum;
                }

                numConstraints = degree + 1 - numConstraints;
                var constraintAdjustedMatrix = new double[numPoints + 1, numConstraints + 1];
                var work = new double[numPoints + 1];
                blas.matrixmatrixmultiply(ref chebyshevMatrix, 1, numPoints, 1, degree + 1, false, ref constraintMatrix, 1,
                    degree + 1, 1, numConstraints, false, 1.0, ref constraintAdjustedMatrix, 1, numPoints, 1,
                    numConstraints, 0.0, ref work);
                blas.copymatrix(ref constraintAdjustedMatrix, 1, numPoints, 1, numConstraints, ref chebyshevMatrix, 1,
                    numPoints, 1, numConstraints);
            }

            // Perform SVD decomposition on the adjusted Chebyshev matrix
            if (!SVD.SVDDecomposition(chebyshevMatrix, numPoints, numConstraints, 1, 1, 2, ref singularValues, ref uMatrix,
                    ref vtMatrix))
                return false;

            var workVector = new double[numConstraints + 1];
            var coefficients = new double[numConstraints + 1];

            for (var i = 1; i <= numConstraints; ++i)
                workVector[i] = 0.0;
            for (var i = 1; i <= numPoints; ++i)
            {
                var value = rhsVector[i];
                for (var j = 1; j <= numConstraints; ++j)
                    workVector[j] += value * uMatrix[i, j];
            }

            for (var i = 1; i <= numConstraints; ++i)
                workVector[i] = !(singularValues[i] != 0.0 &&
                                  singularValues[i] >
                                  5.0000000000000008E-15 * Math.Sqrt(numConstraints) * singularValues[1])
                    ? 0.0
                    : workVector[i] / singularValues[i];
            for (var i = 1; i <= numConstraints; ++i)
                coefficients[i] = 0.0;
            for (var i = 1; i <= numConstraints; ++i)
            {
                var value = workVector[i];
                for (var j = 1; j <= numConstraints; ++j)
                    coefficients[j] += value * vtMatrix[i, j];
            }

            chebyshevCoefficients = new double[degree + 2 + 1];
            for (var i = 1; i <= degree + 1; ++i)
            {
                var sum = 0.0;
                for (var j = 1; j <= numConstraints; ++j)
                    sum += constraintMatrix[i, j] * coefficients[j];
                chebyshevCoefficients[i - 1] = sum + constrainedRhs[i];
            }

            chebyshevCoefficients[degree + 1] = a;
            chebyshevCoefficients[degree + 2] = b;

            return successFlag;
        }


        public static double CalculateChebyshevLeastSquares(int degree, ref double[] coefficients, double x)
        {
            // Normalise x to the range [-1, 1]
            x = 2.0 * (x - coefficients[degree + 1]) / (coefficients[degree + 2] - coefficients[degree + 1]) - 1.0;

            // Initialize variables for the recurrence relation
            var previousTerm = 0.0;
            var currentTerm = 0.0;
            var index = degree;
            double nextTerm;

            // Apply the Chebyshev recurrence relation
            do
            {
                nextTerm = 2.0 * x * currentTerm - previousTerm + coefficients[index];
                previousTerm = currentTerm;
                currentTerm = nextTerm;
                --index;
            } while (index >= 0);

            // Return the calculated Chebyshev least squares value
            return nextTerm - x * previousTerm;
        }
    }
}