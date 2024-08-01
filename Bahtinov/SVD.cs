using System;

namespace CanardConfit.NINA.BahtiFocus.Bahtinov {
    public static class SVD {
        public static bool RMatrixSVD(double[,] inputMatrix, int numRows, int numCols, int uNeeded, int vtNeeded,
            int additionalMemory, ref double[] singularValues, ref double[,] leftSingularVectors,
            ref double[,] rightSingularVectors) {
            var tauQ = new double[0];
            var tauP = new double[0];
            var tau = new double[0];
            var superDiagonal = new double[0];
            var tempMatrix = new double[0, 0];
            var isUpperBidiagonal = false;

            // Clone the input matrix to avoid modifying the original
            inputMatrix = (double[,])inputMatrix.Clone();

            // Check for zero dimensions
            if (numRows == 0 || numCols == 0)
                return true;

            // Determine the smaller dimension between numRows and numCols
            var minDim = Math.Min(numRows, numCols);

            // Initialize singular values array
            singularValues = new double[minDim + 1];

            // Initialize matrices U and VT based on the needed flags
            int uRowCount = 0, uColCount = 0;
            if (uNeeded == 1) {
                uRowCount = numRows;
                uColCount = minDim;
                leftSingularVectors = new double[uRowCount, uColCount];
            } else if (uNeeded == 2) {
                uRowCount = numRows;
                uColCount = numRows;
                leftSingularVectors = new double[uRowCount, uColCount];
            }

            int vtRowCount = 0, vtColCount = 0;
            if (vtNeeded == 1) {
                vtRowCount = minDim;
                vtColCount = numCols;
                rightSingularVectors = new double[vtRowCount, vtColCount];
            } else if (vtNeeded == 2) {
                vtRowCount = numCols;
                vtColCount = numCols;
                rightSingularVectors = new double[vtRowCount, vtColCount];
            }

            // Case when numRows is significantly larger than numCols
            if (numRows > 1.6 * numCols) {
                // Perform QR decomposition on the input matrix
                QR.RMatrixQR(ref inputMatrix, numRows, numCols, ref tau);

                if (uNeeded == 0) {
                    // Transform to bidiagonal form
                    for (int i = 0; i < numCols; ++i)
                    for (int j = 0; j < i; ++j)
                        inputMatrix[i, j] = 0.0;

                    bidiagonal.rmatrixbd(ref inputMatrix, numCols, numCols, ref tauQ, ref tauP);
                    bidiagonal.rmatrixbdunpackpt(ref inputMatrix, numCols, numCols, ref tauP, vtRowCount,
                        ref rightSingularVectors);
                    bidiagonal.rmatrixbdunpackdiagonals(ref inputMatrix, numCols, numCols, ref isUpperBidiagonal,
                        ref singularValues, ref superDiagonal);

                    return bdsvd.rmatrixbdsvd(ref singularValues, superDiagonal, numCols, isUpperBidiagonal, false,
                        ref leftSingularVectors, 0, ref inputMatrix, 0, ref rightSingularVectors, vtColCount);
                }

                // Unpack Q from the QR decomposition
                QR.RMatrixQRUnPackQ(ref inputMatrix, numRows, numCols, ref tau, uColCount, ref leftSingularVectors);

                // Transform to bidiagonal form
                for (int i = 0; i < numCols; ++i)
                for (int j = 0; j < i; ++j)
                    inputMatrix[i, j] = 0.0;

                bidiagonal.rmatrixbd(ref inputMatrix, numCols, numCols, ref tauQ, ref tauP);
                bidiagonal.rmatrixbdunpackpt(ref inputMatrix, numCols, numCols, ref tauP, vtRowCount,
                    ref rightSingularVectors);
                bidiagonal.rmatrixbdunpackdiagonals(ref inputMatrix, numCols, numCols, ref isUpperBidiagonal,
                    ref singularValues, ref superDiagonal);

                bool flag1;
                if (additionalMemory < 1) {
                    // Multiply by Q to form U
                    bidiagonal.rmatrixbdmultiplybyq(ref inputMatrix, numCols, numCols, ref tauQ,
                        ref leftSingularVectors, numRows, numCols, true, false);
                    flag1 = bdsvd.rmatrixbdsvd(ref singularValues, superDiagonal, numCols, isUpperBidiagonal, false,
                        ref leftSingularVectors, numRows, ref inputMatrix, 0, ref rightSingularVectors, vtColCount);
                } else {
                    var workArray = new double[Math.Max(numRows, numCols) + 1];
                    bidiagonal.rmatrixbdunpackq(ref inputMatrix, numCols, numCols, ref tauQ, numCols, ref tempMatrix);
                    blas.copymatrix(ref leftSingularVectors, 0, numRows - 1, 0, numCols - 1, ref inputMatrix, 0,
                        numRows - 1, 0, numCols - 1);
                    blas.inplacetranspose(ref tempMatrix, 0, numCols - 1, 0, numCols - 1, ref workArray);
                    flag1 = bdsvd.rmatrixbdsvd(ref singularValues, superDiagonal, numCols, isUpperBidiagonal, false,
                        ref leftSingularVectors, 0, ref tempMatrix, numCols, ref rightSingularVectors, vtColCount);
                    blas.matrixmatrixmultiply(ref inputMatrix, 0, numRows - 1, 0, numCols - 1, false, ref tempMatrix, 0,
                        numCols - 1, 0, numCols - 1, true, 1.0, ref leftSingularVectors, 0, numRows - 1, 0, numCols - 1,
                        0.0, ref workArray);
                }

                return flag1;
            }

            // Case when numCols is significantly larger than numRows
            if (numCols > 1.6 * numRows) {
                // Perform LQ decomposition on the input matrix
                LQ.RMatrixLQ(ref inputMatrix, numRows, numCols, ref tau);

                if (vtNeeded == 0) {
                    // Transform to bidiagonal form
                    for (int i = 0; i < numRows; ++i)
                    for (int j = i + 1; j < numRows; ++j)
                        inputMatrix[i, j] = 0.0;

                    bidiagonal.rmatrixbd(ref inputMatrix, numRows, numRows, ref tauQ, ref tauP);
                    bidiagonal.rmatrixbdunpackq(ref inputMatrix, numRows, numRows, ref tauQ, uColCount,
                        ref leftSingularVectors);
                    bidiagonal.rmatrixbdunpackdiagonals(ref inputMatrix, numRows, numRows, ref isUpperBidiagonal,
                        ref singularValues, ref superDiagonal);

                    var workArray = new double[numRows + 1];
                    blas.inplacetranspose(ref leftSingularVectors, 0, uRowCount - 1, 0, uColCount - 1, ref workArray);
                    var flag2 = bdsvd.rmatrixbdsvd(ref singularValues, superDiagonal, numRows, isUpperBidiagonal, false,
                        ref inputMatrix, 0, ref leftSingularVectors, uRowCount, ref rightSingularVectors, 0);
                    blas.inplacetranspose(ref leftSingularVectors, 0, uRowCount - 1, 0, uColCount - 1, ref workArray);
                    return flag2;
                }

                // Unpack Q from the LQ decomposition
                LQ.RMatrixLQUnPackQ(ref inputMatrix, numRows, numCols, ref tau, vtRowCount, ref rightSingularVectors);

                // Transform to bidiagonal form
                for (int i = 0; i < numRows; ++i)
                for (int j = i + 1; j < numRows; ++j)
                    inputMatrix[i, j] = 0.0;

                bidiagonal.rmatrixbd(ref inputMatrix, numRows, numRows, ref tauQ, ref tauP);
                bidiagonal.rmatrixbdunpackq(ref inputMatrix, numRows, numRows, ref tauQ, uColCount,
                    ref leftSingularVectors);
                bidiagonal.rmatrixbdunpackdiagonals(ref inputMatrix, numRows, numRows, ref isUpperBidiagonal,
                    ref singularValues, ref superDiagonal);

                var workArray1 = new double[Math.Max(numRows, numCols) + 1];
                blas.inplacetranspose(ref leftSingularVectors, 0, uRowCount - 1, 0, uColCount - 1, ref workArray1);
                bool flag3;
                if (additionalMemory < 1) {
                    bidiagonal.rmatrixbdmultiplybyp(ref inputMatrix, numRows, numRows, ref tauP,
                        ref rightSingularVectors, numRows, numCols, false, true);
                    flag3 = bdsvd.rmatrixbdsvd(ref singularValues, superDiagonal, numRows, isUpperBidiagonal, false,
                        ref inputMatrix, 0, ref leftSingularVectors, uRowCount, ref rightSingularVectors, numCols);
                } else {
                    bidiagonal.rmatrixbdunpackpt(ref inputMatrix, numRows, numRows, ref tauP, numRows, ref tempMatrix);
                    flag3 = bdsvd.rmatrixbdsvd(ref singularValues, superDiagonal, numRows, isUpperBidiagonal, false,
                        ref inputMatrix, 0, ref leftSingularVectors, uRowCount, ref tempMatrix, numRows);
                    blas.copymatrix(ref rightSingularVectors, 0, numRows - 1, 0, numCols - 1, ref inputMatrix, 0,
                        numRows - 1, 0, numCols - 1);
                    blas.matrixmatrixmultiply(ref tempMatrix, 0, numRows - 1, 0, numRows - 1, false, ref inputMatrix, 0,
                        numRows - 1, 0, numCols - 1, false, 1.0, ref rightSingularVectors, 0, numRows - 1, 0,
                        numCols - 1, 0.0, ref workArray1);
                }

                blas.inplacetranspose(ref leftSingularVectors, 0, uRowCount - 1, 0, uColCount - 1, ref workArray1);
                return flag3;
            }

            // General case for other dimensions
            if (numRows <= numCols) {
                bidiagonal.rmatrixbd(ref inputMatrix, numRows, numCols, ref tauQ, ref tauP);
                bidiagonal.rmatrixbdunpackq(ref inputMatrix, numRows, numCols, ref tauQ, uColCount,
                    ref leftSingularVectors);
                bidiagonal.rmatrixbdunpackpt(ref inputMatrix, numRows, numCols, ref tauP, vtRowCount,
                    ref rightSingularVectors);
                bidiagonal.rmatrixbdunpackdiagonals(ref inputMatrix, numRows, numCols, ref isUpperBidiagonal,
                    ref singularValues, ref superDiagonal);

                var workArray = new double[numRows + 1];
                blas.inplacetranspose(ref leftSingularVectors, 0, uRowCount - 1, 0, uColCount - 1, ref workArray);
                var success = bdsvd.rmatrixbdsvd(ref singularValues, superDiagonal, minDim, isUpperBidiagonal, false,
                    ref inputMatrix, 0, ref leftSingularVectors, uRowCount, ref rightSingularVectors, vtColCount);
                blas.inplacetranspose(ref leftSingularVectors, 0, uRowCount - 1, 0, uColCount - 1, ref workArray);
                return success;
            }

            bidiagonal.rmatrixbd(ref inputMatrix, numRows, numCols, ref tauQ, ref tauP);
            bidiagonal.rmatrixbdunpackq(ref inputMatrix, numRows, numCols, ref tauQ, uColCount,
                ref leftSingularVectors);
            bidiagonal.rmatrixbdunpackpt(ref inputMatrix, numRows, numCols, ref tauP, vtRowCount,
                ref rightSingularVectors);
            bidiagonal.rmatrixbdunpackdiagonals(ref inputMatrix, numRows, numCols, ref isUpperBidiagonal,
                ref singularValues, ref superDiagonal);

            bool flag4;
            if (additionalMemory < 2 || uNeeded == 0) {
                flag4 = bdsvd.rmatrixbdsvd(ref singularValues, superDiagonal, minDim, isUpperBidiagonal, false,
                    ref leftSingularVectors, uRowCount, ref inputMatrix, 0, ref rightSingularVectors, vtColCount);
            } else {
                var transposedMatrix = new double[minDim, numRows];
                blas.copyandtranspose(ref leftSingularVectors, 0, numRows - 1, 0, minDim - 1, ref transposedMatrix, 0,
                    minDim - 1, 0, numRows - 1);
                flag4 = bdsvd.rmatrixbdsvd(ref singularValues, superDiagonal, minDim, isUpperBidiagonal, false,
                    ref leftSingularVectors, 0, ref transposedMatrix, numRows, ref rightSingularVectors, vtColCount);
                blas.copyandtranspose(ref transposedMatrix, 0, minDim - 1, 0, numRows - 1, ref leftSingularVectors, 0,
                    numRows - 1, 0, minDim - 1);
            }

            return flag4;
        }

        public static bool SVDDecomposition(double[,] inputMatrix, int numRows, int numCols, int uNeeded, int vtNeeded,
            int additionalMemory, ref double[] singularValues, ref double[,] leftSingularVectors,
            ref double[,] rightSingularVectors) {
            var tauQ = new double[0];
            var tauP = new double[0];
            var tau = new double[0];
            var superDiagonal = new double[0];
            var tempMatrix = new double[0, 0];
            var isUpperBidiagonal = false;

            // Clone the input matrix to avoid modifying the original
            inputMatrix = (double[,])inputMatrix.Clone();

            // Check for zero dimensions
            if (numRows == 0 || numCols == 0)
                return true;

            // Determine the smaller dimension between numRows and numCols
            var minDim = Math.Min(numRows, numCols);
            singularValues = new double[minDim + 1];

            // Initialize matrices U and VT based on the needed flags
            int uRowCount = 0, uColCount = 0;
            if (uNeeded == 1) {
                uRowCount = numRows;
                uColCount = minDim;
                leftSingularVectors = new double[uRowCount + 1, uColCount + 1];
            } else if (uNeeded == 2) {
                uRowCount = numRows;
                uColCount = numRows;
                leftSingularVectors = new double[uRowCount + 1, uColCount + 1];
            }

            int vtRowCount = 0, vtColCount = 0;
            if (vtNeeded == 1) {
                vtRowCount = minDim;
                vtColCount = numCols;
                rightSingularVectors = new double[vtRowCount + 1, vtColCount + 1];
            } else if (vtNeeded == 2) {
                vtRowCount = numCols;
                vtColCount = numCols;
                rightSingularVectors = new double[vtRowCount + 1, vtColCount + 1];
            }

            // Case when numRows is significantly larger than numCols
            if (numRows > 1.6 * numCols) {
                if (uNeeded == 0) {
                    // Perform QR decomposition on the input matrix
                    QR.QRDecomposition(ref inputMatrix, numRows, numCols, ref tau);

                    // Set the lower part of the upper triangular matrix to zero
                    for (var i = 2; i <= numCols; ++i)
                    for (var j = 1; j <= i - 1; ++j)
                        inputMatrix[i, j] = 0.0;

                    // Transform to bidiagonal form
                    bidiagonal.tobidiagonal(ref inputMatrix, numCols, numCols, ref tauQ, ref tauP);
                    bidiagonal.unpackptfrombidiagonal(ref inputMatrix, numCols, numCols, ref tauP, vtRowCount,
                        ref rightSingularVectors);
                    bidiagonal.unpackdiagonalsfrombidiagonal(ref inputMatrix, numCols, numCols, ref isUpperBidiagonal,
                        ref singularValues, ref superDiagonal);

                    return bdsvd.bidiagonalsvddecomposition(ref singularValues, superDiagonal, numCols,
                        isUpperBidiagonal, false, ref leftSingularVectors, 0, ref inputMatrix, 0,
                        ref rightSingularVectors, vtColCount);
                }

                // Perform QR decomposition on the input matrix
                QR.QRDecomposition(ref inputMatrix, numRows, numCols, ref tau);
                QR.UnPackQFromQR(ref inputMatrix, numRows, numCols, ref tau, uColCount, ref leftSingularVectors);

                // Set the lower part of the upper triangular matrix to zero
                for (var i = 2; i <= numCols; ++i)
                for (var j = 1; j <= i - 1; ++j)
                    inputMatrix[i, j] = 0.0;

                // Transform to bidiagonal form
                bidiagonal.tobidiagonal(ref inputMatrix, numCols, numCols, ref tauQ, ref tauP);
                bidiagonal.unpackptfrombidiagonal(ref inputMatrix, numCols, numCols, ref tauP, vtRowCount,
                    ref rightSingularVectors);
                bidiagonal.unpackdiagonalsfrombidiagonal(ref inputMatrix, numCols, numCols, ref isUpperBidiagonal,
                    ref singularValues, ref superDiagonal);

                bool flag1;
                if (additionalMemory < 1) {
                    // Multiply by Q to form U
                    bidiagonal.multiplybyqfrombidiagonal(ref inputMatrix, numCols, numCols, ref tauQ,
                        ref leftSingularVectors, numRows, numCols, true, false);
                    flag1 = bdsvd.bidiagonalsvddecomposition(ref singularValues, superDiagonal, numCols,
                        isUpperBidiagonal, false, ref leftSingularVectors, numRows, ref inputMatrix, 0,
                        ref rightSingularVectors, vtColCount);
                } else {
                    var workArray = new double[Math.Max(numRows, numCols) + 1];
                    bidiagonal.unpackqfrombidiagonal(ref inputMatrix, numCols, numCols, ref tauQ, numCols,
                        ref tempMatrix);
                    blas.copymatrix(ref leftSingularVectors, 1, numRows, 1, numCols, ref inputMatrix, 1, numRows, 1,
                        numCols);
                    blas.inplacetranspose(ref tempMatrix, 1, numCols, 1, numCols, ref workArray);
                    flag1 = bdsvd.bidiagonalsvddecomposition(ref singularValues, superDiagonal, numCols,
                        isUpperBidiagonal, false, ref leftSingularVectors, 0, ref tempMatrix, numCols,
                        ref rightSingularVectors, vtColCount);
                    blas.matrixmatrixmultiply(ref inputMatrix, 1, numRows, 1, numCols, false, ref tempMatrix, 1,
                        numCols, 1, numCols, true, 1.0, ref leftSingularVectors, 1, numRows, 1, numCols, 0.0,
                        ref workArray);
                }

                return flag1;
            }

            // Case when numCols is significantly larger than numRows
            if (numCols > 1.6 * numRows) {
                if (vtNeeded == 0) {
                    // Perform LQ decomposition on the input matrix
                    LQ.LQDecomposition(ref inputMatrix, numRows, numCols, ref tau);

                    // Set the upper part of the lower triangular matrix to zero
                    for (var i = 1; i <= numRows - 1; ++i)
                    for (var j = i + 1; j <= numRows; ++j)
                        inputMatrix[i, j] = 0.0;

                    // Transform to bidiagonal form
                    bidiagonal.tobidiagonal(ref inputMatrix, numRows, numRows, ref tauQ, ref tauP);
                    bidiagonal.unpackqfrombidiagonal(ref inputMatrix, numRows, numRows, ref tauQ, uColCount,
                        ref leftSingularVectors);
                    bidiagonal.unpackdiagonalsfrombidiagonal(ref inputMatrix, numRows, numRows, ref isUpperBidiagonal,
                        ref singularValues, ref superDiagonal);

                    var workArray = new double[numRows + 1];
                    blas.inplacetranspose(ref leftSingularVectors, 1, uRowCount, 1, uColCount, ref workArray);
                    var flag2 = bdsvd.bidiagonalsvddecomposition(ref singularValues, superDiagonal, numRows,
                        isUpperBidiagonal, false, ref inputMatrix, 0, ref leftSingularVectors, uRowCount,
                        ref rightSingularVectors, 0);
                    blas.inplacetranspose(ref leftSingularVectors, 1, uRowCount, 1, uColCount, ref workArray);
                    return flag2;
                }

                // Perform LQ decomposition on the input matrix
                LQ.LQDecomposition(ref inputMatrix, numRows, numCols, ref tau);
                LQ.UnPackQFromLQ(ref inputMatrix, numRows, numCols, ref tau, vtRowCount, ref rightSingularVectors);

                // Set the upper part of the lower triangular matrix to zero
                for (var i = 1; i <= numRows - 1; ++i)
                for (var j = i + 1; j <= numRows; ++j)
                    inputMatrix[i, j] = 0.0;

                // Transform to bidiagonal form
                bidiagonal.tobidiagonal(ref inputMatrix, numRows, numRows, ref tauQ, ref tauP);
                bidiagonal.unpackqfrombidiagonal(ref inputMatrix, numRows, numRows, ref tauQ, uColCount,
                    ref leftSingularVectors);
                bidiagonal.unpackdiagonalsfrombidiagonal(ref inputMatrix, numRows, numRows, ref isUpperBidiagonal,
                    ref singularValues, ref superDiagonal);

                var workArray1 = new double[Math.Max(numRows, numCols) + 1];
                blas.inplacetranspose(ref leftSingularVectors, 1, uRowCount, 1, uColCount, ref workArray1);
                bool flag3;
                if (additionalMemory < 1) {
                    bidiagonal.multiplybypfrombidiagonal(ref inputMatrix, numRows, numRows, ref tauP,
                        ref rightSingularVectors, numRows, numCols, false, true);
                    flag3 = bdsvd.bidiagonalsvddecomposition(ref singularValues, superDiagonal, numRows,
                        isUpperBidiagonal, false, ref inputMatrix, 0, ref leftSingularVectors, uRowCount,
                        ref rightSingularVectors, numCols);
                } else {
                    bidiagonal.unpackptfrombidiagonal(ref inputMatrix, numRows, numRows, ref tauP, numRows,
                        ref tempMatrix);
                    flag3 = bdsvd.bidiagonalsvddecomposition(ref singularValues, superDiagonal, numRows,
                        isUpperBidiagonal, false, ref inputMatrix, 0, ref leftSingularVectors, uRowCount,
                        ref tempMatrix, numRows);
                    blas.copymatrix(ref rightSingularVectors, 1, numRows, 1, numCols, ref inputMatrix, 1, numRows, 1,
                        numCols);
                    blas.matrixmatrixmultiply(ref tempMatrix, 1, numRows, 1, numRows, false, ref inputMatrix, 1,
                        numRows, 1, numCols, false, 1.0, ref rightSingularVectors, 1, numRows, 1, numCols, 0.0,
                        ref workArray1);
                }

                blas.inplacetranspose(ref leftSingularVectors, 1, uRowCount, 1, uColCount, ref workArray1);
                return flag3;
            }

            // General case for other dimensions
            if (numRows <= numCols) {
                // Transform to bidiagonal form
                bidiagonal.tobidiagonal(ref inputMatrix, numRows, numCols, ref tauQ, ref tauP);
                bidiagonal.unpackqfrombidiagonal(ref inputMatrix, numRows, numCols, ref tauQ, uColCount,
                    ref leftSingularVectors);
                bidiagonal.unpackptfrombidiagonal(ref inputMatrix, numRows, numCols, ref tauP, vtRowCount,
                    ref rightSingularVectors);
                bidiagonal.unpackdiagonalsfrombidiagonal(ref inputMatrix, numRows, numCols, ref isUpperBidiagonal,
                    ref singularValues, ref superDiagonal);

                var workArray = new double[numRows + 1];
                blas.inplacetranspose(ref leftSingularVectors, 1, uRowCount, 1, uColCount, ref workArray);
                var flag4 = bdsvd.bidiagonalsvddecomposition(ref singularValues, superDiagonal, minDim,
                    isUpperBidiagonal, false, ref inputMatrix, 0, ref leftSingularVectors, uRowCount,
                    ref rightSingularVectors, vtColCount);
                blas.inplacetranspose(ref leftSingularVectors, 1, uRowCount, 1, uColCount, ref workArray);
                return flag4;
            }

            // Transform to bidiagonal form
            bidiagonal.tobidiagonal(ref inputMatrix, numRows, numCols, ref tauQ, ref tauP);
            bidiagonal.unpackqfrombidiagonal(ref inputMatrix, numRows, numCols, ref tauQ, uColCount,
                ref leftSingularVectors);
            bidiagonal.unpackptfrombidiagonal(ref inputMatrix, numRows, numCols, ref tauP, vtRowCount,
                ref rightSingularVectors);
            bidiagonal.unpackdiagonalsfrombidiagonal(ref inputMatrix, numRows, numCols, ref isUpperBidiagonal,
                ref singularValues, ref superDiagonal);

            bool flag5;
            if (additionalMemory < 2 || uNeeded == 0) {
                flag5 = bdsvd.bidiagonalsvddecomposition(ref singularValues, superDiagonal, minDim, isUpperBidiagonal,
                    false, ref leftSingularVectors, uRowCount, ref inputMatrix, 0, ref rightSingularVectors,
                    vtColCount);
            } else {
                var transposedMatrix = new double[minDim + 1, numRows + 1];
                blas.copyandtranspose(ref leftSingularVectors, 1, numRows, 1, minDim, ref transposedMatrix, 1, minDim,
                    1, numRows);
                flag5 = bdsvd.bidiagonalsvddecomposition(ref singularValues, superDiagonal, minDim, isUpperBidiagonal,
                    false, ref leftSingularVectors, 0, ref transposedMatrix, numRows, ref rightSingularVectors,
                    vtColCount);
                blas.copyandtranspose(ref transposedMatrix, 1, minDim, 1, numRows, ref leftSingularVectors, 1, numRows,
                    1, minDim);
            }
            return flag5;
        }
    }
}