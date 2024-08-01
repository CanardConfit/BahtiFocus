using System;

namespace CanardConfit.NINA.BahtiFocus.Bahtinov {
    public static class LQ {
        public static void RMatrixLQ(ref double[,] matrix, int rows, int cols, ref double[] tau) {
            var reflectionTau = 0.0;
            var minDim = Math.Min(rows, cols);
            var work = new double[rows + 1];
            var columnVector = new double[cols + 1];
            tau = new double[minDim];

            // Loop over each column to compute the LQ decomposition
            for (var k = 0; k < minDim; ++k) {
                // Copy the k-th column of the matrix starting from the diagonal element
                for (var i = 1; i <= cols - k; ++i)
                    columnVector[i] = matrix[k, i + k - 1];

                // Generate the Householder reflection vector
                Reflections.GenerateReflection(ref columnVector, cols - k, ref reflectionTau);
                tau[k] = reflectionTau;

                // Apply the reflection to the k-th column
                for (var i = k; i < cols; ++i)
                    matrix[k, i] = columnVector[i - k + 1];
                columnVector[1] = 1.0;

                // Apply the reflection to the remaining submatrix
                if (k < cols)
                    Reflections.ApplyReflectionFromTheRight(ref matrix, tau[k], ref columnVector, k + 1, rows - 1, k,
                        cols - 1, ref work);
            }
        }

        public static void RMatrixLQUnPackQ(
            ref double[,] lqMatrix,
            int rows,
            int cols,
            ref double[] tau,
            int qRows,
            ref double[,] q) {
            // Check the basic conditions
            if (rows <= 0 || cols <= 0 || qRows <= 0)
                return;

            var minDim = Math.Min(Math.Min(rows, cols), qRows);
            q = new double[qRows, cols];
            var reflectionVector = new double[cols + 1];
            var work = new double[qRows + 1];

            // Initialize the Q matrix as an identity matrix
            for (var i = 0; i < qRows; ++i)
            for (var j = 0; j < cols; ++j)
                q[i, j] = i == j ? 1.0 : 0.0;

            // Applying Householder thinking to rebuild Q
            for (var k = minDim - 1; k >= 0; --k) {
                // Copy the k-th reflection vector from lqMatrix
                for (var i = 1; i <= cols - k; ++i)
                    reflectionVector[i] = lqMatrix[k, i + k - 1];
                reflectionVector[1] = 1.0;

                // Apply the reflection to Q
                Reflections.ApplyReflectionFromTheRight(ref q, tau[k], ref reflectionVector, 0, qRows - 1, k, cols - 1,
                    ref work);
            }
        }

        public static void RMatrixLQUnPackL(ref double[,] sourceMatrix, int rows, int columns,
            ref double[,] lowerMatrix) {
            // Check for invalid matrix dimensions
            if ((rows <= 0) || (columns <= 0))
                return;

            // Initialize the lower triangular matrix with the same dimensions as the source matrix
            lowerMatrix = new double[rows, columns];

            // Set the first row of the lower triangular matrix to zeros
            for (var col = 0; col <= columns - 1; ++col)
                lowerMatrix[0, col] = 0.0;

            // Set the remaining rows of the lower triangular matrix to the same values as the first row
            for (var row = 1; row <= rows - 1; ++row)
            for (var col = 0; col <= columns - 1; ++col)
                lowerMatrix[row, col] = lowerMatrix[0, col];

            // Copy the elements from the source matrix to the lower triangular matrix
            for (var row = 0; row <= rows - 1; ++row) {
                // Determine the boundary for the columns to be copied based on the current row
                var minIndex = Math.Min(row, columns - 1);

                // Copy the elements from the source matrix to the lower triangular matrix
                for (var col = 0; col <= minIndex; ++col)
                    lowerMatrix[row, col] = sourceMatrix[row, col];
            }
        }

        public static void LQDecomposition(ref double[,] matrix, int rows, int cols, ref double[] tau) {
            var reflectionTau = 0.0;
            var minDim = Math.Min(rows, cols);
            var work = new double[rows + 1];
            var reflectionVector = new double[cols + 1];
            tau = new double[minDim + 1];

            // Loop over each column to compute the LQ decomposition
            for (var k = 1; k <= minDim; ++k) {
                var remainingCols = cols - k + 1;
                var colOffset = k - 1;

                // Copy the k-th column of the matrix starting from the diagonal element
                for (var i = 1; i <= remainingCols; ++i)
                    reflectionVector[i] = matrix[k, i + colOffset];

                // Generate the Householder reflection vector
                Reflections.GenerateReflection(ref reflectionVector, remainingCols, ref reflectionTau);
                tau[k] = reflectionTau;

                // Apply the reflection to the k-th column
                var vectorOffset = 1 - k;
                for (var i = k; i <= cols; ++i)
                    matrix[k, i] = reflectionVector[i + vectorOffset];

                reflectionVector[1] = 1.0;

                // Apply the reflection to the remaining submatrix
                if (k < cols)
                    Reflections.ApplyReflectionFromTheRight(ref matrix, tau[k], ref reflectionVector, k + 1, rows, k,
                        cols, ref work);
            }
        }

        public static void UnPackQFromLQ(
            ref double[,] lqMatrix,
            int numRows,
            int numCols,
            ref double[] tau,
            int numQRows,
            ref double[,] qMatrix) {
            // Check for invalid matrix dimensions or zero Q rows
            if (numRows == 0 || numCols == 0 || numQRows == 0)
                return;

            // Determine the minimum dimension for the transformation
            var minDim = Math.Min(Math.Min(numRows, numCols), numQRows);

            // Initialize the Q matrix with dimensions (numQRows + 1, numCols + 1)
            qMatrix = new double[numQRows + 1, numCols + 1];

            // Temporary vectors for reflection application
            var tempVecV = new double[numCols + 1];
            var tempWork = new double[numQRows + 1];

            // Initialize the Q matrix to an identity matrix
            for (var row = 1; row <= numQRows; ++row)
            for (var col = 1; col <= numCols; ++col)
                qMatrix[row, col] = row != col ? 0.0 : 1.0;

            // Apply reflections to form the Q matrix
            for (var reflectionIndex = minDim; reflectionIndex >= 1; --reflectionIndex) {
                // Determine the current column range for reflection
                var colStart = numCols - reflectionIndex + 1;
                var colOffset = reflectionIndex - 1;

                // Copy the current column of the LQ matrix into tempVecV
                for (var index = 1; index <= colStart; ++index)
                    tempVecV[index] = lqMatrix[reflectionIndex, index + colOffset];

                // Set the first element of tempVecV to 1.0
                tempVecV[1] = 1.0;

                // Apply the reflection from the right to form Q
                Reflections.ApplyReflectionFromTheRight(
                    ref qMatrix,
                    tau[reflectionIndex],
                    ref tempVecV,
                    1,
                    numQRows,
                    reflectionIndex,
                    numCols,
                    ref tempWork);
            }
        }

        public static void LQDecompositionUnPacked(
            double[,] inputMatrix,
            int numRows,
            int numCols,
            ref double[,] lowerMatrix,
            ref double[,] orthogonalMatrix) {
            // Temporary array to store tau values for the decomposition
            var tau = new double[0];

            // Clone the input matrix to avoid modifying the original
            inputMatrix = (double[,])inputMatrix.Clone();

            // Check if the number of columns is non-positive
            if (numCols <= 0)
                return;

            // Initialize the orthogonal matrix Q with dimensions (numCols + 1, numCols + 1)
            orthogonalMatrix = new double[numCols + 1, numCols + 1];

            // Initialize the lower triangular matrix L with dimensions (numRows + 1, numCols + 1)
            lowerMatrix = new double[numRows + 1, numCols + 1];

            // Perform the LQ decomposition on the input matrix
            LQDecomposition(ref inputMatrix, numRows, numCols, ref tau);

            // Extract the lower triangular matrix L from the decomposed matrix
            for (var row = 1; row <= numRows; ++row)
            for (var col = 1; col <= numCols; ++col)
                lowerMatrix[row, col] = col <= row ? inputMatrix[row, col] : 0.0;

            // Extract the orthogonal matrix Q from the decomposed matrix
            UnPackQFromLQ(ref inputMatrix, numRows, numCols, ref tau, numCols, ref orthogonalMatrix);
        }
    }
}