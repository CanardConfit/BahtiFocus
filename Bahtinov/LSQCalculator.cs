using System;

namespace CanardConfit.NINA.BahtiFocus.Bahtinov {
    public class LsqCalculator {
        /// <summary>
        /// Calculate the base functions for least squares fitting.
        /// </summary>
        /// <param name="baseFunctions">2D array to store the base functions</param>
        /// <param name="xValues">Array of x-values</param>
        /// <param name="numFunctions">Number of functions to use in fitting</param>
        /// <param name="numPoints">Number of data points</param>
        private void CalculateBaseFunc(ref double[,] baseFunctions, ref float[] xValues, int numFunctions,
            int numPoints) {
            for (var pointIndex = 0; pointIndex < numPoints; ++pointIndex)
            for (var funcIndex = 0; funcIndex < numFunctions; ++funcIndex)
                baseFunctions[pointIndex, funcIndex] = Math.Pow(xValues[pointIndex], numFunctions - 1 - funcIndex);
        }

        /// <summary>
        /// Calculate the peak position using least squares fitting.
        /// </summary>
        /// <param name="yValues">Array of y-values</param>
        /// <param name="estimatedPosition">Estimated position of the peak</param>
        /// <param name="halfFitRange">Half range of the fit</param>
        /// <returns>Calculated peak position</returns>
        public float PeakPosition(float[] yValues, int estimatedPosition, int halfFitRange) {
            try {
                // Define the range of x-values around the estimated peak position
                var xValues = new float[2 * halfFitRange + 1];
                var yFitValues = new float[2 * halfFitRange + 1];
                for (var offset = -halfFitRange; offset <= halfFitRange; ++offset) {
                    xValues[offset + halfFitRange] = offset + estimatedPosition;
                    yFitValues[offset + halfFitRange] = yValues[offset + estimatedPosition];
                }

                // Initialize arrays for least squares calculation
                var numPoints = xValues.Length;
                var weights = new double[numPoints];
                var coefficients = new double[3];
                var baseFunctions = new double[numPoints, 3];

                // Set all weights to 1.0 (equal weighting)
                for (var i = 0; i < numPoints; ++i)
                    weights[i] = 1.0;

                // Calculate base functions
                CalculateBaseFunc(ref baseFunctions, ref xValues, 3, numPoints);

                // Convert y-values to double for least squares calculation
                var yDouble = new double[yFitValues.Length];
                for (var i = 0; i < yFitValues.Length; ++i)
                    yDouble[i] = yFitValues[i];

                // Perform the least squares fitting
                LeastSquares.BuildGeneralLeastSquares(ref yDouble, ref weights, ref baseFunctions, numPoints, 3,
                    ref coefficients);

                // Calculate the peak position from the coefficients
                var peakPosition = (float)(-coefficients[1] / (2.0 * coefficients[0]));
                return peakPosition;
            } catch {
                // Return 0.0f in case of an error
                return 0.0f;
            }
        }
    }
}