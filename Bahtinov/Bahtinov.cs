using NINA.Core.Utility;
using NINA.Core.Utility.Notification;
using System;
using System.Windows;
using System.Windows.Media.Imaging;

namespace CanardConfit.NINA.BahtiFocus.Bahtinov {
    public static class Bahtinov {

        public class Ellipse : BaseINPC {
            public Point Start { get; }
            public int Width { get; }
            public int Height { get; }

            public Ellipse(Point start, int width, int height) {
                Start = start;
                Width = width;
                Height = height;
            }
            
            public Ellipse Scale(double ratioX, double ratioY) {
                return new Ellipse(new Point(Start.X * ratioX, Start.Y * ratioY), (int)(Width * ratioX), (int)(Height * ratioY));
            }
        }
        
        public class Line : BaseINPC {
            public Point P1 { get; }
            public Point P2 { get; }

            public Line(Point p1, Point p2) {
                P1 = p1;
                P2 = p2;
            }
            
            public Line Scale(double ratioX, double ratioY) {
                return new Line(new Point(P1.X * ratioX, P1.Y * ratioY), new Point(P2.X * ratioX, P2.Y * ratioY));
            }
        }
        
        public class BahtinovCalc : BaseINPC {
            public Line LineRight { get; set; }
            public Line LineMiddle { get; set; }
            public Line LineLeft { get; set; }
            public Ellipse EllipseIntersection { get; set; }
            public Ellipse EllipseError { get; set; }
            public Line LineError { get; set; }
            public float FocusError { get; set; }
            public float MaskAngle { get; set; }
            public float Angles1 { get; set; }
            public float Angles2 { get; set; }
            public float Angles3 { get; set; }
            public float AbsoluteFocusError { get; set; }
            public float CriticalFocusThreshold { get; set; }
            public bool CriticalFocus { get; set; }
            public Ellipse[] EllipseCritFocus { get; set; }
            
            public void Scale(double ratioX, double ratioY) {
                LineRight = LineRight?.Scale(ratioX, ratioY);
                LineMiddle = LineMiddle?.Scale(ratioX, ratioY);
                LineLeft = LineLeft?.Scale(ratioX, ratioY);
                EllipseIntersection = EllipseIntersection?.Scale(ratioX, ratioY);
                EllipseError = EllipseError?.Scale(ratioX, ratioY);
                LineError = LineError?.Scale(ratioX, ratioY);
                if (EllipseCritFocus != null) {
                    for (int i = 0; i < EllipseCritFocus.Length; i++) {
                        EllipseCritFocus[i] = EllipseCritFocus[i].Scale(ratioX, ratioY);
                    }
                }
            }
        }
        
        public static (byte[] pixelArray, int stride, int width, int height) GetCenterPixelsFromBitmapSource(BitmapSource bitmapSource, int width, int height)
        {
            if (bitmapSource.Width < width) {
                width = Convert.ToInt32(bitmapSource.Width);
            }
            if (bitmapSource.Height < height) {
                height = Convert.ToInt32(bitmapSource.Height);
            }
            
            // Calculer les coordonnées du rectangle central
            int x = (bitmapSource.PixelWidth - width) / 2;
            int y = (bitmapSource.PixelHeight - height) / 2;
    
            // Définir le rectangle central
            Int32Rect rect = new Int32Rect(x, y, width, height);

            // Calculer le stride pour la région spécifiée
            int stride = width * (bitmapSource.Format.BitsPerPixel / 8);
    
            // Allouer le tableau de pixels pour la région spécifiée
            byte[] pixelArray = new byte[stride * height];
    
            // Copier les pixels de la région spécifiée
            bitmapSource.CopyPixels(rect, pixelArray, stride, 0);
    
            return (pixelArray, stride, width, height);
        }
        
        public static BahtinovCalc CalculateLines(BitmapSource bitmap, ref float[] bahtinovAngles, double diameter, double focalLength, double pixelSize) {
            BahtinovCalc ret = new();
            
            var width = Convert.ToInt32(bitmap.Width);
            var height = Convert.ToInt32(bitmap.Height);
            
            int stride = width * (bitmap.Format.BitsPerPixel / 8);
    
            byte[] pixelArray = new byte[stride * height];
    
            bitmap.CopyPixels(pixelArray, stride, 0);

            // Calculate effective size for the Bahtinov mask analysis
            var effectiveSize = (width < height ? 0.5f * (float)Math.Sqrt(2.0) * width : 0.5f * (float)Math.Sqrt(2.0) * height) - 8f;
            var startX = (int)(0.5 * (width - (double)effectiveSize));
            var endX = (int)(0.5 * (width + (double)effectiveSize));
            var startY = (int)(0.5 * (height - (double)effectiveSize));
            var endY = (int)(0.5 * (height + (double)effectiveSize));

            // Initialize step values for iteration
            var stepX = 1;
            var stepY = stepX;

            // Initialize arrays to store pixel intensities
            var pixelIntensityArray = new float[width, height];

            // Calculate intensity values for each pixel in the effective region
            for (var x = startX; x < endX; ++x)
            {
                for (var y = startY; y < endY; ++y)
                {
                    var pixelIndex = (x + y * width) * stride / width;
                    pixelIntensityArray[x, y] = (float)((double)pixelArray[pixelIndex] +
                                                        pixelArray[pixelIndex + 1] +
                                                        pixelArray[pixelIndex + 2]);
                    pixelIntensityArray[x, y] /= 3f;
                    pixelIntensityArray[x, y] /= byte.MaxValue;
                    pixelIntensityArray[x, y] = (float)Math.Sqrt(pixelIntensityArray[x, y]);
                }
            }

            // Determine center points of the bitmap
            var centerX = (float)((width + 1.0) / 2.0);
            var centerY = (float)((height + 1.0) / 2.0);

            // Arrays to store Bahtinov mask angles and positions
            var bahtinovAngleArray = new float[3];
            var bahtinovPositionArray = new float[3];

            // Variables to store slopes and intercepts of the mask lines
            var angle0Slope = 0.0f;
            var angle0Intercept = 0.0f;

            var angle1X1 = 0.0f;
            var angle1X2 = 0.0f;
            var angle1Y1 = 0.0f;
            var angle1Y2 = 0.0f;

            var angle2Slope = 0.0f;
            var angle2Intercept = 0.0f;

            // Check if Bahtinov angles are uninitialized
            if (bahtinovAngles[0] == 0.0 && bahtinovAngles[1] == 0.0 && bahtinovAngles[2] == 0.0)
            {
                // Number of steps to calculate angles
                var angleStepCount = 180;
                var angleStep = (float)Math.PI / angleStepCount;
                var angleMaxIntensityArray = new float[angleStepCount];
                var anglePositionArray = new float[angleStepCount];
                var interpolatedIntensityArray = new float[width, height];

                // Iterate through each angle step
                for (var angleIndex = 0; angleIndex < angleStepCount; ++angleIndex)
                {
                    var currentAngle = angleStep * angleIndex;
                    var sinAngle = (float)Math.Sin(currentAngle);
                    var cosAngle = (float)Math.Cos(currentAngle);

                    // Rotate the coordinates and interpolate intensity values
                    for (var x = startX; x < endX; x += stepX)
                    {
                        for (var y = startY; y < endY; y += stepY)
                        {
                            var rotateDeltaX = x - centerX;
                            var rotateDeltaY = y - centerY;

                            var rotatedX = (float)(centerX + rotateDeltaX * (double)cosAngle + rotateDeltaY * (double)sinAngle);
                            var rotatedY = (float)(centerY - rotateDeltaX * (double)sinAngle + rotateDeltaY * (double)cosAngle);

                            var floorX = (int)Math.Floor(rotatedX);
                            var ceilX = (int)Math.Ceiling(rotatedX);
                            var floorY = (int)Math.Floor(rotatedY);
                            var ceilY = (int)Math.Ceiling(rotatedY);

                            var fracX = rotatedX - floorX;
                            var fracY = rotatedY - floorY;

                            interpolatedIntensityArray[x, y] = (float)(pixelIntensityArray[floorX, floorY] * (1.0 - fracX) * (1.0 - fracY) +
                                                                      pixelIntensityArray[ceilX, floorY] * (double)fracX * (1.0 - fracY) +
                                                                      pixelIntensityArray[ceilX, ceilY] * (double)fracX * fracY +
                                                                      pixelIntensityArray[floorX, ceilY] * (1.0 - fracX) * fracY);
                        }
                    }

                    // Sum intensities along the vertical axis
                    var verticalSumArray = new float[height];
                    for (var y = 0; y < height; ++y)
                        verticalSumArray[y] = 0.0f;

                    for (var y = startY; y < endY; y += stepY)
                    {
                        var sumCount = 0;
                        for (var x = startX; x < endX; x += stepX)
                        {
                            verticalSumArray[y] += interpolatedIntensityArray[x, y];
                            ++sumCount;
                        }

                        verticalSumArray[y] /= sumCount;
                    }

                    // Smooth the intensity values
                    var smoothedArray = new float[height];
                    for (var y = 0; y < height; ++y)
                        smoothedArray[y] = verticalSumArray[y];

                    for (var y = startY; y < endY; ++y)
                    {
                        smoothedArray[y] = 0.0f;
                        for (var k = -(stepX - 1) / 2; k <= (stepX - 1) / 2; ++k)
                            smoothedArray[y] += verticalSumArray[y + k] / stepX;
                    }

                    for (var y = 0; y < height; ++y)
                        verticalSumArray[y] = smoothedArray[y];

                    // Find the maximum intensity position for the current angle
                    var maxValue = -1f;
                    var maxIndex = -1f;
                    for (var y = startY; y < endY; ++y)
                    {
                        if (verticalSumArray[y] > (double)maxValue)
                        {
                            maxIndex = y;
                            maxValue = verticalSumArray[y];
                        }
                    }

                    try
                    {
                        anglePositionArray[angleIndex] = maxIndex;
                        angleMaxIntensityArray[angleIndex] = maxValue;
                    }
                    catch {
                        Notification.ShowError("Error");
                    }
                }

                // Identify the top 3 angles with maximum intensity
                var angleCount = 0;
                for (var i = 0; i < 3; ++i)
                {
                    var maxAngleValue = -1f;
                    var maxAngleIndex = -1f;
                    var maxAnglePosition = -1f;

                    for (var angleIndex = 0; angleIndex < angleStepCount; ++angleIndex)
                    {
                        if (angleMaxIntensityArray[angleIndex] > (double)maxAngleValue)
                        {
                            maxAngleValue = angleMaxIntensityArray[angleIndex];
                            maxAnglePosition = anglePositionArray[angleIndex];
                            maxAngleIndex = angleIndex * angleStep;
                            angleCount = angleIndex;
                        }
                    }

                    bahtinovPositionArray[i] = maxAnglePosition;
                    bahtinovAngleArray[i] = maxAngleIndex;
                    bahtinovAngles[i] = maxAngleIndex;

                    // Zero out intensity values near the identified angle
                    var angleStepRange = (int)(0.0872664600610733 / angleStep);
                    for (var angleIndex = angleCount - angleStepRange; angleIndex < angleCount + angleStepRange; ++angleIndex)
                    {
                        var index = (angleIndex + angleStepCount) % angleStepCount;
                        try
                        {
                            angleMaxIntensityArray[index] = 0.0f;
                        }
                        catch
                        {
                            Notification.ShowError("Error");
                        }
                    }
                }
            }
            else
            {
                // Process the Bahtinov angles if they are already initialized
                var angleCount = 3;
                var peakPositions = new float[angleCount];
                var peakIntensities = new float[angleCount];
                var interpolatedArray = new float[width, height];

                for (var i = 0; i < angleCount; ++i)
                {
                    var currentAngle = bahtinovAngles[i];
                    var sinAngle = (float)Math.Sin(currentAngle);
                    var cosAngle = (float)Math.Cos(currentAngle);

                    for (var x = startX; x < endX; x += stepX)
                    {
                        for (var y = startY; y < endY; y += stepY)
                        {
                            var rotateDeltaX = x - centerX;
                            var rotateDeltaY = y - centerY;

                            var rotatedX = (float)(centerX + rotateDeltaX * (double)cosAngle + rotateDeltaY * (double)sinAngle);
                            var rotatedY = (float)(centerY - rotateDeltaX * (double)sinAngle + rotateDeltaY * (double)cosAngle);

                            var floorX = (int)Math.Floor(rotatedX);
                            var ceilX = (int)Math.Ceiling(rotatedX);
                            var floorY = (int)Math.Floor(rotatedY);
                            var ceilY = (int)Math.Ceiling(rotatedY);

                            var fracX = rotatedX - floorX;
                            var fracY = rotatedY - floorY;

                            interpolatedArray[x, y] = (float)(pixelIntensityArray[floorX, floorY] * (1.0 - fracX) * (1.0 - fracY) +
                                                              pixelIntensityArray[ceilX, floorY] * (double)fracX * (1.0 - fracY) +
                                                              pixelIntensityArray[ceilX, ceilY] * (double)fracX * fracY +
                                                              pixelIntensityArray[floorX, ceilY] * (1.0 - fracX) * fracY);
                        }
                    }

                    var yvals = new float[height];
                    for (var y = 0; y < height; ++y)
                        yvals[y] = 0.0f;

                    for (var y = startY; y < endY; y += stepY)
                    {
                        var sumCount = 0;
                        for (var x = startX; x < endX; x += stepX)
                        {
                            yvals[y] += interpolatedArray[x, y];
                            ++sumCount;
                        }

                        yvals[y] /= sumCount;
                    }

                    var maxValue = -1f;
                    var estimatedPos = -1f;
                    for (var y = startY; y < endY; ++y)
                    {
                        if (yvals[y] > (double)maxValue)
                        {
                            estimatedPos = y;
                            maxValue = yvals[y];
                        }
                    }

                    var peakPosition = new LsqCalculator().PeakPosition(yvals, (int)estimatedPos, 2);
                    try
                    {
                        peakPositions[i] = peakPosition;
                        peakIntensities[i] = maxValue;
                    }
                    catch
                    {
                        Notification.ShowError("Error");
                    }
                }

                for (var i = 0; i < 3; ++i)
                {
                    bahtinovPositionArray[i] = peakPositions[i];
                    bahtinovAngleArray[i] = bahtinovAngles[i];
                }
            }

            // Sort the angles and positions in ascending order
            for (var i = 0; i < 3; ++i)
            {
                for (var j = i; j < 3; ++j)
                {
                    if (bahtinovAngleArray[j] < (double)bahtinovAngleArray[i])
                    {
                        var tempAngle = bahtinovAngleArray[i];
                        bahtinovAngleArray[i] = bahtinovAngleArray[j];
                        bahtinovAngleArray[j] = tempAngle;

                        var tempPosition = bahtinovPositionArray[i];
                        bahtinovPositionArray[i] = bahtinovPositionArray[j];
                        bahtinovPositionArray[j] = tempPosition;
                    }
                }
            }

            // Adjust angles if the difference between them is too large
            if (bahtinovAngleArray[1] - (double)bahtinovAngleArray[0] > Math.PI / 2.0)
            {
                bahtinovAngleArray[1] -= (float)Math.PI;
                bahtinovAngleArray[2] -= (float)Math.PI;
                bahtinovPositionArray[1] = height - bahtinovPositionArray[1];
                bahtinovPositionArray[2] = height - bahtinovPositionArray[2];
            }

            if (bahtinovAngleArray[2] - (double)bahtinovAngleArray[1] > 1.5707963705062866)
            {
                bahtinovAngleArray[2] -= (float)Math.PI;
                bahtinovPositionArray[2] = height - bahtinovPositionArray[2];
            }

            // Sort the angles and positions again after adjustment
            for (var i = 0; i < 3; ++i)
            {
                for (var j = i; j < 3; ++j)
                {
                    if (bahtinovAngleArray[j] < (double)bahtinovAngleArray[i])
                    {
                        var tempAngle = bahtinovAngleArray[i];
                        bahtinovAngleArray[i] = bahtinovAngleArray[j];
                        bahtinovAngleArray[j] = tempAngle;

                        var tempPosition = bahtinovPositionArray[i];
                        bahtinovPositionArray[i] = bahtinovPositionArray[j];
                        bahtinovPositionArray[j] = tempPosition;
                    }
                }
            }

            // Draw lines based on the calculated angles and positions
            for (var i = 0; i < 3; ++i)
            {
                var minDimension = Math.Min(centerX, centerY);
                var x1 = centerX + -minDimension * (float)Math.Cos(bahtinovAngleArray[i]) +
                         (bahtinovPositionArray[i] - centerY) * (float)Math.Sin(bahtinovAngleArray[i]);
                var x2 = centerX + minDimension * (float)Math.Cos(bahtinovAngleArray[i]) +
                         (bahtinovPositionArray[i] - centerY) * (float)Math.Sin(bahtinovAngleArray[i]);
                var y1 = centerY + -minDimension * (float)Math.Sin(bahtinovAngleArray[i]) +
                         (float)-(bahtinovPositionArray[i] - (double)centerY) * (float)Math.Cos(bahtinovAngleArray[i]);
                var y2 = centerY + minDimension * (float)Math.Sin(bahtinovAngleArray[i]) +
                         (float)-(bahtinovPositionArray[i] - (double)centerY) * (float)Math.Cos(bahtinovAngleArray[i]);

                if (i == 0)
                {
                    var slope = (float)((y2 - (double)y1) / (x2 - (double)x1));
                    var intercept = (float)(-(double)x1 * ((y2 - (double)y1) / (x2 - (double)x1))) + y1;
                    angle0Slope = slope;
                    angle0Intercept = intercept;
                }
                else if (i == 1)
                {
                    angle1X1 = x1;
                    angle1X2 = x2;
                    angle1Y1 = y1;
                    angle1Y2 = y2;
                }
                else if (i == 2)
                {
                    var slope = (float)((y2 - (double)y1) / (x2 - (double)x1));
                    var intercept = (float)(-(double)x1 * ((y2 - (double)y1) / (x2 - (double)x1))) + y1;
                    angle2Slope = slope;
                    angle2Intercept = intercept;
                }

                Point p1 = new Point(x1, height - y1);
                Point p2 = new Point(x2, height - y2);
                
                switch (i) {
                    case 0:
                        ret.LineRight = new Line(p1, p2);
                        break;
                    case 1:
                        ret.LineMiddle = new Line(p1, p2);
                        break;
                    case 2:
                        ret.LineLeft = new Line(p1, p2);
                        break;
                }
            }

            // Calculate the intersection of the lines
            var intersectionX = (float)(-(angle0Intercept - (double)angle2Intercept) / (angle0Slope - (double)angle2Slope));
            var intersectionY = angle0Slope * intersectionX + angle0Intercept;

            var ellipseRadius = 8;
            
            ret.EllipseIntersection = new Ellipse(new Point(intersectionX - ellipseRadius, height - intersectionY - ellipseRadius), ellipseRadius * 2, ellipseRadius * 2);
            
            // Calculate projection factor
            var projectionFactor = (float)((intersectionX - (double)angle1X1) * (angle1X2 - (double)angle1X1) +
                                            (intersectionY - (double)angle1Y1) * (angle1Y2 - (double)angle1Y1)) /
                                   (float)((angle1X2 - (double)angle1X1) * (angle1X2 - (double)angle1X1) +
                                           (angle1Y2 - (double)angle1Y1) * (angle1Y2 - (double)angle1Y1));

            var projectedX = angle1X1 + projectionFactor * (angle1X2 - angle1X1);
            var projectedY = angle1Y1 + projectionFactor * (angle1Y2 - angle1Y1);

            // Calculate error distance and sign
            var errorDistance = (float)Math.Sqrt((intersectionX - (double)projectedX) * (intersectionX - (double)projectedX) +
                                                 (intersectionY - (double)projectedY) * (intersectionY - (double)projectedY));

            var deltaX = intersectionX - projectedX;
            var deltaY = intersectionY - projectedY;
            var directionX = angle1X2 - angle1X1;
            var directionY = angle1Y2 - angle1Y1;
            
            float errorSign = -Math.Sign((float)(deltaX * (double)directionY - deltaY * (double)directionX));

            // Draw projected point and error line
            var projectedXLong = intersectionX + (float)((projectedX - (double)intersectionX) * 20.0);
            var projectedYLong = intersectionY + (float)((projectedY - (double)intersectionY) * 20.0);

            var ellipseLongRadius = 8;
            
            ret.EllipseError = new Ellipse(new Point(projectedXLong - ellipseLongRadius, height - projectedYLong - ellipseLongRadius), ellipseLongRadius * 2, ellipseLongRadius * 2);

            ret.LineError = new Line(new Point(projectedXLong, height - projectedYLong), new Point(intersectionX, height - intersectionY));
            
            // Display focus error information

            ret.FocusError = errorSign * errorDistance;

            // Calculate and display Bahtinov mask information
            var degreeFactor = 57.2957764f;
            var radianFactor = (float)Math.PI / 180f;
            var averageAngle = Math.Abs((float)((bahtinovAngleArray[2] - (double)bahtinovAngleArray[0]) / 2.0));
            var calculatedError = (float)(9.0 / 32.0 * ((float)diameter / ((float)focalLength * pixelSize)) * 
                                          (1.0 + Math.Cos(45.0 * radianFactor) * (1.0 + Math.Tan(averageAngle))));

            ret.MaskAngle = averageAngle * degreeFactor;
            
            ret.Angles1 = degreeFactor * bahtinovAngleArray[0];
            ret.Angles2 = degreeFactor * bahtinovAngleArray[1];
            ret.Angles3 = degreeFactor * bahtinovAngleArray[2];

            ret.AbsoluteFocusError = errorSign * errorDistance / calculatedError;

            // Calculate critical focus threshold and display information
            ret.CriticalFocusThreshold = (float)(8.9999997499035089E-07 *
                                                 (focalLength / diameter) *
                                                 (focalLength / diameter));
            
            ret.CriticalFocus = Math.Abs(ret.AbsoluteFocusError * 1E-06) < Math.Abs(ret.CriticalFocusThreshold);

            // Highlight the center point if within critical focus
            if (ret.CriticalFocus) {
                
                Ellipse[] list = new Ellipse[3];
                
                for (var i = 32; i < 128; i += 32)
                {
                    list[(i / 32) - 1] = new Ellipse(new Point(intersectionX - i, height - intersectionY - i), i * 2, i * 2);
                }

                ret.EllipseCritFocus = list;
            }

            return ret;
        }
    }
}