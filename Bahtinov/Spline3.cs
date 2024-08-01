
using System;

namespace CanardConfit.NINA.BahtiFocus.Bahtinov {
    public static class Spline3
    {
        public static void BuildLinearSpline(double[] x, double[] y, int n, ref double[] c) {
            x = (double[])x.Clone();
            y = (double[])y.Clone();
            HeapSortPoints(ref x, ref y, n);
            var num = 3 + n + (n - 1) * 4;
            c = new double[num - 1 + 1];
            c[0] = num;
            c[1] = 3.0;
            c[2] = n;
            for (var index = 0; index <= n - 1; ++index)
                c[3 + index] = x[index];
            for (var index = 0; index <= n - 2; ++index) {
                c[3 + n + 4 * index] = y[index];
                c[3 + n + 4 * index + 1] = (y[index + 1] - y[index]) / (x[index + 1] - x[index]);
                c[3 + n + 4 * index + 2] = 0.0;
                c[3 + n + 4 * index + 3] = 0.0;
            }
        }

        public static void BuildCubicSpline(double[] x, double[] y, int n, int boundltype, double boundl, int boundrtype, double boundr, ref double[] c) {
            var x1 = new double[0];
            x = (double[])x.Clone();
            y = (double[])y.Clone();
            var a = new double[n - 1 + 1];
            var b = new double[n - 1 + 1];
            var c1 = new double[n - 1 + 1];
            var d = new double[n - 1 + 1];
            if (n == 2 && boundltype == 0 && boundrtype == 0) {
                boundltype = 2;
                boundl = 0.0;
                boundrtype = 2;
                boundr = 0.0;
            }

            HeapSortPoints(ref x, ref y, n);
            if (boundltype == 0) {
                a[0] = 0.0;
                b[0] = 1.0;
                c1[0] = 1.0;
                d[0] = 2.0 * (y[1] - y[0]) / (x[1] - x[0]);
            }

            if (boundltype == 1) {
                a[0] = 0.0;
                b[0] = 1.0;
                c1[0] = 0.0;
                d[0] = boundl;
            }

            if (boundltype == 2) {
                a[0] = 0.0;
                b[0] = 2.0;
                c1[0] = 1.0;
                d[0] = 3.0 * (y[1] - y[0]) / (x[1] - x[0]) - 0.5 * boundl * (x[1] - x[0]);
            }

            for (var index = 1; index <= n - 2; ++index) {
                a[index] = x[index + 1] - x[index];
                b[index] = 2.0 * (x[index + 1] - x[index - 1]);
                c1[index] = x[index] - x[index - 1];
                d[index] = 3.0 * (y[index] - y[index - 1]) / (x[index] - x[index - 1]) * (x[index + 1] - x[index]) +
                           3.0 * (y[index + 1] - y[index]) / (x[index + 1] - x[index]) * (x[index] - x[index - 1]);
            }

            if (boundrtype == 0) {
                a[n - 1] = 1.0;
                b[n - 1] = 1.0;
                c1[n - 1] = 0.0;
                d[n - 1] = 2.0 * (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
            }

            if (boundrtype == 1) {
                a[n - 1] = 0.0;
                b[n - 1] = 1.0;
                c1[n - 1] = 0.0;
                d[n - 1] = boundr;
            }

            if (boundrtype == 2) {
                a[n - 1] = 1.0;
                b[n - 1] = 2.0;
                c1[n - 1] = 0.0;
                d[n - 1] = 3.0 * (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]) + 0.5 * boundr * (x[n - 1] - x[n - 2]);
            }

            SolveTriDiagonal(a, b, c1, d, n, ref x1);
            BuildHermiteSpline(x, y, x1, n, ref c);
        }

        public static void BuildHermiteSpline(double[] x, double[] y, double[] d, int n, ref double[] c) {
            x = (double[])x.Clone();
            y = (double[])y.Clone();
            d = (double[])d.Clone();
            HeapSortDPoints(ref x, ref y, ref d, n);
            var num1 = 3 + n + (n - 1) * 4;
            c = new double[num1 - 1 + 1];
            c[0] = num1;
            c[1] = 3.0;
            c[2] = n;
            for (var index = 0; index <= n - 1; ++index)
                c[3 + index] = x[index];
            for (var index = 0; index <= n - 2; ++index) {
                var X = x[index + 1] - x[index];
                var num2 = CustomMath.Sqr(X);
                var num3 = X * num2;
                c[3 + n + 4 * index] = y[index];
                c[3 + n + 4 * index + 1] = d[index];
                c[3 + n + 4 * index + 2] = (3.0 * (y[index + 1] - y[index]) - 2.0 * d[index] * X - d[index + 1] * X) / num2;
                c[3 + n + 4 * index + 3] = (2.0 * (y[index] - y[index + 1]) + d[index] * X + d[index + 1] * X) / num3;
            }
        }

        public static void BuildAkimaSpline(double[] x, double[] y, int n, ref double[] c) {
            x = (double[])x.Clone();
            y = (double[])y.Clone();
            HeapSortPoints(ref x, ref y, n);
            var numArray4 = new double[n - 2 + 1];
            var numArray5 = new double[n - 2 + 1];
            for (var index = 0; index <= n - 2; ++index)
                numArray5[index] = (y[index + 1] - y[index]) / (x[index + 1] - x[index]);
            for (var index = 1; index <= n - 2; ++index)
                numArray4[index] = Math.Abs(numArray5[index] - numArray5[index - 1]);
            var d = new double[n - 1 + 1];
            for (var index = 2; index <= n - 3; ++index)
                d[index] = Math.Abs(numArray4[index - 1]) + Math.Abs(numArray4[index + 1]) == 0.0
                    ? ((x[index + 1] - x[index]) * numArray5[index - 1] + (x[index] - x[index - 1]) * numArray5[index]) /
                      (x[index + 1] - x[index - 1])
                    : (numArray4[index + 1] * numArray5[index - 1] + numArray4[index - 1] * numArray5[index]) /
                      (numArray4[index + 1] + numArray4[index - 1]);
            d[0] = DiffThreePoint(x[0], x[0], y[0], x[1], y[1], x[2], y[2]);
            d[1] = DiffThreePoint(x[1], x[0], y[0], x[1], y[1], x[2], y[2]);
            d[n - 2] = DiffThreePoint(x[n - 2], x[n - 3], y[n - 3], x[n - 2], y[n - 2], x[n - 1], y[n - 1]);
            d[n - 1] = DiffThreePoint(x[n - 1], x[n - 3], y[n - 3], x[n - 2], y[n - 2], x[n - 1], y[n - 1]);
            BuildHermiteSpline(x, y, d, n, ref c);
        }

        public static double SplineInterpolation(ref double[] c, double x) {
            var num1 = (int)Math.Round(c[2]);
            var index1 = 3;
            var num2 = 3 + num1 - 2 + 1;
            while (index1 != num2 - 1) {
                var index2 = (index1 + num2) / 2;
                if (c[index2] >= x)
                    num2 = index2;
                else
                    index1 = index2;
            }

            x -= c[index1];
            var index3 = 3 + num1 + 4 * (index1 - 3);
            return c[index3] + x * (c[index3 + 1] + x * (c[index3 + 2] + x * c[index3 + 3]));
        }

        public static void SplineDifferentiation(ref double[] c, double x, ref double s, ref double ds, ref double d2s) {
            var num1 = (int)Math.Round(c[2]);
            var index1 = 3;
            var num2 = 3 + num1 - 2 + 1;
            while (index1 != num2 - 1) {
                var index2 = (index1 + num2) / 2;
                if (c[index2] >= x)
                    num2 = index2;
                else
                    index1 = index2;
            }

            x -= c[index1];
            var index3 = 3 + num1 + 4 * (index1 - 3);
            s = c[index3] + x * (c[index3 + 1] + x * (c[index3 + 2] + x * c[index3 + 3]));
            ds = c[index3 + 1] + 2.0 * x * c[index3 + 2] + 3.0 * CustomMath.Sqr(x) * c[index3 + 3];
            d2s = 2.0 * c[index3 + 2] + 6.0 * x * c[index3 + 3];
        }

        public static void SplineCopy(ref double[] c, ref double[] cc) {
            var num = (int)Math.Round(c[0]);
            cc = new double[num - 1 + 1];
            for (var index = 0; index <= num - 1; ++index)
                cc[index] = c[index];
        }

        public static void SplineUnPack(ref double[] c, ref int n, ref double[,] tbl) {
            n = (int)Math.Round(c[2]);
            tbl = new double[n - 2 + 1, 6];
            for (var index = 0; index <= n - 2; ++index) {
                tbl[index, 0] = c[3 + index];
                tbl[index, 1] = c[3 + index + 1];
                tbl[index, 2] = c[3 + n + 4 * index];
                tbl[index, 3] = c[3 + n + 4 * index + 1];
                tbl[index, 4] = c[3 + n + 4 * index + 2];
                tbl[index, 5] = c[3 + n + 4 * index + 3];
            }
        }

        public static void SplineLinTransX(ref double[] c, double a, double b) {
            var s = 0.0;
            var ds = 0.0;
            var d2s = 0.0;
            var n = (int)Math.Round(c[2]);
            if (a == 0.0) {
                var num = SplineInterpolation(ref c, b);
                for (var index = 0; index <= n - 2; ++index) {
                    c[3 + n + 4 * index] = num;
                    c[3 + n + 4 * index + 1] = 0.0;
                    c[3 + n + 4 * index + 2] = 0.0;
                    c[3 + n + 4 * index + 3] = 0.0;
                }
            } else {
                var x = new double[n - 1 + 1];
                var y = new double[n - 1 + 1];
                var d = new double[n - 1 + 1];
                for (var index = 0; index <= n - 1; ++index) {
                    x[index] = c[3 + index];
                    SplineDifferentiation(ref c, x[index], ref s, ref ds, ref d2s);
                    x[index] = (x[index] - b) / a;
                    y[index] = s;
                    d[index] = a * ds;
                }

                BuildHermiteSpline(x, y, d, n, ref c);
            }
        }

        public static void SplineLinTransY(ref double[] c, double a, double b) {
            var num = (int)Math.Round(c[2]);
            for (var index = 0; index <= num - 2; ++index) {
                c[3 + num + 4 * index] = a * c[3 + num + 4 * index] + b;
                c[3 + num + 4 * index + 1] = a * c[3 + num + 4 * index + 1];
                c[3 + num + 4 * index + 2] = a * c[3 + num + 4 * index + 2];
                c[3 + num + 4 * index + 3] = a * c[3 + num + 4 * index + 3];
            }
        }

        public static double SplineIntegration(ref double[] c, double x) {
            var num1 = (int)Math.Round(c[2]);
            var index1 = 3;
            var num2 = 3 + num1 - 2 + 1;
            while (index1 != num2 - 1) {
                var index2 = (index1 + num2) / 2;
                if (c[index2] >= x)
                    num2 = index2;
                else
                    index1 = index2;
            }

            var num3 = 0.0;
            for (var index3 = 3; index3 <= index1 - 1; ++index3) {
                var X = c[index3 + 1] - c[index3];
                var index4 = 3 + num1 + 4 * (index3 - 3);
                num3 = num3 + c[index4] * X + c[index4 + 1] * CustomMath.Sqr(X) / 2.0 + c[index4 + 2] * CustomMath.Sqr(X) * X / 3.0 +
                       c[index4 + 3] * CustomMath.Sqr(CustomMath.Sqr(X)) / 4.0;
            }

            var X1 = x - c[index1];
            var index5 = 3 + num1 + 4 * (index1 - 3);
            return num3 + c[index5] * X1 + c[index5 + 1] * CustomMath.Sqr(X1) / 2.0 + c[index5 + 2] * CustomMath.Sqr(X1) * X1 / 3.0 +
                   c[index5 + 3] * CustomMath.Sqr(CustomMath.Sqr(X1)) / 4.0;
        }

        public static void Spline3BuildTable(int n, int diffn, double[] x, double[] y, double boundl, double boundr, ref double[,] ctbl) {
            x = (double[])x.Clone();
            y = (double[])y.Clone();
            --n;
            var num1 = (n + 1) / 2;
            do {
                var num2 = num1;
                do {
                    var index = num2 - num1;
                    var flag = true;
                    do {
                        if (x[index] <= x[index + num1]) {
                            flag = false;
                        } else {
                            var num3 = x[index];
                            x[index] = x[index + num1];
                            x[index + num1] = num3;
                            var num4 = y[index];
                            y[index] = y[index + num1];
                            y[index + num1] = num4;
                        }

                        --index;
                    } while (index >= 0 && flag);

                    ++num2;
                } while (num2 <= n);

                num1 /= 2;
            } while (num1 > 0);

            ctbl = new double[5, n + 1];
            ++n;
            double num5;
            double num6;
            double num7;
            double num8;
            if (diffn == 1) {
                num5 = 1.0;
                num6 = 6.0 / (x[1] - x[0]) * ((y[1] - y[0]) / (x[1] - x[0]) - boundl);
                num7 = 1.0;
                num8 = 6.0 / (x[n - 1] - x[n - 2]) * (boundr - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
            } else {
                num5 = 0.0;
                num6 = 2.0 * boundl;
                num7 = 0.0;
                num8 = 2.0 * boundr;
            }

            var num9 = n - 1;
            if (n < 2)
                return;
            if (n > 2) {
                var num10 = x[1] - x[0];
                var num11 = y[1] - y[0];
                for (var index = 2; index <= num9; ++index) {
                    var num12 = x[index] - x[index - 1];
                    var num13 = y[index] - y[index - 1];
                    var num14 = num10 + num12;
                    ctbl[1, index - 1] = num12 / num14;
                    ctbl[2, index - 1] = 1.0 - ctbl[1, index - 1];
                    ctbl[3, index - 1] = 6.0 * (num13 / num12 - num11 / num10) / num14;
                    num10 = num12;
                    num11 = num13;
                }
            }

            ctbl[1, 0] = -(num5 / 2.0);
            ctbl[2, 0] = num6 / 2.0;
            if (n != 2) {
                for (var index = 2; index <= num9; ++index) {
                    var num15 = ctbl[2, index - 1] * ctbl[1, index - 2] + 2.0;
                    ctbl[1, index - 1] = -(ctbl[1, index - 1] / num15);
                    ctbl[2, index - 1] = (ctbl[3, index - 1] - ctbl[2, index - 1] * ctbl[2, index - 2]) / num15;
                }
            }

            var num16 = (num8 - num7 * ctbl[2, num9 - 1]) / (num7 * ctbl[1, num9 - 1] + 2.0);
            for (var index1 = 1; index1 <= num9; ++index1) {
                var index2 = n - index1;
                var num17 = ctbl[1, index2 - 1] * num16 + ctbl[2, index2 - 1];
                var num18 = x[index2] - x[index2 - 1];
                ctbl[3, index2 - 1] = (num16 - num17) / num18 / 6.0;
                ctbl[2, index2 - 1] = num17 / 2.0;
                ctbl[1, index2 - 1] = (y[index2] - y[index2 - 1]) / num18 -
                                      (ctbl[2, index2 - 1] + ctbl[3, index2 - 1] * num18) * num18;
                num16 = num17;
            }

            for (var index = 1; index <= n; ++index) {
                ctbl[0, index - 1] = y[index - 1];
                ctbl[4, index - 1] = x[index - 1];
            }
        }

        public static double Spline3Interpolate(int n, ref double[,] c, double x) {
            --n;
            var num1 = n;
            var num2 = 0;
            while (num1 > 0) {
                var num3 = num1 / 2;
                var index = num2 + num3;
                if (c[4, index] < x) {
                    num2 = index + 1;
                    num1 = num1 - num3 - 1;
                }
                else
                {
                    num1 = num3;
                }
            }

            var index1 = num2 - 1;
            if (index1 < 0)
                index1 = 0;
            return c[0, index1] + (x - c[4, index1]) *
                (c[1, index1] + (x - c[4, index1]) * (c[2, index1] + c[3, index1] * (x - c[4, index1])));
        }

        private static void HeapSortPoints(ref double[] x, ref double[] y, int n)
        {
            var flag1 = true;
            var flag2 = true;
            for (var index = 1; index <= n - 1; ++index)
            {
                flag1 &= x[index] > x[index - 1];
                flag2 &= x[index] < x[index - 1];
            }

            if (flag1)
                return;
            if (flag2)
            {
                for (var index1 = 0; index1 <= n - 1; ++index1)
                {
                    var index2 = n - 1 - index1;
                    if (index2 <= index1)
                        break;
                    var num1 = x[index1];
                    x[index1] = x[index2];
                    x[index2] = num1;
                    var num2 = y[index1];
                    y[index1] = y[index2];
                    y[index2] = num2;
                }
            }
            else
            {
                if (n == 1)
                    return;
                var num3 = 2;
                do
                {
                    var num4 = num3;
                    while (num4 != 1)
                    {
                        var num5 = num4 / 2;
                        if (x[num5 - 1] >= x[num4 - 1])
                        {
                            num4 = 1;
                        }
                        else
                        {
                            var num6 = x[num5 - 1];
                            x[num5 - 1] = x[num4 - 1];
                            x[num4 - 1] = num6;
                            var num7 = y[num5 - 1];
                            y[num5 - 1] = y[num4 - 1];
                            y[num4 - 1] = num7;
                            num4 = num5;
                        }
                    }

                    ++num3;
                } while (num3 <= n);

                var index3 = n - 1;
                do
                {
                    var num8 = x[index3];
                    x[index3] = x[0];
                    x[0] = num8;
                    var num9 = y[index3];
                    y[index3] = y[0];
                    y[0] = num9;
                    var num10 = 1;
                    while (num10 != 0)
                    {
                        var index4 = 2 * num10;
                        if (index4 > index3)
                        {
                            num10 = 0;
                        }
                        else
                        {
                            if (index4 < index3 && x[index4] > x[index4 - 1])
                                ++index4;
                            if (x[num10 - 1] >= x[index4 - 1])
                            {
                                num10 = 0;
                            }
                            else
                            {
                                var num11 = x[index4 - 1];
                                x[index4 - 1] = x[num10 - 1];
                                x[num10 - 1] = num11;
                                var num12 = y[index4 - 1];
                                y[index4 - 1] = y[num10 - 1];
                                y[num10 - 1] = num12;
                                num10 = index4;
                            }
                        }
                    }

                    --index3;
                } while (index3 >= 1);
            }
        }

        private static void HeapSortDPoints(ref double[] x, ref double[] y, ref double[] d, int n)
        {
            var flag1 = true;
            var flag2 = true;
            for (var index = 1; index <= n - 1; ++index)
            {
                flag1 &= x[index] > x[index - 1];
                flag2 &= x[index] < x[index - 1];
            }

            if (flag1)
                return;
            if (flag2)
            {
                for (var index1 = 0; index1 <= n - 1; ++index1)
                {
                    var index2 = n - 1 - index1;
                    if (index2 <= index1)
                        break;
                    var num1 = x[index1];
                    x[index1] = x[index2];
                    x[index2] = num1;
                    var num2 = y[index1];
                    y[index1] = y[index2];
                    y[index2] = num2;
                    var num3 = d[index1];
                    d[index1] = d[index2];
                    d[index2] = num3;
                }
            }
            else
            {
                if (n == 1)
                    return;
                var num4 = 2;
                do
                {
                    var num5 = num4;
                    while (num5 != 1)
                    {
                        var num6 = num5 / 2;
                        if (x[num6 - 1] >= x[num5 - 1])
                        {
                            num5 = 1;
                        }
                        else
                        {
                            var num7 = x[num6 - 1];
                            x[num6 - 1] = x[num5 - 1];
                            x[num5 - 1] = num7;
                            var num8 = y[num6 - 1];
                            y[num6 - 1] = y[num5 - 1];
                            y[num5 - 1] = num8;
                            var num9 = d[num6 - 1];
                            d[num6 - 1] = d[num5 - 1];
                            d[num5 - 1] = num9;
                            num5 = num6;
                        }
                    }

                    ++num4;
                } while (num4 <= n);

                var index3 = n - 1;
                do
                {
                    var num10 = x[index3];
                    x[index3] = x[0];
                    x[0] = num10;
                    var num11 = y[index3];
                    y[index3] = y[0];
                    y[0] = num11;
                    var num12 = d[index3];
                    d[index3] = d[0];
                    d[0] = num12;
                    var num13 = 1;
                    while (num13 != 0)
                    {
                        var index4 = 2 * num13;
                        if (index4 > index3)
                        {
                            num13 = 0;
                        }
                        else
                        {
                            if (index4 < index3 && x[index4] > x[index4 - 1])
                                ++index4;
                            if (x[num13 - 1] >= x[index4 - 1])
                            {
                                num13 = 0;
                            }
                            else
                            {
                                var num14 = x[index4 - 1];
                                x[index4 - 1] = x[num13 - 1];
                                x[num13 - 1] = num14;
                                var num15 = y[index4 - 1];
                                y[index4 - 1] = y[num13 - 1];
                                y[num13 - 1] = num15;
                                var num16 = d[index4 - 1];
                                d[index4 - 1] = d[num13 - 1];
                                d[num13 - 1] = num16;
                                num13 = index4;
                            }
                        }
                    }

                    --index3;
                } while (index3 >= 1);
            }
        }

        private static void SolveTriDiagonal(
            double[] a,
            double[] b,
            double[] c,
            double[] d,
            int n,
            ref double[] x)
        {
            a = (double[])a.Clone();
            b = (double[])b.Clone();
            c = (double[])c.Clone();
            d = (double[])d.Clone();
            x = new double[n - 1 + 1];
            a[0] = 0.0;
            c[n - 1] = 0.0;
            for (var index = 1; index <= n - 1; ++index)
            {
                var num = a[index] / b[index - 1];
                b[index] = b[index] - num * c[index - 1];
                d[index] = d[index] - num * d[index - 1];
            }

            x[n - 1] = d[n - 1] / b[n - 1];
            for (var index = n - 2; index >= 0; --index)
                x[index] = (d[index] - c[index] * x[index + 1]) / b[index];
        }

        private static double DiffThreePoint(
            double t,
            double x0,
            double f0,
            double x1,
            double f1,
            double x2,
            double f2)
        {
            t -= x0;
            x1 -= x0;
            x2 -= x0;
            var num1 = (f2 - f0 - x2 / x1 * (f1 - f0)) / (CustomMath.Sqr(x2) - x1 * x2);
            var num2 = (f1 - f0 - num1 * CustomMath.Sqr(x1)) / x1;
            return 2.0 * num1 * t + num2;
        }
    }
}