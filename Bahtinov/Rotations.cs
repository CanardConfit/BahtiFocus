using System;

namespace CanardConfit.NINA.BahtiFocus.Bahtinov {
    public static class Rotations {
        public static void ApplyRotationsFromTheLeft(bool isforward, int m1, int m2, int n1, int n2, ref double[] c,
            ref double[] s, ref double[,] a, ref double[] work) {
            if (m1 > m2 || n1 > n2)
                return;
            if (isforward) {
                if (n1 != n2)
                    for (var index1 = m1; index1 <= m2 - 1; ++index1) {
                        var num1 = c[index1 - m1 + 1];
                        var num2 = s[index1 - m1 + 1];
                        if (num1 != 1.0 || num2 != 0.0) {
                            var index2 = index1 + 1;
                            for (var index3 = n1; index3 <= n2; ++index3)
                                work[index3] = num1 * a[index2, index3];
                            for (var index4 = n1; index4 <= n2; ++index4)
                                work[index4] = work[index4] - num2 * a[index1, index4];
                            for (var index5 = n1; index5 <= n2; ++index5)
                                a[index1, index5] = num1 * a[index1, index5];
                            for (var index6 = n1; index6 <= n2; ++index6)
                                a[index1, index6] = a[index1, index6] + num2 * a[index2, index6];
                            for (var index7 = n1; index7 <= n2; ++index7)
                                a[index2, index7] = work[index7];
                        }
                    }
                else
                    for (var index = m1; index <= m2 - 1; ++index) {
                        var num3 = c[index - m1 + 1];
                        var num4 = s[index - m1 + 1];
                        if (num3 != 1.0 || num4 != 0.0) {
                            var num5 = a[index + 1, n1];
                            a[index + 1, n1] = num3 * num5 - num4 * a[index, n1];
                            a[index, n1] = num4 * num5 + num3 * a[index, n1];
                        }
                    }
            } else if (n1 != n2) {
                for (var index8 = m2 - 1; index8 >= m1; --index8) {
                    var num6 = c[index8 - m1 + 1];
                    var num7 = s[index8 - m1 + 1];
                    if (num6 != 1.0 || num7 != 0.0) {
                        var index9 = index8 + 1;
                        for (var index10 = n1; index10 <= n2; ++index10)
                            work[index10] = num6 * a[index9, index10];
                        for (var index11 = n1; index11 <= n2; ++index11)
                            work[index11] = work[index11] - num7 * a[index8, index11];
                        for (var index12 = n1; index12 <= n2; ++index12)
                            a[index8, index12] = num6 * a[index8, index12];
                        for (var index13 = n1; index13 <= n2; ++index13)
                            a[index8, index13] = a[index8, index13] + num7 * a[index9, index13];
                        for (var index14 = n1; index14 <= n2; ++index14)
                            a[index9, index14] = work[index14];
                    }
                }
            } else {
                for (var index = m2 - 1; index >= m1; --index) {
                    var num8 = c[index - m1 + 1];
                    var num9 = s[index - m1 + 1];
                    if (num8 != 1.0 || num9 != 0.0) {
                        var num10 = a[index + 1, n1];
                        a[index + 1, n1] = num8 * num10 - num9 * a[index, n1];
                        a[index, n1] = num9 * num10 + num8 * a[index, n1];
                    }
                }
            }
        }

        public static void ApplyRotationsFromTheRight(bool isforward, int m1, int m2, int n1, int n2, ref double[] c,
            ref double[] s, ref double[,] a, ref double[] work) {
            if (isforward) {
                if (m1 != m2)
                    for (var index1 = n1; index1 <= n2 - 1; ++index1) {
                        var num1 = c[index1 - n1 + 1];
                        var num2 = s[index1 - n1 + 1];
                        if (num1 != 1.0 || num2 != 0.0) {
                            var index2 = index1 + 1;
                            for (var index3 = m1; index3 <= m2; ++index3)
                                work[index3] = num1 * a[index3, index2];
                            for (var index4 = m1; index4 <= m2; ++index4)
                                work[index4] = work[index4] - num2 * a[index4, index1];
                            for (var index5 = m1; index5 <= m2; ++index5)
                                a[index5, index1] = num1 * a[index5, index1];
                            for (var index6 = m1; index6 <= m2; ++index6)
                                a[index6, index1] = a[index6, index1] + num2 * a[index6, index2];
                            for (var index7 = m1; index7 <= m2; ++index7)
                                a[index7, index2] = work[index7];
                        }
                    }
                else
                    for (var index = n1; index <= n2 - 1; ++index) {
                        var num3 = c[index - n1 + 1];
                        var num4 = s[index - n1 + 1];
                        if (num3 != 1.0 || num4 != 0.0) {
                            var num5 = a[m1, index + 1];
                            a[m1, index + 1] = num3 * num5 - num4 * a[m1, index];
                            a[m1, index] = num4 * num5 + num3 * a[m1, index];
                        }
                    }
            } else if (m1 != m2) {
                for (var index8 = n2 - 1; index8 >= n1; --index8) {
                    var num6 = c[index8 - n1 + 1];
                    var num7 = s[index8 - n1 + 1];
                    if (num6 != 1.0 || num7 != 0.0) {
                        var index9 = index8 + 1;
                        for (var index10 = m1; index10 <= m2; ++index10)
                            work[index10] = num6 * a[index10, index9];
                        for (var index11 = m1; index11 <= m2; ++index11)
                            work[index11] = work[index11] - num7 * a[index11, index8];
                        for (var index12 = m1; index12 <= m2; ++index12)
                            a[index12, index8] = num6 * a[index12, index8];
                        for (var index13 = m1; index13 <= m2; ++index13)
                            a[index13, index8] = a[index13, index8] + num7 * a[index13, index9];
                        for (var index14 = m1; index14 <= m2; ++index14)
                            a[index14, index9] = work[index14];
                    }
                }
            } else {
                for (var index = n2 - 1; index >= n1; --index) {
                    var num8 = c[index - n1 + 1];
                    var num9 = s[index - n1 + 1];
                    if (num8 != 1.0 || num9 != 0.0) {
                        var num10 = a[m1, index + 1];
                        a[m1, index + 1] = num8 * num10 - num9 * a[m1, index];
                        a[m1, index] = num9 * num10 + num8 * a[m1, index];
                    }
                }
            }
        }

        public static void GenerateRotation(double f, double g, ref double cs, ref double sn, ref double r) {
            if (g == 0.0) {
                cs = 1.0;
                sn = 0.0;
                r = f;
            } else if (f == 0.0) {
                cs = 0.0;
                sn = 1.0;
                r = g;
            } else {
                var X1 = f;
                var X2 = g;
                r = Math.Sqrt(CustomMath.Sqr(X1) + CustomMath.Sqr(X2));
                cs = X1 / r;
                sn = X2 / r;
                if (!((Math.Abs(f) > Math.Abs(g)) & (cs < 0.0)))
                    return;
                cs = -cs;
                sn = -sn;
                r = -r;
            }
        }
    }
}