using System;

namespace CanardConfit.NINA.BahtiFocus.Bahtinov {
    public static class blas {
        public static double vectornorm2(ref double[] x, int i1, int i2) {
            var num1 = i2 - i1 + 1;
            if (num1 < 1)
                return 0.0;
            if (num1 == 1)
                return Math.Abs(x[i1]);
            var num2 = 0.0;
            var d = 1.0;
            for (var index = i1; index <= i2; ++index)
                if (x[index] != 0.0) {
                    var num3 = Math.Abs(x[index]);
                    if (num2 < num3) {
                        d = 1.0 + d * CustomMath.Sqr(num2 / num3);
                        num2 = num3;
                    } else {
                        d += CustomMath.Sqr(num3 / num2);
                    }
                }

            return num2 * Math.Sqrt(d);
        }

        public static int vectoridxabsmax(ref double[] x, int i1, int i2) {
            var index1 = i1;
            Math.Abs(x[index1]);
            for (var index2 = i1 + 1; index2 <= i2; ++index2)
                if (Math.Abs(x[index2]) > Math.Abs(x[index1]))
                    index1 = index2;
            return index1;
        }

        public static int columnidxabsmax(ref double[,] x, int i1, int i2, int j) {
            var index1 = i1;
            Math.Abs(x[index1, j]);
            for (var index2 = i1 + 1; index2 <= i2; ++index2)
                if (Math.Abs(x[index2, j]) > Math.Abs(x[index1, j]))
                    index1 = index2;
            return index1;
        }

        public static int rowidxabsmax(ref double[,] x, int j1, int j2, int i) {
            var index1 = j1;
            Math.Abs(x[i, index1]);
            for (var index2 = j1 + 1; index2 <= j2; ++index2)
                if (Math.Abs(x[i, index2]) > Math.Abs(x[i, index1]))
                    index1 = index2;
            return index1;
        }

        public static double upperhessenberg1norm(
            ref double[,] a,
            int i1,
            int i2,
            int j1,
            int j2,
            ref double[] work) {
            for (var index = j1; index <= j2; ++index)
                work[index] = 0.0;
            for (var index1 = i1; index1 <= i2; ++index1)
            for (var index2 = Math.Max(j1, j1 + index1 - i1 - 1); index2 <= j2; ++index2)
                work[index2] = work[index2] + Math.Abs(a[index1, index2]);
            var val1 = 0.0;
            for (var index = j1; index <= j2; ++index)
                val1 = Math.Max(val1, work[index]);
            return val1;
        }

        public static void copymatrix(
            ref double[,] a,
            int is1,
            int is2,
            int js1,
            int js2,
            ref double[,] b,
            int id1,
            int id2,
            int jd1,
            int jd2) {
            if ((is1 > is2) | (js1 > js2))
                return;
            for (var index1 = is1; index1 <= is2; ++index1) {
                var index2 = index1 - is1 + id1;
                var num = js1 - jd1;
                for (var index3 = jd1; index3 <= jd2; ++index3)
                    b[index2, index3] = a[index1, index3 + num];
            }
        }

        public static void inplacetranspose(
            ref double[,] a,
            int i1,
            int i2,
            int j1,
            int j2,
            ref double[] work) {
            if ((i1 > i2) | (j1 > j2))
                return;
            for (var index1 = i1; index1 <= i2 - 1; ++index1) {
                var index2 = j1 + index1 - i1;
                var num1 = index1 + 1;
                var num2 = j1 + num1 - i1;
                var num3 = i2 - index1;
                var num4 = num1 - 1;
                for (var index3 = 1; index3 <= num3; ++index3)
                    work[index3] = a[index3 + num4, index2];
                var num5 = num2 - num1;
                for (var index4 = num1; index4 <= i2; ++index4)
                    a[index4, index2] = a[index1, index4 + num5];
                var num6 = 1 - num2;
                for (var index5 = num2; index5 <= j2; ++index5)
                    a[index1, index5] = work[index5 + num6];
            }
        }

        public static void copyandtranspose(
            ref double[,] a,
            int is1,
            int is2,
            int js1,
            int js2,
            ref double[,] b,
            int id1,
            int id2,
            int jd1,
            int jd2) {
            if ((is1 > is2) | (js1 > js2))
                return;
            for (var index1 = is1; index1 <= is2; ++index1) {
                var index2 = index1 - is1 + jd1;
                var num = js1 - id1;
                for (var index3 = id1; index3 <= id2; ++index3)
                    b[index3, index2] = a[index1, index3 + num];
            }
        }

        public static void matrixvectormultiply(
            ref double[,] a,
            int i1,
            int i2,
            int j1,
            int j2,
            bool trans,
            ref double[] x,
            int ix1,
            int ix2,
            double alpha,
            ref double[] y,
            int iy1,
            int iy2,
            double beta) {
            if (!trans) {
                if ((i1 > i2) | (j1 > j2))
                    return;
                if (beta == 0.0)
                    for (var index = iy1; index <= iy2; ++index)
                        y[index] = 0.0;
                else
                    for (var index = iy1; index <= iy2; ++index)
                        y[index] = beta * y[index];
                for (var index1 = i1; index1 <= i2; ++index1) {
                    var num1 = ix1 - j1;
                    var num2 = 0.0;
                    for (var index2 = j1; index2 <= j2; ++index2)
                        num2 += a[index1, index2] * x[index2 + num1];
                    y[iy1 + index1 - i1] = y[iy1 + index1 - i1] + alpha * num2;
                }
            } else {
                if ((i1 > i2) | (j1 > j2))
                    return;
                if (beta == 0.0)
                    for (var index = iy1; index <= iy2; ++index)
                        y[index] = 0.0;
                else
                    for (var index = iy1; index <= iy2; ++index)
                        y[index] = beta * y[index];
                for (var index3 = i1; index3 <= i2; ++index3) {
                    var num3 = alpha * x[ix1 + index3 - i1];
                    var num4 = j1 - iy1;
                    for (var index4 = iy1; index4 <= iy2; ++index4)
                        y[index4] = y[index4] + num3 * a[index3, index4 + num4];
                }
            }
        }

        public static double pythag2(double x, double y) {
            var val1 = Math.Abs(x);
            var val2 = Math.Abs(y);
            var num1 = Math.Max(val1, val2);
            var num2 = Math.Min(val1, val2);
            return num2 != 0.0 ? num1 * Math.Sqrt(1.0 + CustomMath.Sqr(num2 / num1)) : num1;
        }

        public static void matrixmatrixmultiply(
            ref double[,] a,
            int ai1,
            int ai2,
            int aj1,
            int aj2,
            bool transa,
            ref double[,] b,
            int bi1,
            int bi2,
            int bj1,
            int bj2,
            bool transb,
            double alpha,
            ref double[,] c,
            int ci1,
            int ci2,
            int cj1,
            int cj2,
            double beta,
            ref double[] work) {
            var index1 = 0;
            int val1_1;
            int val2_1;
            if (!transa) {
                val1_1 = ai2 - ai1 + 1;
                val2_1 = aj2 - aj1 + 1;
            } else {
                val1_1 = aj2 - aj1 + 1;
                val2_1 = ai2 - ai1 + 1;
            }

            int val1_2;
            int val2_2;
            if (!transb) {
                val1_2 = bi2 - bi1 + 1;
                val2_2 = bj2 - bj1 + 1;
            } else {
                val1_2 = bj2 - bj1 + 1;
                val2_2 = bi2 - bi1 + 1;
            }

            if ((val1_1 <= 0) | (val2_1 <= 0) | (val1_2 <= 0) | (val2_2 <= 0))
                return;
            var num1 = val1_1;
            var val2_3 = Math.Max(val1_1, val2_1);
            var index2 = Math.Max(Math.Max(val1_2, val2_3), val2_2);
            work[1] = 0.0;
            work[index2] = 0.0;
            if (beta == 0.0)
                for (var index3 = ci1; index3 <= ci2; ++index3)
                for (var index4 = cj1; index4 <= cj2; ++index4)
                    c[index3, index4] = 0.0;
            else
                for (var index5 = ci1; index5 <= ci2; ++index5)
                for (var index6 = cj1; index6 <= cj2; ++index6)
                    c[index5, index6] = beta * c[index5, index6];
            if (!transa & !transb) {
                for (var index7 = ai1; index7 <= ai2; ++index7)
                for (var index8 = bi1; index8 <= bi2; ++index8) {
                    var num2 = alpha * a[index7, aj1 + index8 - bi1];
                    var index9 = ci1 + index7 - ai1;
                    var num3 = bj1 - cj1;
                    for (var index10 = cj1; index10 <= cj2; ++index10)
                        c[index9, index10] = c[index9, index10] + num2 * b[index8, index10 + num3];
                }
            } else if (!transa & transb) {
                if (val1_1 * val2_1 < val1_2 * val2_2)
                    for (var index11 = bi1; index11 <= bi2; ++index11)
                    for (var index12 = ai1; index12 <= ai2; ++index12) {
                        var num4 = bj1 - aj1;
                        var num5 = 0.0;
                        for (var index13 = aj1; index13 <= aj2; ++index13)
                            num5 += a[index12, index13] * b[index11, index13 + num4];
                        c[ci1 + index12 - ai1, cj1 + index11 - bi1] =
                            c[ci1 + index12 - ai1, cj1 + index11 - bi1] + alpha * num5;
                    }
                else
                    for (var index14 = ai1; index14 <= ai2; ++index14)
                    for (var index15 = bi1; index15 <= bi2; ++index15) {
                        var num6 = bj1 - aj1;
                        var num7 = 0.0;
                        for (var index16 = aj1; index16 <= aj2; ++index16)
                            num7 += a[index14, index16] * b[index15, index16 + num6];
                        c[ci1 + index14 - ai1, cj1 + index15 - bi1] =
                            c[ci1 + index14 - ai1, cj1 + index15 - bi1] + alpha * num7;
                    }
            } else if (transa & !transb) {
                for (var index17 = aj1; index17 <= aj2; ++index17)
                for (var index18 = bi1; index18 <= bi2; ++index18) {
                    var num8 = alpha * a[ai1 + index18 - bi1, index17];
                    var index19 = ci1 + index17 - aj1;
                    var num9 = bj1 - cj1;
                    for (var index20 = cj1; index20 <= cj2; ++index20)
                        c[index19, index20] = c[index19, index20] + num8 * b[index18, index20 + num9];
                }
            } else {
                if (!(transa & transb))
                    return;
                if (val1_1 * val2_1 < val1_2 * val2_2)
                    for (var index21 = bi1; index21 <= bi2; ++index21) {
                        for (var index22 = 1; index22 <= num1; ++index22)
                            work[index22] = 0.0;
                        for (var index23 = ai1; index23 <= ai2; ++index23) {
                            var num10 = alpha * b[index21, bj1 + index23 - ai1];
                            index1 = cj1 + index21 - bi1;
                            var num11 = aj1 - 1;
                            for (var index24 = 1; index24 <= num1; ++index24)
                                work[index24] = work[index24] + num10 * a[index23, index24 + num11];
                        }

                        var num12 = 1 - ci1;
                        for (var index25 = ci1; index25 <= ci2; ++index25)
                            c[index25, index1] = c[index25, index1] + work[index25 + num12];
                    }
                else
                    for (var index26 = aj1; index26 <= aj2; ++index26) {
                        var num13 = ai2 - ai1 + 1;
                        var num14 = ai1 - 1;
                        for (var index27 = 1; index27 <= num13; ++index27)
                            work[index27] = a[index27 + num14, index26];
                        for (var index28 = bi1; index28 <= bi2; ++index28) {
                            var num15 = bj1 - 1;
                            var num16 = 0.0;
                            for (var index29 = 1; index29 <= num13; ++index29)
                                num16 += work[index29] * b[index28, index29 + num15];
                            c[ci1 + index26 - aj1, cj1 + index28 - bi1] =
                                c[ci1 + index26 - aj1, cj1 + index28 - bi1] + alpha * num16;
                        }
                    }
            }
        }
    }
}