using System;

namespace CanardConfit.NINA.BahtiFocus.Bahtinov {
    public static class Reflections {
        public static void GenerateReflection(ref double[] x, int n, ref double tau) {
            if (n <= 1) {
                tau = 0.0;
            } else {
                var num1 = x[1];
                var val2 = 0.0;
                for (var index = 2; index <= n; ++index)
                    val2 = Math.Max(Math.Abs(x[index]), val2);
                var d = 0.0;
                if (val2 != 0.0) {
                    for (var index = 2; index <= n; ++index)
                        d += CustomMath.Sqr(x[index] / val2);
                    d = Math.Sqrt(d) * val2;
                }

                if (d == 0.0) {
                    tau = 0.0;
                } else {
                    var num2 = Math.Max(Math.Abs(num1), Math.Abs(d));
                    var num3 = -(num2 * Math.Sqrt(CustomMath.Sqr(num1 / num2) + CustomMath.Sqr(d / num2)));
                    if (num1 < 0.0)
                        num3 = -num3;
                    tau = (num3 - num1) / num3;
                    var num4 = 1.0 / (num1 - num3);
                    for (var index = 2; index <= n; ++index)
                        x[index] = num4 * x[index];
                    x[1] = num3;
                }
            }
        }

        public static void ApplyReflectionFromTheLeft(ref double[,] c, double tau, ref double[] v, int m1, int m2,
            int n1, int n2, ref double[] work) {
            if (tau == 0.0 || n1 > n2 || m1 > m2)
                return;
            for (var index = n1; index <= n2; ++index)
                work[index] = 0.0;
            for (var index1 = m1; index1 <= m2; ++index1) {
                var num = v[index1 + 1 - m1];
                for (var index2 = n1; index2 <= n2; ++index2)
                    work[index2] += num * c[index1, index2];
            }

            for (var index3 = m1; index3 <= m2; ++index3) {
                var num = v[index3 - m1 + 1] * tau;
                for (var index4 = n1; index4 <= n2; ++index4)
                    c[index3, index4] -= num * work[index4];
            }
        }

        public static void ApplyReflectionFromTheRight(ref double[,] c, double tau, ref double[] v, int m1, int m2,
            int n1, int n2, ref double[] work) {
            if (tau == 0.0 || n1 > n2 || m1 > m2)
                return;
            for (var index1 = m1; index1 <= m2; ++index1) {
                var num1 = 1 - n1;
                var num2 = 0.0;
                for (var index2 = n1; index2 <= n2; ++index2)
                    num2 += c[index1, index2] * v[index2 + num1];
                work[index1] = num2;
            }

            for (var index3 = m1; index3 <= m2; ++index3) {
                var num3 = work[index3] * tau;
                var num4 = 1 - n1;
                for (var index4 = n1; index4 <= n2; ++index4)
                    c[index3, index4] -= num3 * v[index4 + num4];
            }
        }
    }
}