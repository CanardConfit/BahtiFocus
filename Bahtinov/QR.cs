using System;

namespace CanardConfit.NINA.BahtiFocus.Bahtinov {
    public static class QR {
        public static void RMatrixQR(ref double[,] a, int m, int n, ref double[] tau) {
            var tau1 = 0.0;
            if (m <= 0 || n <= 0)
                return;
            var num1 = Math.Min(m, n);
            var work = new double[n - 1 + 1];
            var numArray3 = new double[m + 1];
            tau = new double[num1 - 1 + 1];
            var num2 = num1;
            for (var m1 = 0; m1 <= num2 - 1; ++m1) {
                var num3 = m1 - 1;
                for (var index = 1; index <= m - m1; ++index)
                    numArray3[index] = a[index + num3, m1];
                Reflections.GenerateReflection(ref numArray3, m - m1, ref tau1);
                tau[m1] = tau1;
                var num4 = 1 - m1;
                for (var index = m1; index <= m - 1; ++index)
                    a[index, m1] = numArray3[index + num4];
                numArray3[1] = 1.0;
                if (m1 < n)
                    Reflections.ApplyReflectionFromTheLeft(ref a, tau[m1], ref numArray3, m1, m - 1, m1 + 1, n - 1,
                        ref work);
            }
        }

        public static void RMatrixQRUnPackQ(
            ref double[,] a,
            int m,
            int n,
            ref double[] tau,
            int qcolumns,
            ref double[,] q) {
            if (m <= 0 || n <= 0 || qcolumns <= 0)
                return;
            var num1 = Math.Min(Math.Min(m, n), qcolumns);
            q = new double[m - 1 + 1, qcolumns - 1 + 1];
            var v = new double[m + 1];
            var work = new double[qcolumns - 1 + 1];
            for (var index1 = 0; index1 <= m - 1; ++index1)
            for (var index2 = 0; index2 <= qcolumns - 1; ++index2)
                q[index1, index2] = index1 != index2 ? 0.0 : 1.0;
            for (var m1 = num1 - 1; m1 >= 0; --m1) {
                var num2 = m1 - 1;
                for (var index = 1; index <= m - m1; ++index)
                    v[index] = a[index + num2, m1];
                v[1] = 1.0;
                Reflections.ApplyReflectionFromTheLeft(ref q, tau[m1], ref v, m1, m - 1, 0, qcolumns - 1, ref work);
            }
        }

        public static void RMatrixQRUnPackR(ref double[,] a, int m, int n, ref double[,] r) {
            if (m <= 0 || n <= 0)
                return;
            var num = Math.Min(m, n);
            r = new double[m - 1 + 1, n - 1 + 1];
            for (var index = 0; index <= n - 1; ++index)
                r[0, index] = 0.0;
            for (var index1 = 1; index1 <= m - 1; ++index1)
            for (var index2 = 0; index2 <= n - 1; ++index2)
                r[index1, index2] = r[0, index2];
            for (var index3 = 0; index3 <= num - 1; ++index3)
            for (var index4 = index3; index4 <= n - 1; ++index4)
                r[index3, index4] = a[index3, index4];
        }

        public static void QRDecomposition(ref double[,] a, int m, int n, ref double[] tau) {
            var tau1 = 0.0;
            var num1 = Math.Min(m, n);
            var work = new double[n + 1];
            var numArray3 = new double[m + 1];
            tau = new double[num1 + 1];
            var num2 = Math.Min(m, n);
            for (var m1 = 1; m1 <= num2; ++m1) {
                var n1 = m - m1 + 1;
                var num3 = m1 - 1;
                for (var index = 1; index <= n1; ++index)
                    numArray3[index] = a[index + num3, m1];
                Reflections.GenerateReflection(ref numArray3, n1, ref tau1);
                tau[m1] = tau1;
                var num4 = 1 - m1;
                for (var index = m1; index <= m; ++index)
                    a[index, m1] = numArray3[index + num4];
                numArray3[1] = 1.0;
                if (m1 < n)
                    Reflections.ApplyReflectionFromTheLeft(ref a, tau[m1], ref numArray3, m1, m, m1 + 1, n, ref work);
            }
        }

        public static void UnPackQFromQR(
            ref double[,] a,
            int m,
            int n,
            ref double[] tau,
            int qcolumns,
            ref double[,] q) {
            if (m == 0 || n == 0 || qcolumns == 0)
                return;
            var num1 = Math.Min(Math.Min(m, n), qcolumns);
            q = new double[m + 1, qcolumns + 1];
            var v = new double[m + 1];
            var work = new double[qcolumns + 1];
            for (var index1 = 1; index1 <= m; ++index1)
            for (var index2 = 1; index2 <= qcolumns; ++index2)
                q[index1, index2] = index1 != index2 ? 0.0 : 1.0;
            for (var m1 = num1; m1 >= 1; --m1) {
                var num2 = m - m1 + 1;
                var num3 = m1 - 1;
                for (var index = 1; index <= num2; ++index)
                    v[index] = a[index + num3, m1];
                v[1] = 1.0;
                Reflections.ApplyReflectionFromTheLeft(ref q, tau[m1], ref v, m1, m, 1, qcolumns, ref work);
            }
        }

        public static void QRDecompositionUnPacked(
            double[,] a,
            int m,
            int n,
            ref double[,] q,
            ref double[,] r) {
            var tau = new double[0];
            a = (double[,])a.Clone();
            var num = Math.Min(m, n);
            if (n <= 0)
                return;
            q = new double[m + 1, m + 1];
            r = new double[m + 1, n + 1];
            QRDecomposition(ref a, m, n, ref tau);
            for (var index = 1; index <= n; ++index)
                r[1, index] = 0.0;
            for (var index1 = 2; index1 <= m; ++index1)
            for (var index2 = 1; index2 <= n; ++index2)
                r[index1, index2] = r[1, index2];
            for (var index3 = 1; index3 <= num; ++index3)
            for (var index4 = index3; index4 <= n; ++index4)
                r[index3, index4] = a[index3, index4];
            UnPackQFromQR(ref a, m, n, ref tau, m, ref q);
        }
    }
}