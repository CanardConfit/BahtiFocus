using System;

namespace CanardConfit.NINA.BahtiFocus.Bahtinov {
    public static class bidiagonal
    {
        public static void rmatrixbd(
            ref double[,] a,
            int m,
            int n,
            ref double[] tauq,
            ref double[] taup)
        {
            var numArray1 = new double[0];
            var numArray2 = new double[0];
            var tau = 0.0;
            if ((n <= 0) | (m <= 0))
                return;
            Math.Min(m, n);
            var num1 = Math.Max(m, n);
            var work = new double[num1 + 1];
            var numArray3 = new double[num1 + 1];
            if (m >= n)
            {
                tauq = new double[n - 1 + 1];
                taup = new double[n - 1 + 1];
            }
            else
            {
                tauq = new double[m - 1 + 1];
                taup = new double[m - 1 + 1];
            }

            if (m >= n)
                for (var m1 = 0; m1 <= n - 1; ++m1)
                {
                    var num2 = m1 - 1;
                    for (var index = 1; index <= m - m1; ++index)
                        numArray3[index] = a[index + num2, m1];
                    Reflections.GenerateReflection(ref numArray3, m - m1, ref tau);
                    tauq[m1] = tau;
                    var num3 = 1 - m1;
                    for (var index = m1; index <= m - 1; ++index)
                        a[index, m1] = numArray3[index + num3];
                    numArray3[1] = 1.0;
                    Reflections.ApplyReflectionFromTheLeft(ref a, tau, ref numArray3, m1, m - 1, m1 + 1, n - 1, ref work);
                    if (m1 < n - 1)
                    {
                        var num4 = m1 + 1 - 1;
                        for (var index = 1; index <= n - m1 - 1; ++index)
                            numArray3[index] = a[m1, index + num4];
                        Reflections.GenerateReflection(ref numArray3, n - 1 - m1, ref tau);
                        taup[m1] = tau;
                        var num5 = 1 - (m1 + 1);
                        for (var index = m1 + 1; index <= n - 1; ++index)
                            a[m1, index] = numArray3[index + num5];
                        numArray3[1] = 1.0;
                        Reflections.ApplyReflectionFromTheRight(ref a, tau, ref numArray3, m1 + 1, m - 1, m1 + 1, n - 1,
                            ref work);
                    }
                    else
                    {
                        taup[m1] = 0.0;
                    }
                }
            else
                for (var n1 = 0; n1 <= m - 1; ++n1)
                {
                    var num6 = n1 - 1;
                    for (var index = 1; index <= n - n1; ++index)
                        numArray3[index] = a[n1, index + num6];
                    Reflections.GenerateReflection(ref numArray3, n - n1, ref tau);
                    taup[n1] = tau;
                    var num7 = 1 - n1;
                    for (var index = n1; index <= n - 1; ++index)
                        a[n1, index] = numArray3[index + num7];
                    numArray3[1] = 1.0;
                    Reflections.ApplyReflectionFromTheRight(ref a, tau, ref numArray3, n1 + 1, m - 1, n1, n - 1, ref work);
                    if (n1 < m - 1)
                    {
                        var num8 = n1 + 1 - 1;
                        for (var index = 1; index <= m - 1 - n1; ++index)
                            numArray3[index] = a[index + num8, n1];
                        Reflections.GenerateReflection(ref numArray3, m - 1 - n1, ref tau);
                        tauq[n1] = tau;
                        var num9 = 1 - (n1 + 1);
                        for (var index = n1 + 1; index <= m - 1; ++index)
                            a[index, n1] = numArray3[index + num9];
                        numArray3[1] = 1.0;
                        Reflections.ApplyReflectionFromTheLeft(ref a, tau, ref numArray3, n1 + 1, m - 1, n1 + 1, n - 1,
                            ref work);
                    }
                    else
                    {
                        tauq[n1] = 0.0;
                    }
                }
        }

        public static void rmatrixbdunpackq(
            ref double[,] qp,
            int m,
            int n,
            ref double[] tauq,
            int qcolumns,
            ref double[,] q)
        {
            if ((m == 0) | (n == 0) | (qcolumns == 0))
                return;
            q = new double[m - 1 + 1, qcolumns - 1 + 1];
            for (var index1 = 0; index1 <= m - 1; ++index1)
            for (var index2 = 0; index2 <= qcolumns - 1; ++index2)
                q[index1, index2] = index1 != index2 ? 0.0 : 1.0;
            rmatrixbdmultiplybyq(ref qp, m, n, ref tauq, ref q, m, qcolumns, false, false);
        }

        public static void rmatrixbdmultiplybyq(
            ref double[,] qp,
            int m,
            int n,
            ref double[] tauq,
            ref double[,] z,
            int zrows,
            int zcolumns,
            bool fromtheright,
            bool dotranspose)
        {
            var numArray1 = new double[0];
            var numArray2 = new double[0];
            if ((m <= 0) | (n <= 0) | (zrows <= 0) | (zcolumns <= 0))
                return;
            var num1 = Math.Max(Math.Max(Math.Max(m, n), zrows), zcolumns);
            var v = new double[num1 + 1];
            var work = new double[num1 + 1];
            if (m >= n)
            {
                int num2;
                int num3;
                int num4;
                if (fromtheright)
                {
                    num2 = 0;
                    num3 = n - 1;
                    num4 = 1;
                }
                else
                {
                    num2 = n - 1;
                    num3 = 0;
                    num4 = -1;
                }

                if (dotranspose)
                {
                    var num5 = num2;
                    num2 = num3;
                    num3 = num5;
                    num4 = -num4;
                }

                var index1 = num2;
                do
                {
                    var num6 = index1 - 1;
                    for (var index2 = 1; index2 <= m - index1; ++index2)
                        v[index2] = qp[index2 + num6, index1];
                    v[1] = 1.0;
                    if (fromtheright)
                        Reflections.ApplyReflectionFromTheRight(ref z, tauq[index1], ref v, 0, zrows - 1, index1, m - 1,
                            ref work);
                    else
                        Reflections.ApplyReflectionFromTheLeft(ref z, tauq[index1], ref v, index1, m - 1, 0, zcolumns - 1,
                            ref work);
                    index1 += num4;
                } while (index1 != num3 + num4);
            }
            else
            {
                int num7;
                int num8;
                int num9;
                if (fromtheright)
                {
                    num7 = 0;
                    num8 = m - 2;
                    num9 = 1;
                }
                else
                {
                    num7 = m - 2;
                    num8 = 0;
                    num9 = -1;
                }

                if (dotranspose)
                {
                    var num10 = num7;
                    num7 = num8;
                    num8 = num10;
                    num9 = -num9;
                }

                if (m - 1 <= 0)
                    return;
                var index3 = num7;
                do
                {
                    var num11 = index3 + 1 - 1;
                    for (var index4 = 1; index4 <= m - index3 - 1; ++index4)
                        v[index4] = qp[index4 + num11, index3];
                    v[1] = 1.0;
                    if (fromtheright)
                        Reflections.ApplyReflectionFromTheRight(ref z, tauq[index3], ref v, 0, zrows - 1, index3 + 1, m - 1,
                            ref work);
                    else
                        Reflections.ApplyReflectionFromTheLeft(ref z, tauq[index3], ref v, index3 + 1, m - 1, 0,
                            zcolumns - 1, ref work);
                    index3 += num9;
                } while (index3 != num8 + num9);
            }
        }

        public static void rmatrixbdunpackpt(
            ref double[,] qp,
            int m,
            int n,
            ref double[] taup,
            int ptrows,
            ref double[,] pt)
        {
            if ((m == 0) | (n == 0) | (ptrows == 0))
                return;
            pt = new double[ptrows - 1 + 1, n - 1 + 1];
            for (var index1 = 0; index1 <= ptrows - 1; ++index1)
            for (var index2 = 0; index2 <= n - 1; ++index2)
                pt[index1, index2] = index1 != index2 ? 0.0 : 1.0;
            rmatrixbdmultiplybyp(ref qp, m, n, ref taup, ref pt, ptrows, n, true, true);
        }

        public static void rmatrixbdmultiplybyp(
            ref double[,] qp,
            int m,
            int n,
            ref double[] taup,
            ref double[,] z,
            int zrows,
            int zcolumns,
            bool fromtheright,
            bool dotranspose)
        {
            var numArray1 = new double[0];
            var numArray2 = new double[0];
            if ((m <= 0) | (n <= 0) | (zrows <= 0) | (zcolumns <= 0))
                return;
            var num1 = Math.Max(Math.Max(Math.Max(m, n), zrows), zcolumns);
            numArray1 = new double[num1 + 1];
            numArray2 = new double[num1 + 1];
            var v = new double[num1 + 1];
            var work = new double[num1 + 1];
            if (m >= n)
            {
                int num2;
                int num3;
                int num4;
                if (fromtheright)
                {
                    num2 = n - 2;
                    num3 = 0;
                    num4 = -1;
                }
                else
                {
                    num2 = 0;
                    num3 = n - 2;
                    num4 = 1;
                }

                if (!dotranspose)
                {
                    var num5 = num2;
                    num2 = num3;
                    num3 = num5;
                    num4 = -num4;
                }

                if (n - 1 <= 0)
                    return;
                var index1 = num2;
                do
                {
                    var num6 = index1 + 1 - 1;
                    for (var index2 = 1; index2 <= n - 1 - index1; ++index2)
                        v[index2] = qp[index1, index2 + num6];
                    v[1] = 1.0;
                    if (fromtheright)
                        Reflections.ApplyReflectionFromTheRight(ref z, taup[index1], ref v, 0, zrows - 1, index1 + 1, n - 1,
                            ref work);
                    else
                        Reflections.ApplyReflectionFromTheLeft(ref z, taup[index1], ref v, index1 + 1, n - 1, 0,
                            zcolumns - 1, ref work);
                    index1 += num4;
                } while (index1 != num3 + num4);
            }
            else
            {
                int num7;
                int num8;
                int num9;
                if (fromtheright)
                {
                    num7 = m - 1;
                    num8 = 0;
                    num9 = -1;
                }
                else
                {
                    num7 = 0;
                    num8 = m - 1;
                    num9 = 1;
                }

                if (!dotranspose)
                {
                    var num10 = num7;
                    num7 = num8;
                    num8 = num10;
                    num9 = -num9;
                }

                var index3 = num7;
                do
                {
                    var num11 = index3 - 1;
                    for (var index4 = 1; index4 <= n - index3; ++index4)
                        v[index4] = qp[index3, index4 + num11];
                    v[1] = 1.0;
                    if (fromtheright)
                        Reflections.ApplyReflectionFromTheRight(ref z, taup[index3], ref v, 0, zrows - 1, index3, n - 1,
                            ref work);
                    else
                        Reflections.ApplyReflectionFromTheLeft(ref z, taup[index3], ref v, index3, n - 1, 0, zcolumns - 1,
                            ref work);
                    index3 += num9;
                } while (index3 != num8 + num9);
            }
        }

        public static void rmatrixbdunpackdiagonals(
            ref double[,] b,
            int m,
            int n,
            ref bool isupper,
            ref double[] d,
            ref double[] e)
        {
            isupper = m >= n;
            if ((m <= 0) | (n <= 0))
                return;
            if (isupper)
            {
                d = new double[n - 1 + 1];
                e = new double[n - 1 + 1];
                for (var index = 0; index <= n - 2; ++index)
                {
                    d[index] = b[index, index];
                    e[index] = b[index, index + 1];
                }

                d[n - 1] = b[n - 1, n - 1];
            }
            else
            {
                d = new double[m - 1 + 1];
                e = new double[m - 1 + 1];
                for (var index = 0; index <= m - 2; ++index)
                {
                    d[index] = b[index, index];
                    e[index] = b[index + 1, index];
                }

                d[m - 1] = b[m - 1, m - 1];
            }
        }

        public static void tobidiagonal(
            ref double[,] a,
            int m,
            int n,
            ref double[] tauq,
            ref double[] taup)
        {
            var numArray1 = new double[0];
            var numArray2 = new double[0];
            var tau = 0.0;
            var num1 = Math.Min(m, n);
            var num2 = Math.Max(m, n);
            var work = new double[num2 + 1];
            var numArray3 = new double[num2 + 1];
            taup = new double[num1 + 1];
            tauq = new double[num1 + 1];
            if (m >= n)
                for (var m1 = 1; m1 <= n; ++m1)
                {
                    var n1 = m - m1 + 1;
                    var num3 = m1 - 1;
                    for (var index = 1; index <= n1; ++index)
                        numArray3[index] = a[index + num3, m1];
                    Reflections.GenerateReflection(ref numArray3, n1, ref tau);
                    tauq[m1] = tau;
                    var num4 = 1 - m1;
                    for (var index = m1; index <= m; ++index)
                        a[index, m1] = numArray3[index + num4];
                    numArray3[1] = 1.0;
                    Reflections.ApplyReflectionFromTheLeft(ref a, tau, ref numArray3, m1, m, m1 + 1, n, ref work);
                    if (m1 < n)
                    {
                        var n2 = n - m1;
                        var num5 = m1 + 1;
                        var num6 = num5 - 1;
                        for (var index = 1; index <= n2; ++index)
                            numArray3[index] = a[m1, index + num6];
                        Reflections.GenerateReflection(ref numArray3, n2, ref tau);
                        taup[m1] = tau;
                        var num7 = 1 - num5;
                        for (var index = num5; index <= n; ++index)
                            a[m1, index] = numArray3[index + num7];
                        numArray3[1] = 1.0;
                        Reflections.ApplyReflectionFromTheRight(ref a, tau, ref numArray3, m1 + 1, m, m1 + 1, n, ref work);
                    }
                    else
                    {
                        taup[m1] = 0.0;
                    }
                }
            else
                for (var n1 = 1; n1 <= m; ++n1)
                {
                    var n3 = n - n1 + 1;
                    var num8 = n1 - 1;
                    for (var index = 1; index <= n3; ++index)
                        numArray3[index] = a[n1, index + num8];
                    Reflections.GenerateReflection(ref numArray3, n3, ref tau);
                    taup[n1] = tau;
                    var num9 = 1 - n1;
                    for (var index = n1; index <= n; ++index)
                        a[n1, index] = numArray3[index + num9];
                    numArray3[1] = 1.0;
                    Reflections.ApplyReflectionFromTheRight(ref a, tau, ref numArray3, n1 + 1, m, n1, n, ref work);
                    if (n1 < m)
                    {
                        var n4 = m - n1;
                        var num10 = n1 + 1;
                        var num11 = num10 - 1;
                        for (var index = 1; index <= n4; ++index)
                            numArray3[index] = a[index + num11, n1];
                        Reflections.GenerateReflection(ref numArray3, n4, ref tau);
                        tauq[n1] = tau;
                        var num12 = 1 - num10;
                        for (var index = num10; index <= m; ++index)
                            a[index, n1] = numArray3[index + num12];
                        numArray3[1] = 1.0;
                        Reflections.ApplyReflectionFromTheLeft(ref a, tau, ref numArray3, n1 + 1, m, n1 + 1, n, ref work);
                    }
                    else
                    {
                        tauq[n1] = 0.0;
                    }
                }
        }

        public static void unpackqfrombidiagonal(
            ref double[,] qp,
            int m,
            int n,
            ref double[] tauq,
            int qcolumns,
            ref double[,] q)
        {
            var numArray1 = new double[0];
            var numArray2 = new double[0];
            if ((m == 0) | (n == 0) | (qcolumns == 0))
                return;
            q = new double[m + 1, qcolumns + 1];
            var v = new double[m + 1];
            var work = new double[qcolumns + 1];
            for (var index1 = 1; index1 <= m; ++index1)
            for (var index2 = 1; index2 <= qcolumns; ++index2)
                q[index1, index2] = index1 != index2 ? 0.0 : 1.0;
            if (m >= n)
                for (var m1 = Math.Min(n, qcolumns); m1 >= 1; --m1)
                {
                    var num1 = m - m1 + 1;
                    var num2 = m1 - 1;
                    for (var index = 1; index <= num1; ++index)
                        v[index] = qp[index + num2, m1];
                    v[1] = 1.0;
                    Reflections.ApplyReflectionFromTheLeft(ref q, tauq[m1], ref v, m1, m, 1, qcolumns, ref work);
                }
            else
                for (var index3 = Math.Min(m - 1, qcolumns - 1); index3 >= 1; --index3)
                {
                    var num3 = m - index3;
                    var num4 = index3 + 1 - 1;
                    for (var index4 = 1; index4 <= num3; ++index4)
                        v[index4] = qp[index4 + num4, index3];
                    v[1] = 1.0;
                    Reflections.ApplyReflectionFromTheLeft(ref q, tauq[index3], ref v, index3 + 1, m, 1, qcolumns,
                        ref work);
                }
        }

        public static void multiplybyqfrombidiagonal(
            ref double[,] qp,
            int m,
            int n,
            ref double[] tauq,
            ref double[,] z,
            int zrows,
            int zcolumns,
            bool fromtheright,
            bool dotranspose)
        {
            var numArray1 = new double[0];
            var numArray2 = new double[0];
            if ((m <= 0) | (n <= 0) | (zrows <= 0) | (zcolumns <= 0))
                return;
            var num1 = Math.Max(Math.Max(Math.Max(m, n), zrows), zcolumns);
            var v = new double[num1 + 1];
            var work = new double[num1 + 1];
            if (m >= n)
            {
                int num2;
                int num3;
                int num4;
                if (fromtheright)
                {
                    num2 = 1;
                    num3 = n;
                    num4 = 1;
                }
                else
                {
                    num2 = n;
                    num3 = 1;
                    num4 = -1;
                }

                if (dotranspose)
                {
                    var num5 = num2;
                    num2 = num3;
                    num3 = num5;
                    num4 = -num4;
                }

                var index1 = num2;
                do
                {
                    var num6 = m - index1 + 1;
                    var num7 = index1 - 1;
                    for (var index2 = 1; index2 <= num6; ++index2)
                        v[index2] = qp[index2 + num7, index1];
                    v[1] = 1.0;
                    if (fromtheright)
                        Reflections.ApplyReflectionFromTheRight(ref z, tauq[index1], ref v, 1, zrows, index1, m, ref work);
                    else
                        Reflections.ApplyReflectionFromTheLeft(ref z, tauq[index1], ref v, index1, m, 1, zcolumns,
                            ref work);
                    index1 += num4;
                } while (index1 != num3 + num4);
            }
            else
            {
                int num8;
                int num9;
                int num10;
                if (fromtheright)
                {
                    num8 = 1;
                    num9 = m - 1;
                    num10 = 1;
                }
                else
                {
                    num8 = m - 1;
                    num9 = 1;
                    num10 = -1;
                }

                if (dotranspose)
                {
                    var num11 = num8;
                    num8 = num9;
                    num9 = num11;
                    num10 = -num10;
                }

                if (m - 1 <= 0)
                    return;
                var index3 = num8;
                do
                {
                    var num12 = m - index3;
                    var num13 = index3 + 1 - 1;
                    for (var index4 = 1; index4 <= num12; ++index4)
                        v[index4] = qp[index4 + num13, index3];
                    v[1] = 1.0;
                    if (fromtheright)
                        Reflections.ApplyReflectionFromTheRight(ref z, tauq[index3], ref v, 1, zrows, index3 + 1, m,
                            ref work);
                    else
                        Reflections.ApplyReflectionFromTheLeft(ref z, tauq[index3], ref v, index3 + 1, m, 1, zcolumns,
                            ref work);
                    index3 += num10;
                } while (index3 != num9 + num10);
            }
        }

        public static void unpackptfrombidiagonal(
            ref double[,] qp,
            int m,
            int n,
            ref double[] taup,
            int ptrows,
            ref double[,] pt)
        {
            var numArray1 = new double[0];
            var numArray2 = new double[0];
            if ((m == 0) | (n == 0) | (ptrows == 0))
                return;
            pt = new double[ptrows + 1, n + 1];
            var v = new double[n + 1];
            var work = new double[ptrows + 1];
            for (var index1 = 1; index1 <= ptrows; ++index1)
            for (var index2 = 1; index2 <= n; ++index2)
                pt[index1, index2] = index1 != index2 ? 0.0 : 1.0;
            if (m >= n)
                for (var index3 = Math.Min(n - 1, ptrows - 1); index3 >= 1; --index3)
                {
                    var num1 = n - index3;
                    var num2 = index3 + 1 - 1;
                    for (var index4 = 1; index4 <= num1; ++index4)
                        v[index4] = qp[index3, index4 + num2];
                    v[1] = 1.0;
                    Reflections.ApplyReflectionFromTheRight(ref pt, taup[index3], ref v, 1, ptrows, index3 + 1, n,
                        ref work);
                }
            else
                for (var n1 = Math.Min(m, ptrows); n1 >= 1; --n1)
                {
                    var num3 = n - n1 + 1;
                    var num4 = n1 - 1;
                    for (var index = 1; index <= num3; ++index)
                        v[index] = qp[n1, index + num4];
                    v[1] = 1.0;
                    Reflections.ApplyReflectionFromTheRight(ref pt, taup[n1], ref v, 1, ptrows, n1, n, ref work);
                }
        }

        public static void multiplybypfrombidiagonal(
            ref double[,] qp,
            int m,
            int n,
            ref double[] taup,
            ref double[,] z,
            int zrows,
            int zcolumns,
            bool fromtheright,
            bool dotranspose)
        {
            var numArray1 = new double[0];
            var numArray2 = new double[0];
            if ((m <= 0) | (n <= 0) | (zrows <= 0) | (zcolumns <= 0))
                return;
            var num1 = Math.Max(Math.Max(Math.Max(m, n), zrows), zcolumns);
            numArray1 = new double[num1 + 1];
            numArray2 = new double[num1 + 1];
            var v = new double[num1 + 1];
            var work = new double[num1 + 1];
            if (m >= n)
            {
                int num2;
                int num3;
                int num4;
                if (fromtheright)
                {
                    num2 = n - 1;
                    num3 = 1;
                    num4 = -1;
                }
                else
                {
                    num2 = 1;
                    num3 = n - 1;
                    num4 = 1;
                }

                if (!dotranspose)
                {
                    var num5 = num2;
                    num2 = num3;
                    num3 = num5;
                    num4 = -num4;
                }

                if (n - 1 <= 0)
                    return;
                var index1 = num2;
                do
                {
                    var num6 = n - index1;
                    var num7 = index1 + 1 - 1;
                    for (var index2 = 1; index2 <= num6; ++index2)
                        v[index2] = qp[index1, index2 + num7];
                    v[1] = 1.0;
                    if (fromtheright)
                        Reflections.ApplyReflectionFromTheRight(ref z, taup[index1], ref v, 1, zrows, index1 + 1, n,
                            ref work);
                    else
                        Reflections.ApplyReflectionFromTheLeft(ref z, taup[index1], ref v, index1 + 1, n, 1, zcolumns,
                            ref work);
                    index1 += num4;
                } while (index1 != num3 + num4);
            }
            else
            {
                int num8;
                int num9;
                int num10;
                if (fromtheright)
                {
                    num8 = m;
                    num9 = 1;
                    num10 = -1;
                }
                else
                {
                    num8 = 1;
                    num9 = m;
                    num10 = 1;
                }

                if (!dotranspose)
                {
                    var num11 = num8;
                    num8 = num9;
                    num9 = num11;
                    num10 = -num10;
                }

                var index3 = num8;
                do
                {
                    var num12 = n - index3 + 1;
                    var num13 = index3 - 1;
                    for (var index4 = 1; index4 <= num12; ++index4)
                        v[index4] = qp[index3, index4 + num13];
                    v[1] = 1.0;
                    if (fromtheright)
                        Reflections.ApplyReflectionFromTheRight(ref z, taup[index3], ref v, 1, zrows, index3, n, ref work);
                    else
                        Reflections.ApplyReflectionFromTheLeft(ref z, taup[index3], ref v, index3, n, 1, zcolumns,
                            ref work);
                    index3 += num10;
                } while (index3 != num9 + num10);
            }
        }

        public static void unpackdiagonalsfrombidiagonal(
            ref double[,] b,
            int m,
            int n,
            ref bool isupper,
            ref double[] d,
            ref double[] e)
        {
            isupper = m >= n;
            if ((m == 0) | (n == 0))
                return;
            if (isupper)
            {
                d = new double[n + 1];
                e = new double[n + 1];
                for (var index = 1; index <= n - 1; ++index)
                {
                    d[index] = b[index, index];
                    e[index] = b[index, index + 1];
                }

                d[n] = b[n, n];
            }
            else
            {
                d = new double[m + 1];
                e = new double[m + 1];
                for (var index = 1; index <= m - 1; ++index)
                {
                    d[index] = b[index, index];
                    e[index] = b[index + 1, index];
                }

                d[m] = b[m, m];
            }
        }
    }
}