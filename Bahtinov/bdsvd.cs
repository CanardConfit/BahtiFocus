using System;

namespace CanardConfit.NINA.BahtiFocus.Bahtinov {
    public static class bdsvd {
        public static bool rmatrixbdsvd(
            ref double[] d,
            double[] e,
            int n,
            bool isupper,
            bool isfractionalaccuracyrequired,
            ref double[,] u,
            int nru,
            ref double[,] c,
            int ncc,
            ref double[,] vt,
            int ncvt) {
            var numArray = new double[0];
            var e1 = new double[0];
            e = (double[])e.Clone();
            var d1 = new double[n + 1];
            var num1 = -1;
            for (var index = 1; index <= n; ++index)
                d1[index] = d[index + num1];
            if (n > 1) {
                e1 = new double[n - 1 + 1];
                var num2 = -1;
                for (var index = 1; index <= n - 1; ++index)
                    e1[index] = e[index + num2];
            }

            var flag = bidiagonalsvddecompositioninternal(ref d1, e1, n, isupper, isfractionalaccuracyrequired, ref u,
                0,
                nru, ref c, 0, ncc, ref vt, 0, ncvt);
            var num3 = 1;
            for (var index = 0; index <= n - 1; ++index)
                d[index] = d1[index + num3];
            return flag;
        }

        public static bool bidiagonalsvddecomposition(
            ref double[] d,
            double[] e,
            int n,
            bool isupper,
            bool isfractionalaccuracyrequired,
            ref double[,] u,
            int nru,
            ref double[,] c,
            int ncc,
            ref double[,] vt,
            int ncvt) {
            e = (double[])e.Clone();
            return bidiagonalsvddecompositioninternal(ref d, e, n, isupper, isfractionalaccuracyrequired, ref u, 1, nru,
                ref c, 1, ncc, ref vt, 1, ncvt);
        }

        private static bool bidiagonalsvddecompositioninternal(
            ref double[] d,
            double[] e,
            int n,
            bool isupper,
            bool isfractionalaccuracyrequired,
            ref double[,] u,
            int ustart,
            int nru,
            ref double[,] c,
            int cstart,
            int ncc,
            ref double[,] vt,
            int vstart,
            int ncvt) {
            var index1 = 0;
            var num1 = 0.0;
            var num2 = 0.0;
            var cs1 = 0.0;
            var sn1 = 0.0;
            var num3 = 0.0;
            var ssmin1 = 0.0;
            var ssmin2 = 0.0;
            var ssmax = 0.0;
            var num4 = 0.0;
            var num5 = 0.0;
            var sn2 = 0.0;
            var numArray1 = new double[0];
            var numArray2 = new double[0];
            var numArray3 = new double[0];
            var numArray4 = new double[0];
            var numArray5 = new double[0];
            var numArray6 = new double[0];
            var numArray7 = new double[0];
            var numArray8 = new double[0];
            var r = 0.0;
            e = (double[])e.Clone();
            var flag1 = true;
            switch (n) {
                case 0:
                    return flag1;
                case 1:
                    if (d[1] < 0.0) {
                        d[1] = -d[1];
                        if (ncvt > 0)
                            for (var index2 = vstart; index2 <= vstart + ncvt - 1; ++index2)
                                vt[vstart, index2] = -1.0 * vt[vstart, index2];
                    }

                    return flag1;
                default:
                    var c1 = new double[n - 1 + 1];
                    var s1 = new double[n - 1 + 1];
                    var c2 = new double[n - 1 + 1];
                    var s2 = new double[n - 1 + 1];
                    var m2 = ustart + Math.Max(nru - 1, 0);
                    var n2_1 = vstart + Math.Max(ncvt - 1, 0);
                    var n2_2 = cstart + Math.Max(ncc - 1, 0);
                    var work1 = new double[m2 + 1];
                    var work2 = new double[n2_1 + 1];
                    var work3 = new double[n2_2 + 1];
                    var num6 = 12;
                    var isforward = true;
                    var numArray9 = new double[n + 1];
                    for (var index3 = 1; index3 <= n - 1; ++index3)
                        numArray9[index3] = e[index3];
                    e = new double[n + 1];
                    for (var index4 = 1; index4 <= n - 1; ++index4)
                        e[index4] = numArray9[index4];
                    e[n] = 0.0;
                    var num7 = 0;
                    var num8 = 5E-16;
                    var num9 = 1E-300;
                    if (!isupper) {
                        for (var index5 = 1; index5 <= n - 1; ++index5) {
                            Rotations.GenerateRotation(d[index5], e[index5], ref cs1, ref sn2, ref num3);
                            d[index5] = num3;
                            e[index5] = sn2 * d[index5 + 1];
                            d[index5 + 1] = cs1 * d[index5 + 1];
                            c1[index5] = cs1;
                            s1[index5] = sn2;
                        }

                        if (nru > 0)
                            Rotations.ApplyRotationsFromTheRight(isforward, ustart, m2, 1 + ustart - 1, n + ustart - 1,
                                ref c1, ref s1, ref u, ref work1);
                        if (ncc > 0)
                            Rotations.ApplyRotationsFromTheLeft(isforward, 1 + cstart - 1, n + cstart - 1, cstart, n2_2,
                                ref c1, ref s1, ref c, ref work3);
                    }

                    var num10 = Math.Max(10.0, Math.Min(100.0, Math.Pow(num8, -0.125))) * num8;
                    if (!isfractionalaccuracyrequired)
                        num10 = -num10;
                    var val1_1 = 0.0;
                    for (var index6 = 1; index6 <= n; ++index6)
                        val1_1 = Math.Max(val1_1, Math.Abs(d[index6]));
                    for (var index7 = 1; index7 <= n - 1; ++index7)
                        val1_1 = Math.Max(val1_1, Math.Abs(e[index7]));
                    var val1_2 = 0.0;
                    double num11;
                    if (num10 >= 0.0) {
                        var val1_3 = Math.Abs(d[1]);
                        if (val1_3 != 0.0) {
                            var val2 = val1_3;
                            for (var index8 = 2; index8 <= n; ++index8) {
                                val2 = Math.Abs(d[index8]) * (val2 / (val2 + Math.Abs(e[index8 - 1])));
                                val1_3 = Math.Min(val1_3, val2);
                                if (val1_3 == 0.0)
                                    break;
                            }
                        }

                        var num12 = val1_3 / Math.Sqrt(n);
                        num11 = Math.Max(num10 * num12, num6 * n * n * num9);
                    } else {
                        num11 = Math.Max(Math.Abs(num10) * val1_1, num6 * n * n * num9);
                    }

                    var num13 = num6 * n * n;
                    var num14 = 0;
                    var num15 = -1;
                    var num16 = -1;
                    var index9 = n;
                    while (index9 > 1) {
                        if (num14 > num13)
                            return false;
                        if ((num10 < 0.0) & (Math.Abs(d[index9]) <= num11))
                            d[index9] = 0.0;
                        var val1_4 = Math.Abs(d[index9]);
                        var val1_5 = val1_4;
                        var flag2 = false;
                        for (var index10 = 1; index10 <= index9 - 1; ++index10) {
                            index1 = index9 - index10;
                            var num17 = Math.Abs(d[index1]);
                            var val2 = Math.Abs(e[index1]);
                            if ((num10 < 0.0) & (num17 <= num11))
                                d[index1] = 0.0;
                            if (val2 <= num11) {
                                flag2 = true;
                                break;
                            }

                            val1_5 = Math.Min(val1_5, num17);
                            val1_4 = Math.Max(val1_4, Math.Max(num17, val2));
                        }

                        if (!flag2) {
                            index1 = 0;
                        } else {
                            e[index1] = 0.0;
                            if (index1 == index9 - 1) {
                                --index9;
                                continue;
                            }
                        }

                        ++index1;
                        if (index1 == index9 - 1) {
                            svdv2x2(d[index9 - 1], e[index9 - 1], d[index9], ref ssmin2, ref ssmax, ref num5, ref num2,
                                ref num4, ref num1);
                            d[index9 - 1] = ssmax;
                            e[index9 - 1] = 0.0;
                            d[index9] = ssmin2;
                            if (ncvt > 0) {
                                var index11 = index9 + (vstart - 1);
                                var index12 = index9 - 1 + (vstart - 1);
                                for (var index13 = vstart; index13 <= n2_1; ++index13)
                                    work2[index13] = num2 * vt[index12, index13];
                                for (var index14 = vstart; index14 <= n2_1; ++index14)
                                    work2[index14] = work2[index14] + num5 * vt[index11, index14];
                                for (var index15 = vstart; index15 <= n2_1; ++index15)
                                    vt[index11, index15] = num2 * vt[index11, index15];
                                for (var index16 = vstart; index16 <= n2_1; ++index16)
                                    vt[index11, index16] = vt[index11, index16] - num5 * vt[index12, index16];
                                for (var index17 = vstart; index17 <= n2_1; ++index17)
                                    vt[index12, index17] = work2[index17];
                            }

                            if (nru > 0) {
                                var index18 = index9 + ustart - 1;
                                var index19 = index9 - 1 + ustart - 1;
                                for (var index20 = ustart; index20 <= m2; ++index20)
                                    work1[index20] = num1 * u[index20, index19];
                                for (var index21 = ustart; index21 <= m2; ++index21)
                                    work1[index21] = work1[index21] + num4 * u[index21, index18];
                                for (var index22 = ustart; index22 <= m2; ++index22)
                                    u[index22, index18] = num1 * u[index22, index18];
                                for (var index23 = ustart; index23 <= m2; ++index23)
                                    u[index23, index18] = u[index23, index18] - num4 * u[index23, index19];
                                for (var index24 = ustart; index24 <= m2; ++index24)
                                    u[index24, index19] = work1[index24];
                            }

                            if (ncc > 0) {
                                var index25 = index9 + cstart - 1;
                                var index26 = index9 - 1 + cstart - 1;
                                for (var index27 = cstart; index27 <= n2_2; ++index27)
                                    work3[index27] = num1 * c[index26, index27];
                                for (var index28 = cstart; index28 <= n2_2; ++index28)
                                    work3[index28] = work3[index28] + num4 * c[index25, index28];
                                for (var index29 = cstart; index29 <= n2_2; ++index29)
                                    c[index25, index29] = num1 * c[index25, index29];
                                for (var index30 = cstart; index30 <= n2_2; ++index30)
                                    c[index25, index30] = c[index25, index30] - num4 * c[index26, index30];
                                for (var index31 = cstart; index31 <= n2_2; ++index31)
                                    c[index26, index31] = work3[index31];
                            }

                            index9 -= 2;
                        } else {
                            var flag3 = false;
                            if ((num7 == 1) & (Math.Abs(d[index1]) < 0.001 * Math.Abs(d[index9])))
                                flag3 = true;
                            if ((num7 == 2) & (Math.Abs(d[index9]) < 0.001 * Math.Abs(d[index1])))
                                flag3 = true;
                            if ((index1 != num15) | (index9 != num16) | flag3)
                                num7 = Math.Abs(d[index1]) < Math.Abs(d[index9]) ? 2 : 1;
                            if (num7 == 1) {
                                if ((Math.Abs(e[index9 - 1]) <= Math.Abs(num10) * Math.Abs(d[index9])) |
                                    ((num10 < 0.0) & (Math.Abs(e[index9 - 1]) <= num11))) {
                                    e[index9 - 1] = 0.0;
                                    continue;
                                }

                                if (num10 >= 0.0) {
                                    var val2 = Math.Abs(d[index1]);
                                    val1_2 = val2;
                                    var flag4 = false;
                                    for (var index32 = index1; index32 <= index9 - 1; ++index32) {
                                        if (Math.Abs(e[index32]) <= num10 * val2) {
                                            e[index32] = 0.0;
                                            flag4 = true;
                                            break;
                                        }

                                        val2 = Math.Abs(d[index32 + 1]) * (val2 / (val2 + Math.Abs(e[index32])));
                                        val1_2 = Math.Min(val1_2, val2);
                                    }

                                    if (flag4)
                                        continue;
                                }
                            } else {
                                if ((Math.Abs(e[index1]) <= Math.Abs(num10) * Math.Abs(d[index1])) |
                                    ((num10 < 0.0) & (Math.Abs(e[index1]) <= num11))) {
                                    e[index1] = 0.0;
                                    continue;
                                }

                                if (num10 >= 0.0) {
                                    var val2 = Math.Abs(d[index9]);
                                    val1_2 = val2;
                                    var flag5 = false;
                                    for (var index33 = index9 - 1; index33 >= index1; --index33) {
                                        if (Math.Abs(e[index33]) <= num10 * val2) {
                                            e[index33] = 0.0;
                                            flag5 = true;
                                            break;
                                        }

                                        val2 = Math.Abs(d[index33]) * (val2 / (val2 + Math.Abs(e[index33])));
                                        val1_2 = Math.Min(val1_2, val2);
                                    }

                                    if (flag5)
                                        continue;
                                }
                            }

                            num15 = index1;
                            num16 = index9;
                            if ((num10 >= 0.0) & (n * num10 * (val1_2 / val1_4) <= Math.Max(num8, 0.01 * num10))) {
                                ssmin1 = 0.0;
                            } else {
                                double num18;
                                if (num7 == 1) {
                                    num18 = Math.Abs(d[index1]);
                                    svd2x2(d[index9 - 1], e[index9 - 1], d[index9], ref ssmin1, ref num3);
                                } else {
                                    num18 = Math.Abs(d[index9]);
                                    svd2x2(d[index1], e[index1], d[index1 + 1], ref ssmin1, ref num3);
                                }

                                if (num18 > 0.0 && CustomMath.Sqr(ssmin1 / num18) < num8)
                                    ssmin1 = 0.0;
                            }

                            num14 = num14 + index9 - index1;
                            double cs2;
                            if (ssmin1 == 0.0) {
                                if (num7 == 1) {
                                    cs1 = 1.0;
                                    cs2 = 1.0;
                                    for (var index34 = index1; index34 <= index9 - 1; ++index34) {
                                        Rotations.GenerateRotation(d[index34] * cs1, e[index34], ref cs1, ref sn2,
                                            ref num3);
                                        if (index34 > index1)
                                            e[index34 - 1] = sn1 * num3;
                                        Rotations.GenerateRotation(cs2 * num3, d[index34 + 1] * sn2, ref cs2, ref sn1,
                                            ref r);
                                        d[index34] = r;
                                        c1[index34 - index1 + 1] = cs1;
                                        s1[index34 - index1 + 1] = sn2;
                                        c2[index34 - index1 + 1] = cs2;
                                        s2[index34 - index1 + 1] = sn1;
                                    }

                                    var num19 = d[index9] * cs1;
                                    d[index9] = num19 * cs2;
                                    e[index9 - 1] = num19 * sn1;
                                    if (ncvt > 0)
                                        Rotations.ApplyRotationsFromTheLeft(isforward, index1 + vstart - 1,
                                            index9 + vstart - 1, vstart, n2_1, ref c1, ref s1, ref vt, ref work2);
                                    if (nru > 0)
                                        Rotations.ApplyRotationsFromTheRight(isforward, ustart, m2, index1 + ustart - 1,
                                            index9 + ustart - 1, ref c2, ref s2, ref u, ref work1);
                                    if (ncc > 0)
                                        Rotations.ApplyRotationsFromTheLeft(isforward, index1 + cstart - 1,
                                            index9 + cstart - 1, cstart, n2_2, ref c2, ref s2, ref c, ref work3);
                                    if (Math.Abs(e[index9 - 1]) <= num11)
                                        e[index9 - 1] = 0.0;
                                } else {
                                    cs1 = 1.0;
                                    cs2 = 1.0;
                                    for (var index35 = index9; index35 >= index1 + 1; --index35) {
                                        Rotations.GenerateRotation(d[index35] * cs1, e[index35 - 1], ref cs1, ref sn2,
                                            ref num3);
                                        if (index35 < index9)
                                            e[index35] = sn1 * num3;
                                        Rotations.GenerateRotation(cs2 * num3, d[index35 - 1] * sn2, ref cs2, ref sn1,
                                            ref r);
                                        d[index35] = r;
                                        c1[index35 - index1] = cs1;
                                        s1[index35 - index1] = -sn2;
                                        c2[index35 - index1] = cs2;
                                        s2[index35 - index1] = -sn1;
                                    }

                                    var num20 = d[index1] * cs1;
                                    d[index1] = num20 * cs2;
                                    e[index1] = num20 * sn1;
                                    if (ncvt > 0)
                                        Rotations.ApplyRotationsFromTheLeft(!isforward, index1 + vstart - 1,
                                            index9 + vstart - 1, vstart, n2_1, ref c2, ref s2, ref vt, ref work2);
                                    if (nru > 0)
                                        Rotations.ApplyRotationsFromTheRight(!isforward, ustart, m2,
                                            index1 + ustart - 1,
                                            index9 + ustart - 1, ref c1, ref s1, ref u, ref work1);
                                    if (ncc > 0)
                                        Rotations.ApplyRotationsFromTheLeft(!isforward, index1 + cstart - 1,
                                            index9 + cstart - 1, cstart, n2_2, ref c1, ref s1, ref c, ref work3);
                                    if (Math.Abs(e[index1]) <= num11)
                                        e[index1] = 0.0;
                                }
                            } else if (num7 == 1) {
                                var f1 = (Math.Abs(d[index1]) - ssmin1) *
                                         (extsignbdsqr(1.0, d[index1]) + ssmin1 / d[index1]);
                                var g = e[index1];
                                for (var index36 = index1; index36 <= index9 - 1; ++index36) {
                                    Rotations.GenerateRotation(f1, g, ref num2, ref num5, ref num3);
                                    if (index36 > index1)
                                        e[index36 - 1] = num3;
                                    var f2 = num2 * d[index36] + num5 * e[index36];
                                    e[index36] = num2 * e[index36] - num5 * d[index36];
                                    g = num5 * d[index36 + 1];
                                    d[index36 + 1] = num2 * d[index36 + 1];
                                    Rotations.GenerateRotation(f2, g, ref num1, ref num4, ref num3);
                                    d[index36] = num3;
                                    f1 = num1 * e[index36] + num4 * d[index36 + 1];
                                    d[index36 + 1] = num1 * d[index36 + 1] - num4 * e[index36];
                                    if (index36 < index9 - 1) {
                                        g = num4 * e[index36 + 1];
                                        e[index36 + 1] = num1 * e[index36 + 1];
                                    }

                                    c1[index36 - index1 + 1] = num2;
                                    s1[index36 - index1 + 1] = num5;
                                    c2[index36 - index1 + 1] = num1;
                                    s2[index36 - index1 + 1] = num4;
                                }

                                e[index9 - 1] = f1;
                                if (ncvt > 0)
                                    Rotations.ApplyRotationsFromTheLeft(isforward, index1 + vstart - 1,
                                        index9 + vstart - 1,
                                        vstart, n2_1, ref c1, ref s1, ref vt, ref work2);
                                if (nru > 0)
                                    Rotations.ApplyRotationsFromTheRight(isforward, ustart, m2, index1 + ustart - 1,
                                        index9 + ustart - 1, ref c2, ref s2, ref u, ref work1);
                                if (ncc > 0)
                                    Rotations.ApplyRotationsFromTheLeft(isforward, index1 + cstart - 1,
                                        index9 + cstart - 1,
                                        cstart, n2_2, ref c2, ref s2, ref c, ref work3);
                                if (Math.Abs(e[index9 - 1]) <= num11)
                                    e[index9 - 1] = 0.0;
                            } else {
                                var f3 = (Math.Abs(d[index9]) - ssmin1) *
                                         (extsignbdsqr(1.0, d[index9]) + ssmin1 / d[index9]);
                                var g = e[index9 - 1];
                                for (var index37 = index9; index37 >= index1 + 1; --index37) {
                                    Rotations.GenerateRotation(f3, g, ref num2, ref num5, ref num3);
                                    if (index37 < index9)
                                        e[index37] = num3;
                                    var f4 = num2 * d[index37] + num5 * e[index37 - 1];
                                    e[index37 - 1] = num2 * e[index37 - 1] - num5 * d[index37];
                                    g = num5 * d[index37 - 1];
                                    d[index37 - 1] = num2 * d[index37 - 1];
                                    Rotations.GenerateRotation(f4, g, ref num1, ref num4, ref num3);
                                    d[index37] = num3;
                                    f3 = num1 * e[index37 - 1] + num4 * d[index37 - 1];
                                    d[index37 - 1] = num1 * d[index37 - 1] - num4 * e[index37 - 1];
                                    if (index37 > index1 + 1) {
                                        g = num4 * e[index37 - 2];
                                        e[index37 - 2] = num1 * e[index37 - 2];
                                    }

                                    c1[index37 - index1] = num2;
                                    s1[index37 - index1] = -num5;
                                    c2[index37 - index1] = num1;
                                    s2[index37 - index1] = -num4;
                                }

                                e[index1] = f3;
                                if (Math.Abs(e[index1]) <= num11)
                                    e[index1] = 0.0;
                                if (ncvt > 0)
                                    Rotations.ApplyRotationsFromTheLeft(!isforward, index1 + vstart - 1,
                                        index9 + vstart - 1, vstart, n2_1, ref c2, ref s2, ref vt, ref work2);
                                if (nru > 0)
                                    Rotations.ApplyRotationsFromTheRight(!isforward, ustart, m2, index1 + ustart - 1,
                                        index9 + ustart - 1, ref c1, ref s1, ref u, ref work1);
                                if (ncc > 0)
                                    Rotations.ApplyRotationsFromTheLeft(!isforward, index1 + cstart - 1,
                                        index9 + cstart - 1, cstart, n2_2, ref c1, ref s1, ref c, ref work3);
                            }
                        }
                    }

                    for (var index38 = 1; index38 <= n; ++index38)
                        if (d[index38] < 0.0) {
                            d[index38] = -d[index38];
                            if (ncvt > 0)
                                for (var index39 = vstart; index39 <= n2_1; ++index39)
                                    vt[index38 + vstart - 1, index39] = -1.0 * vt[index38 + vstart - 1, index39];
                        }

                    for (var index40 = 1; index40 <= n - 1; ++index40) {
                        var index41 = 1;
                        var num21 = d[1];
                        for (var index42 = 2; index42 <= n + 1 - index40; ++index42)
                            if (d[index42] <= num21) {
                                index41 = index42;
                                num21 = d[index42];
                            }

                        if (index41 != n + 1 - index40) {
                            d[index41] = d[n + 1 - index40];
                            d[n + 1 - index40] = num21;
                            if (ncvt > 0) {
                                var num22 = n + 1 - index40;
                                for (var index43 = vstart; index43 <= n2_1; ++index43)
                                    work2[index43] = vt[index41 + vstart - 1, index43];
                                for (var index44 = vstart; index44 <= n2_1; ++index44)
                                    vt[index41 + vstart - 1, index44] = vt[num22 + vstart - 1, index44];
                                for (var index45 = vstart; index45 <= n2_1; ++index45)
                                    vt[num22 + vstart - 1, index45] = work2[index45];
                            }

                            if (nru > 0) {
                                var num23 = n + 1 - index40;
                                for (var index46 = ustart; index46 <= m2; ++index46)
                                    work1[index46] = u[index46, index41 + ustart - 1];
                                for (var index47 = ustart; index47 <= m2; ++index47)
                                    u[index47, index41 + ustart - 1] = u[index47, num23 + ustart - 1];
                                for (var index48 = ustart; index48 <= m2; ++index48)
                                    u[index48, num23 + ustart - 1] = work1[index48];
                            }

                            if (ncc > 0) {
                                var num24 = n + 1 - index40;
                                for (var index49 = cstart; index49 <= n2_2; ++index49)
                                    work3[index49] = c[index41 + cstart - 1, index49];
                                for (var index50 = cstart; index50 <= n2_2; ++index50)
                                    c[index41 + cstart - 1, index50] = c[num24 + cstart - 1, index50];
                                for (var index51 = cstart; index51 <= n2_2; ++index51)
                                    c[num24 + cstart - 1, index51] = work3[index51];
                            }
                        }
                    }

                    return flag1;
            }
        }

        private static double extsignbdsqr(double a, double b) {
            return b < 0.0 ? -Math.Abs(a) : Math.Abs(a);
        }

        private static void svd2x2(double f, double g, double h, ref double ssmin, ref double ssmax) {
            var val1_1 = Math.Abs(f);
            var val2_1 = Math.Abs(g);
            var val2_2 = Math.Abs(h);
            var num1 = Math.Min(val1_1, val2_2);
            var val1_2 = Math.Max(val1_1, val2_2);
            if (num1 == 0.0) {
                ssmin = 0.0;
                if (val1_2 == 0.0)
                    ssmax = val2_1;
                else
                    ssmax = Math.Max(val1_2, val2_1) *
                            Math.Sqrt(1.0 + CustomMath.Sqr(Math.Min(val1_2, val2_1) / Math.Max(val1_2, val2_1)));
            } else if (val2_1 < val1_2) {
                var num2 = 1.0 + num1 / val1_2;
                var num3 = (val1_2 - num1) / val1_2;
                var num4 = CustomMath.Sqr(val2_1 / val1_2);
                var num5 = 2.0 / (Math.Sqrt(num2 * num2 + num4) + Math.Sqrt(num3 * num3 + num4));
                ssmin = num1 * num5;
                ssmax = val1_2 / num5;
            } else {
                var num6 = val1_2 / val2_1;
                if (num6 == 0.0) {
                    ssmin = num1 * val1_2 / val2_1;
                    ssmax = val2_1;
                } else {
                    var num7 = 1.0 + num1 / val1_2;
                    var num8 = (val1_2 - num1) / val1_2;
                    var num9 = 1.0 / (Math.Sqrt(1.0 + CustomMath.Sqr(num7 * num6)) +
                                      Math.Sqrt(1.0 + CustomMath.Sqr(num8 * num6)));
                    ssmin = num1 * num9 * num6;
                    ssmin += ssmin;
                    ssmax = val2_1 / (num9 + num9);
                }
            }
        }

        private static void svdv2x2(
            double f,
            double g,
            double h,
            ref double ssmin,
            ref double ssmax,
            ref double snr,
            ref double csr,
            ref double snl,
            ref double csl) {
            var num1 = 0.0;
            var num2 = 0.0;
            var num3 = 0.0;
            var num4 = 0.0;
            var b1 = 0.0;
            var b2 = f;
            var num5 = Math.Abs(b2);
            var num6 = h;
            var num7 = Math.Abs(h);
            var num8 = 1;
            var flag1 = num7 > num5;
            if (flag1) {
                num8 = 3;
                var num9 = b2;
                b2 = num6;
                num6 = num9;
                var num10 = num5;
                num5 = num7;
                num7 = num10;
            }

            var b3 = g;
            var num11 = Math.Abs(b3);
            if (num11 == 0.0) {
                ssmin = num7;
                ssmax = num5;
                num1 = 1.0;
                num2 = 1.0;
                num3 = 0.0;
                num4 = 0.0;
            } else {
                var flag2 = true;
                if (num11 > num5) {
                    num8 = 2;
                    if (num5 / num11 < 5E-16) {
                        flag2 = false;
                        ssmax = num11;
                        if (num7 > 1.0) {
                            var num12 = num11 / num7;
                            ssmin = num5 / num12;
                        } else {
                            var num13 = num5 / num11;
                            ssmin = num13 * num7;
                        }

                        num1 = 1.0;
                        num3 = num6 / b3;
                        num4 = 1.0;
                        num2 = b2 / b3;
                    }
                }

                if (flag2) {
                    var a = num5 - num7;
                    var num14 = a != num5 ? a / num5 : 1.0;
                    var num15 = b3 / b2;
                    var num16 = 2.0 - num14;
                    var num17 = num15 * num15;
                    var num18 = Math.Sqrt(num16 * num16 + num17);
                    var num19 = num14 != 0.0 ? Math.Sqrt(num14 * num14 + num17) : Math.Abs(num15);
                    var num20 = 0.5 * (num18 + num19);
                    ssmin = num7 / num20;
                    ssmax = num5 * num20;
                    var num21 = num17 != 0.0 ? (num15 / (num18 + num16) + num15 / (num19 + num14)) * (1.0 + num20) :
                        num14 != 0.0 ? b3 / extsignbdsqr(a, b2) + num15 / num16 :
                        extsignbdsqr(2.0, b2) * extsignbdsqr(1.0, b3);
                    var num22 = Math.Sqrt(num21 * num21 + 4.0);
                    num2 = 2.0 / num22;
                    num4 = num21 / num22;
                    num1 = (num2 + num4 * num15) / num20;
                    num3 = num6 / b2 * num4 / num20;
                }
            }

            if (flag1) {
                csl = num4;
                snl = num2;
                csr = num3;
                snr = num1;
            } else {
                csl = num1;
                snl = num3;
                csr = num2;
                snr = num4;
            }

            if (num8 == 1)
                b1 = extsignbdsqr(1.0, csr) * extsignbdsqr(1.0, csl) * extsignbdsqr(1.0, f);
            if (num8 == 2)
                b1 = extsignbdsqr(1.0, snr) * extsignbdsqr(1.0, csl) * extsignbdsqr(1.0, g);
            if (num8 == 3)
                b1 = extsignbdsqr(1.0, snr) * extsignbdsqr(1.0, snl) * extsignbdsqr(1.0, h);
            ssmax = extsignbdsqr(ssmax, b1);
            ssmin = extsignbdsqr(ssmin, b1 * extsignbdsqr(1.0, f) * extsignbdsqr(1.0, h));
        }
    }
}