using System;

public static class Reflections
{
    public static void GenerateReflection(ref double[] x, int n, ref double tau)
    {
        if (n <= 1)
        {
            tau = 0.0;
        }
        else
        {
            var num1 = x[1];
            var val2 = 0.0;
            for (var index = 2; index <= n; ++index)
                val2 = Math.Max(Math.Abs(x[index]), val2);
            var d = 0.0;
            if (val2 != 0.0)
            {
                for (var index = 2; index <= n; ++index)
                    d += AP.Math.Sqr(x[index] / val2);
                d = Math.Sqrt(d) * val2;
            }

            if (d == 0.0)
            {
                tau = 0.0;
            }
            else
            {
                var num2 = Math.Max(Math.Abs(num1), Math.Abs(d));
                var num3 = -(num2 * Math.Sqrt(AP.Math.Sqr(num1 / num2) + AP.Math.Sqr(d / num2)));
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

    public static void ApplyReflectionFromTheLeft(ref double[,] c, double tau, ref double[] v, int m1, int m2, int n1, int n2, ref double[] work) {
        if (tau == 0.0 || n1 > n2 || m1 > m2)
            return;
        for (var index = n1; index <= n2; ++index)
            work[index] = 0.0;
        for (var index1 = m1; index1 <= m2; ++index1)
        {
            var num = v[index1 + 1 - m1];
            for (var index2 = n1; index2 <= n2; ++index2)
                work[index2] += num * c[index1, index2];
        }

        for (var index3 = m1; index3 <= m2; ++index3)
        {
            var num = v[index3 - m1 + 1] * tau;
            for (var index4 = n1; index4 <= n2; ++index4)
                c[index3, index4] -= num * work[index4];
        }
    }

    public static void ApplyReflectionFromTheRight(ref double[,] c, double tau, ref double[] v, int m1, int m2, int n1, int n2, ref double[] work) {
        if (tau == 0.0 || n1 > n2 || m1 > m2)
            return;
        for (var index1 = m1; index1 <= m2; ++index1)
        {
            var num1 = 1 - n1;
            var num2 = 0.0;
            for (var index2 = n1; index2 <= n2; ++index2)
                num2 += c[index1, index2] * v[index2 + num1];
            work[index1] = num2;
        }

        for (var index3 = m1; index3 <= m2; ++index3)
        {
            var num3 = work[index3] * tau;
            var num4 = 1 - n1;
            for (var index4 = n1; index4 <= n2; ++index4)
                c[index3, index4] -= num3 * v[index4 + num4];
        }
    }

    private static void TestReflections()
    {
        var tau = 0.0;
        var num2 = 1000;
        var val1_1 = 0.0;
        var val1_2 = 0.0;
        var val1_3 = 0.0;
        for (var index1 = 1; index1 <= num2; ++index1)
        {
            var num3 = 1 + AP.Math.RandomInteger(10);
            var num4 = 1 + AP.Math.RandomInteger(10);
            var num5 = Math.Max(num4, num3);
            var numArray8 = new double[num5 + 1];
            var numArray9 = new double[num5 + 1];
            var work = new double[num5 + 1];
            var numArray10 = new double[num5 + 1, num5 + 1];
            var numArray11 = new double[num5 + 1, num5 + 1];
            var c = new double[num5 + 1, num5 + 1];
            var numArray12 = new double[num5 + 1, num5 + 1];
            for (var index2 = 1; index2 <= num3; ++index2)
            {
                numArray8[index2] = 2.0 * AP.Math.RandomReal() - 1.0;
                numArray9[index2] = numArray8[index2];
            }

            GenerateReflection(ref numArray9, num3, ref tau);
            var num6 = numArray9[1];
            numArray9[1] = 1.0;
            for (var index3 = 1; index3 <= num3; ++index3)
                for (var index4 = 1; index4 <= num3; ++index4)
                    numArray10[index3, index4] = index3 != index4
                        ? -(tau * numArray9[index3] * numArray9[index4])
                        : 1.0 - tau * numArray9[index3] * numArray9[index4];
            var num7 = 0.0;
            for (var index5 = 1; index5 <= num3; ++index5)
            {
                var num8 = 0.0;
                for (var index6 = 1; index6 <= num3; ++index6)
                    num8 += numArray10[index5, index6] * numArray8[index6];
                num7 = index5 != 1 ? Math.Max(num7, Math.Abs(num8)) : Math.Max(num7, Math.Abs(num8 - num6));
            }

            val1_3 = Math.Max(val1_3, num7);
            for (var index7 = 1; index7 <= num4; ++index7)
            {
                numArray8[index7] = 2.0 * AP.Math.RandomReal() - 1.0;
                numArray9[index7] = numArray8[index7];
            }

            for (var index8 = 1; index8 <= num4; ++index8)
            {
                for (var index9 = 1; index9 <= num3; ++index9)
                {
                    numArray11[index8, index9] = 2.0 * AP.Math.RandomReal() - 1.0;
                    c[index8, index9] = numArray11[index8, index9];
                }
            }

            GenerateReflection(ref numArray9, num4, ref tau);
            numArray9[1] = 1.0;
            ApplyReflectionFromTheLeft(ref c, tau, ref numArray9, 1, num4, 1, num3, ref work);
            for (var index10 = 1; index10 <= num4; ++index10)
                for (var index11 = 1; index11 <= num4; ++index11)
                    numArray10[index10, index11] = index10 != index11
                        ? -(tau * numArray9[index10] * numArray9[index11])
                        : 1.0 - tau * numArray9[index10] * numArray9[index11];
            for (var index12 = 1; index12 <= num4; ++index12)
            {
                for (var index13 = 1; index13 <= num3; ++index13)
                {
                    var num9 = 0.0;
                    for (var index14 = 1; index14 <= num4; ++index14)
                        num9 += numArray10[index12, index14] * numArray11[index14, index13];
                    numArray12[index12, index13] = num9;
                }
            }

            var num10 = 0.0;
            for (var index15 = 1; index15 <= num4; ++index15)
                for (var index16 = 1; index16 <= num3; ++index16)
                    num10 = Math.Max(num10, Math.Abs(c[index15, index16] - numArray12[index15, index16]));
            val1_2 = Math.Max(val1_2, num10);
            for (var index17 = 1; index17 <= num3; ++index17)
            {
                numArray8[index17] = 2.0 * AP.Math.RandomReal() - 1.0;
                numArray9[index17] = numArray8[index17];
            }

            for (var index18 = 1; index18 <= num4; ++index18)
            {
                for (var index19 = 1; index19 <= num3; ++index19)
                {
                    numArray11[index18, index19] = 2.0 * AP.Math.RandomReal() - 1.0;
                    c[index18, index19] = numArray11[index18, index19];
                }
            }

            GenerateReflection(ref numArray9, num3, ref tau);
            numArray9[1] = 1.0;
            ApplyReflectionFromTheRight(ref c, tau, ref numArray9, 1, num4, 1, num3, ref work);
            for (var index20 = 1; index20 <= num3; ++index20)
                for (var index21 = 1; index21 <= num3; ++index21)
                    numArray10[index20, index21] = index20 != index21
                        ? -(tau * numArray9[index20] * numArray9[index21])
                        : 1.0 - tau * numArray9[index20] * numArray9[index21];
            for (var index22 = 1; index22 <= num4; ++index22)
            {
                for (var index23 = 1; index23 <= num3; ++index23)
                {
                    var num11 = 0.0;
                    for (var index24 = 1; index24 <= num3; ++index24)
                        num11 += numArray11[index22, index24] * numArray10[index24, index23];
                    numArray12[index22, index23] = num11;
                }
            }

            var num12 = 0.0;
            for (var index25 = 1; index25 <= num4; ++index25)
            for (var index26 = 1; index26 <= num3; ++index26)
                num12 = Math.Max(num12, Math.Abs(c[index25, index26] - numArray12[index25, index26]));
            val1_1 = Math.Max(val1_1, num12);
        }

        var x = new double[11];
        for (var index = 1; index <= 10; ++index)
            x[index] = 1.0000000000000001E+298 * (2.0 * AP.Math.RandomReal() - 1.0);
        GenerateReflection(ref x, 10, ref tau);
        Console.Write("TESTING REFLECTIONS");
        Console.WriteLine();
        Console.Write("Pass count is ");
        Console.Write("{0,0:d}", num2);
        Console.WriteLine();
        Console.Write("Generate     absolute error is       ");
        Console.Write("{0,5:E3}", val1_3);
        Console.WriteLine();
        Console.Write("Apply(Left)  absolute error is       ");
        Console.Write("{0,5:E3}", val1_2);
        Console.WriteLine();
        Console.Write("Apply(Right) absolute error is       ");
        Console.Write("{0,5:E3}", val1_1);
        Console.WriteLine();
        Console.Write("Overflow crash test passed");
        Console.WriteLine();
    }
}