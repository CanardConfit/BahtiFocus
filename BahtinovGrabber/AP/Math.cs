// Decompiled with JetBrains decompiler
// Type: AP.Math
// Assembly: Bahtinov_grabber, Version=1.0.0.1, Culture=neutral, PublicKeyToken=null
// MVID: 76A55480-8958-463E-9EFB-D618BEF833F5
// Assembly location: D:\Users\Tom\Downloads\Bahtinov_grabber64.exe

using System;

namespace AP
{
    public class Math
    {
        public const double MachineEpsilon = 5E-16;
        public const double MaxRealNumber = 1E+300;
        public const double MinRealNumber = 1E-300;
        public static Random RndObject = new Random(DateTime.Now.Millisecond);

        public static double RandomReal()
        {
            lock (RndObject)
            {
                return RndObject.NextDouble();
            }
        }

        public static int RandomInteger(int N)
        {
            lock (RndObject)
            {
                return RndObject.Next(N);
            }
        }

        public static double Sqr(double X)
        {
            return X * X;
        }

        public static double AbsComplex(Complex z)
        {
            var num1 = System.Math.Abs(z.x);
            var num2 = System.Math.Abs(z.y);
            var num3 = num1 > num2 ? num1 : num2;
            var num4 = num1 < num2 ? num1 : num2;
            if (num4 == 0.0)
                return num3;
            var num5 = num4 / num3;
            return num3 * System.Math.Sqrt(1.0 + num5 * num5);
        }

        public static Complex Conj(Complex z)
        {
            return new Complex(z.x, -z.y);
        }

        public static Complex CSqr(Complex z)
        {
            return new Complex(z.x * z.x - z.y * z.y, 2.0 * z.x * z.y);
        }
    }
}