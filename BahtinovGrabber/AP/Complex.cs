// Decompiled with JetBrains decompiler
// Type: AP.Complex
// Assembly: Bahtinov_grabber, Version=1.0.0.1, Culture=neutral, PublicKeyToken=null
// MVID: 76A55480-8958-463E-9EFB-D618BEF833F5
// Assembly location: D:\Users\Tom\Downloads\Bahtinov_grabber64.exe

namespace AP
{
    public struct Complex
    {
        public double x;
        public double y;

        public Complex(double _x)
        {
            x = _x;
            y = 0.0;
        }

        public Complex(double _x, double _y)
        {
            x = _x;
            y = _y;
        }

        public static implicit operator Complex(double _x)
        {
            return new Complex(_x);
        }

        public static bool operator ==(Complex lhs, Complex rhs)
        {
            return (lhs.x == rhs.x) & (lhs.y == rhs.y);
        }

        public static bool operator !=(Complex lhs, Complex rhs)
        {
            return (lhs.x != rhs.x) | (lhs.y != rhs.y);
        }

        public static Complex operator +(Complex lhs)
        {
            return lhs;
        }

        public static Complex operator -(Complex lhs)
        {
            return new Complex(-lhs.x, -lhs.y);
        }

        public static Complex operator +(Complex lhs, Complex rhs)
        {
            return new Complex(lhs.x + rhs.x, lhs.y + rhs.y);
        }

        public static Complex operator -(Complex lhs, Complex rhs)
        {
            return new Complex(lhs.x - rhs.x, lhs.y - rhs.y);
        }

        public static Complex operator *(Complex lhs, Complex rhs)
        {
            return new Complex(lhs.x * rhs.x - lhs.y * rhs.y, lhs.x * rhs.y + lhs.y * rhs.x);
        }

        public static Complex operator /(Complex lhs, Complex rhs)
        {
            Complex complex;
            if (System.Math.Abs(rhs.y) < System.Math.Abs(rhs.x))
            {
                var num1 = rhs.y / rhs.x;
                var num2 = rhs.x + rhs.y * num1;
                complex.x = (lhs.x + lhs.y * num1) / num2;
                complex.y = (lhs.y - lhs.x * num1) / num2;
            }
            else
            {
                var num3 = rhs.x / rhs.y;
                var num4 = rhs.y + rhs.x * num3;
                complex.x = (lhs.y + lhs.x * num3) / num4;
                complex.y = (-lhs.x + lhs.y * num3) / num4;
            }

            return complex;
        }
    }
}