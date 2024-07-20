using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

using System.IO;

namespace Assets
{
    
        static class FlMath
        {
            static double[] Hai = { 0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 25000, 26000 };
            static double[] ai = { 340.3, 336.4, 332.5, 328.6, 324.6, 320.5, 316.4, 312.3, 308.0, 303.8, 299.5, 295.1, 295.1, 296.9 };
            static double[] Hroi = { 0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000, 22000, 24000, 26000, 28000, 30000, 32000, 34000, 36000, 38000, 40000, 50000 };
            static double[] roi = { 0.1249, 0.1134, 0.1027, 0.09273, 0.08355, 0.07508, 0.06728, 0.06012, 0.05355, 0.04755, 0.04208, 0.0371, 0.03169, 0.02706, 0.02311, 0.01974, 0.01686, 0.01440, 0.0123, 0.0105, 0.00897, 0.00654, 0.00477, 0.00348, 0.00254, 0.00185, 0.001365, 0.000996, 0.000727, 0.000530, 0.000387, 0.000080 };
            public readonly static double g = 9.81;
        
            public static double a(double H, double V) // проабгрейдить изменив температуру
            {
                return V / Interp(Hai, ai, UseEndValue(Hai, H));
            }

            public static double q(double H, double V)
            {
                return 0.5 * V * V * Interp(Hroi, roi, UseEndValue(Hroi, H)) * g;
            }

            public static double Ve(double H, double V)
            {
                return Math.Sqrt(8.0064 * Interp(Hroi, roi, UseEndValue(Hroi, H))) * V * 3.6;
            }

            public static void ToSvSK(double[] vector, double ugol)
            {
            double x0, y0, cosT, sinT;
            x0 = vector[0]; y0 = vector[1];
            cosT = Math.Cos(ugol / 180 * Math.PI);
            sinT = Math.Sin(ugol / 180 * Math.PI);
            vector[0] = cosT * x0 + sinT * y0;
            vector[1] = -sinT * x0 + cosT * y0;
            }
        #region // ИНТЕРПОЛЯЦИЯ
        public static double Interp(double[] Xm, double[] Ym, double xi) // функция одномерной интерполяции
            {
                double yi; // результат интерполяции

                if (xi >= Xm[Xm.Length - 1]) { yi = Ym[Ym.Length - 1]; } // превышение верхней границы массива аргументов берем последнее значение функции
                else if (xi <= Xm[0]) { yi = Ym[0]; } // значение точки меньше нижней границы массива аргументов берем первое значение функции
                else
                {
                    for (int i = 0; i < Xm.Length - 1; i++)
                    {
                        if (Xm[i] <= xi && xi < Xm[i + 1])
                        {
                            yi = Ym[i] + (Ym[i + 1] - Ym[i]) / (Xm[i + 1] - Xm[i]) * (xi - Xm[i]);
                            return yi;
                        }
                    }
                    yi = 0;
                }
                return yi;
            }

            public static double Interp(double[] Xm1, double[] Xm2, double[,] Ym, double x1i, double x2i) // функция двухмерной интерполяции
            {
                double y;
                x1i = UseEndValue(Xm1, x1i); x2i = UseEndValue(Xm2, x2i);

                for (int i = 0; i < Xm1.Length - 1; i++)
                {
                    if (Xm1[i] <= x1i && x1i <= Xm1[i + 1])
                    {
                        for (int j = 0; j < Xm2.Length - 1; j++)
                        {
                            if (Xm2[j] <= x2i && x2i <= Xm2[j + 1])
                            {
                                double k = (Xm1[i + 1] - Xm1[i]) * (Xm2[j + 1] - Xm2[j]);

                                y = Ym[i, j] * (Xm1[i + 1] - x1i) * (Xm2[j + 1] - x2i) / k +
                                    Ym[i + 1, j] * (x1i - Xm1[i]) * (Xm2[j + 1] - x2i) / k +
                                    Ym[i, j + 1] * (Xm1[i + 1] - x1i) * (x2i - Xm2[j]) / k +
                                    Ym[i + 1, j + 1] * (x1i - Xm1[i]) * (x2i - Xm2[j]) / k;
                                return y;
                            }

                        }
                    }

                }

                return 0;
            }

            public static double Interp(double[] Xm1, double[] Xm2, double[] Xm3, double[,,] Ym, double xi1, double xi2, double xi3) // функция трехмерной интерполяции
            {
                double y; // результат
                xi1 = UseEndValue(Xm1, xi1); xi2 = UseEndValue(Xm2, xi2); xi3 = UseEndValue(Xm3, xi3); //ограничиваем максимальные значения 
                                                                                                       // Console.WriteLine("xi1= " + xi1 + " Xm1[0]= " + Xm1[0]);
                for (int i = 0; i < Xm1.Length - 1; i++)
                {
                    if (Xm1[i] <= xi1 && xi1 <= Xm1[i + 1]) //проверяем вхождение по первому аргументу
                    {
                        for (int j = 0; j < Xm2.Length - 1; j++)
                        {
                            if (Xm2[j] <= xi2 && xi2 <= Xm2[j + 1]) //проверяем вхождение по второму аргументу
                            {
                                for (int k = 0; k < Xm3.Length - 1; k++)
                                {
                                    if (Xm3[k] <= xi3 && xi3 <= Xm3[k + 1]) //проверяем вхождение по третьему аргументу
                                    {
                                        double xd, yd, zd, c00, c01, c10, c11, c0, c1;
                                        xd = (xi1 - Xm1[i]) / (Xm1[i + 1] - Xm1[i]);
                                        yd = (xi2 - Xm2[j]) / (Xm2[j + 1] - Xm2[j]);
                                        zd = (xi3 - Xm3[k]) / (Xm3[k + 1] - Xm3[k]);

                                        c00 = Ym[i, j, k] * (1 - xd) + Ym[i + 1, j, k] * xd;
                                        c01 = Ym[i, j, k + 1] * (1 - xd) + Ym[i + 1, j, k + 1] * xd;
                                        c10 = Ym[i, j + 1, k] * (1 - xd) + Ym[i + 1, j + 1, k] * xd;
                                        c11 = Ym[i, j + 1, k + 1] * (1 - xd) + Ym[i + 1, j + 1, k + 1] * xd;

                                        c0 = c00 * (1 - yd) + c10 * yd;
                                        c1 = c01 * (1 - yd) + c11 * yd;

                                        y = c0 * (1 - zd) + c1 * zd;
                                        return y;
                                    }
                                }
                            }
                        }
                    }
                }
                return 0;
            }
        #endregion
        public static double UseEndValue(double[] X, double x)
            {
                if (x >= X[X.Length - 1]) { x = X[X.Length - 1]; }
                else if (x <= X[0]) { x = X[0]; }
                else { return x; }
                return x;
            }
        }

    static class Prodol
    {
        #region // ИНИЦИАЛИЗАЦИЯ ПЕРЕМЕННЫХ БАЗЫ
        public static readonly double[][] argM = new double[22][];
        public static readonly double[][] fcya = new double[22][];
        public static readonly double[][] fmza = new double[22][];
        public static readonly double[][] fmzFi = new double[22][];

        public static readonly double[][] demfCase = new double[22][]; // сопоставление номера самолета и его масивов демфирования и mzDB
        public static readonly double[][] argMdemf = new double[16][];
        public static readonly double[][] fmz_at = new double[16][];
        public static readonly double[][] fmz_wz = new double[16][];

        public static readonly double[] L = new double[22];
        public static readonly double[] S = new double[22];
        public static readonly double[] ba = new double[22];
        public static readonly double[] m0 = new double[22];
        public static readonly double[] Jz0 = new double[22];
        public static readonly int[] KtoK = new int[22];

        public static string[] PlansNames = new string[22];
        #endregion

        public static void WakeBasa()
        {
            //ИМЕНА САМОЛЕТОВ
            for (int kk = 0; kk < 1; kk++)
            {
                PlansNames[0] = "FighterPlane";
            }
            #region ГЕОМЕТРИЯ
            //РАЗМАХ КРЫЛА
            for (int kk = 0; kk < 1; kk++)
            {
                L[0] = 10.5;
            }

            //ПЛОЩАДЬ КРЫЛА
            for (int kk = 0; kk < 1; kk++)
            {
                S[0] = 50;
            }

            //ХОРДА КРЫЛА
            for (int kk = 0; kk < 1; kk++)
            {
                ba[0] = 5.7;
            }
#endregion
            #region //БАНК АЭРОДИНАМИКА
            // число М
            for (int kk = 0; kk < 1; kk++)
            {
                argM[0] = new double[] { 0.2, 0.6, 0.81, 1.2, 1.5, 1.7 };
            }
            // cy_a
            for (int kk = 0; kk < 1; kk++)
            {
                fcya[0] = new double[] { 0.050793575, 0.053519927, 0.06283023, 0.102436037, 0.077096351, 0.06153199 };
            }
            // mz_a
            for (int kk = 0; kk < 1; kk++)
            {
                fmza[0] = new double[] { -0.080802828, -0.088566275, -0.126035966, -0.1666957, -0.091361021, -0.075134367 };
            }
            // mz_DB
            for (int kk = 0; kk < 1; kk++)
            {
                fmzFi[0] = new double[] { -0.0019190912, -0.0021445574, -0.0024638309, -0.0020812476, -0.0017341263, -0.0015320706 };
            }

            // параметры сопоставления самолетов и их массивов демпфирования
            for (int kk = 0; kk < 1; kk++)
            {
                KtoK[0] = 4;
            }
            // число M для параметров демпфирования
            for (int kk = 0; kk < 1; kk++)
            {
                argMdemf[0] = new double[] { 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6 };
            }
            // mz_AT
            for (int kk = 0; kk < 1; kk++)
            {
                fmz_at[0] = new double[] { -2.2, -2.2, -2.2, -2.2, -2.2, -2, -1.85 };
            }
            //mz_WZ
            for (int kk = 0; kk < 1; kk++)
            {
                fmz_wz[0] = new double[] { -9.21, -9.21, -9.21, -9.21, -9.21, -9.21, -9.32 };
            }
            #endregion
            #region// начальные масса и момент инерции
            //масса
            for (int kk = 0; kk < 1; kk++)
            {
                m0[0] = 18000;
            }
            //момент инерции
            for (int kk = 0; kk < 1; kk++)
            {
                Jz0[0] = 267866;
            }
            #endregion
        }

        public static void sim(
        double dT,
        double m, double Jz, double DB,
        double V, double Kx, double Ky, int PlaneId,
        double windX, double windY,
        ref double alpha,
        ref double tang,
        ref double a_t,
        ref double wz,
        ref double Tet,
        ref double H,
        ref double Lx
            )
        {
            double a, q;


            a = FlMath.a(H, V);
            q = FlMath.q(H, V);

            // вычисление табличных коэффициентов
            AeroBank(PlaneId, a, out double mz_at, out double mz_a, out double mz_wz, out double mz_DB, out double cy_a);
            // вычисление РАЗМЕРНЫХ коэффициентов
            Kf_Pr(m, FlMath.g, q, S[PlaneId], ba[PlaneId], V, mz_at, mz_a, mz_wz, mz_DB, cy_a,
            out double Mz_at, out double Mz_wz, out double Mz_DB, out double Mz_a, out double Yao_a, out double Yao_V);
            // интегрирование
            PfProdInt(Mz_at, Mz_wz, Mz_DB, Mz_a, Yao_a, Yao_V, windX, windY, DB, V, m, Jz, dT,
                ref alpha, ref a_t, ref wz, ref Tet, ref H, ref Lx, ref tang, Kx, Ky);
        }

        public static void AeroBank(int PlaneId, double a, out double mz_at, out double mz_a, out double mz_wz, out double mz_DB, out double cy_a)
        {
            double a1, a2; // крайние значения по Маху таблиц данных
            int demfId;
            a1 = FlMath.UseEndValue(argM[PlaneId], a);
            cy_a = FlMath.Interp(argM[PlaneId], fcya[PlaneId], a1);
            mz_a = FlMath.Interp(argM[PlaneId], fmza[PlaneId], a1);
            mz_DB = FlMath.Interp(argM[PlaneId], fmzFi[PlaneId], a1);

            //fmz_at = new double[16][];
            // public static readonly double[][] fmz_wz

            demfId = KtoK[PlaneId];
            a2 = FlMath.UseEndValue(argMdemf[demfId], a);
            mz_at = FlMath.Interp(argMdemf[demfId], fmz_at[demfId], a2);
            mz_wz = FlMath.Interp(argMdemf[demfId], fmz_wz[demfId], a2);
        }

        public static void Kf_Pr(double m, double g, double q, double S, double ba, double V0,
            double mz_at, double mz_a, double mz_wz, double mz_DB, double cy_a,

            out double Mz_at,
            out double Mz_wz,
            out double Mz_DB,
            out double Mz_a,
            out double Yao_a,
            out double Yao_V
            ) // вычисляем для интегрирования коэффициенты
        {
            double Kz = q * S * ba;
            Mz_at = Kz * ba / V0 * mz_at;
            Mz_wz = Kz * ba / V0 * mz_wz;
            Mz_DB = Kz * mz_DB;
            Mz_a = Kz * mz_a;
            Yao_a = q * S * cy_a;
            Yao_V = 2 * m * g / V0;

        }

        public static void PfProdInt(
         double Mz_at, double Mz_wz, double Mz_DB, double Mz_a,
         double Yao_a, double Yao_V,
         double Windx, double Windy, double DB, double V0,
         double m, double Jz,
         double dt, // ШАГ интегрирования
         ref double alpha, ref double a_t, ref double wz, ref double Tet,
         ref double H, ref double Lx, ref double tang,
         double Kx, double Ky // параметры изменения скорости изменения координат
            )
        {
            double toRad=Math.PI/180;
            double toDegree = 1 / toRad;
            double aT = Windy / V0 * toDegree;
            double dTet_dt = (Yao_V * (-Windx) + Yao_a * (alpha + aT)) / (m * V0) * toDegree;

            //Math.Abs(Kx*Math.Cos(tang* toRad)+Ky*Math.Sin(tang*toRad))

            //Math.Abs(Kx * Math.Cos(Tet * toRad) + Ky * Math.Sin(Tet * toRad)) вариант с ограничением путем перехода из одних с.к. в другие(вышло не осень)
            //Math.Abs(-Kx * Math.Sin(Tet * toRad) + Ky * Math.Cos(Tet * toRad))
            double dX_dt =  Kx* (V0 * Math.Cos(Tet * toRad) - Windx * 1);
            double dH_dt =  Ky* (V0 * Math.Sin(Tet * toRad) - Windy * 1);

            double dwz_dt = (Mz_a * (alpha + aT) + Mz_at * a_t *toRad + Mz_wz * wz * toRad + Mz_DB * (DB + wz * 0.0 + 0.0 * alpha)) / Jz * toDegree;

            // интегрирование эйлером
            wz = wz + dt * dwz_dt;
            tang = tang + dt * wz;
            a_t = wz - dTet_dt;

            if (tang >= 360)
            { tang = tang - 360; Tet = Tet - 360;} // tang = FlMath.UseEndValue(radTang, tang); }
            else if (tang <= -360)
            { tang = tang + 360; Tet = Tet + 360; } // tang = FlMath.UseEndValue(radTang, tang); }

            Tet = Tet + dt * dTet_dt;

            alpha = tang - Tet;
            Lx = Lx + dt * dX_dt;
            H = H + dt * dH_dt;
        }

        static void riter(string Name1, double value1,
            string Name2, double value2,
            string Name3, double value3,
            string Name4, double value4,
            string Name5, double value5,
            string Name6, double value6,
            string Name7, double value7,
            string Name8, double value8,
            string Name9, double value9,
            string Name10, double value10
            )
        {

            string writePath = @"C:\Users\Паааппв\Desktop\GameDev_Analize\anal1.txt";

            try
            {
                using (StreamWriter sw = new StreamWriter(writePath, true, System.Text.Encoding.Default))
                {
                    sw.WriteLine(Name1 + " "+value1.ToString() + " " +
                        Name2 + " " + value2.ToString() + " " +
                        Name3 + " " + value3.ToString() + " " +
                        Name4 + " " + value4.ToString() + " " +
                        Name5 + " " + value5.ToString() + " " +
                        Name6 + " " + value6.ToString() + " " +
                        Name7 + " " + value7.ToString() + " " +
                        Name8 + " " + value8.ToString() + " " +
                        Name9 + " " + value9.ToString() + " "+
                        Name10 + " " + value10.ToString() + " ");
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
        }


    }
}
