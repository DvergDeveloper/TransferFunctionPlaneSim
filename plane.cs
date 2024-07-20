using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using Assets;

public class plane : MonoBehaviour
{

    public Transform ThisPlaneTr;
    public Transform DBtransform;
    public float gTime;
    public int PlaneId = 0;
    public double m, Jz;
    double dt;
    public int play = 1;
    public double DB = 0;
    public double V = 100;
    public double Vi = 50;

    public double Vmin = 50;  // (на высоте 0) задать массив в базе для всех самолетов и брать оттуда
    public double Vmax = 700; // (на высоте 0) задать массив в базе для всех самолетов и брать оттуда
    public double maxdV = 20; // максимальное изменение скорости

    public double Kx = 0.2; //коэффициент умножения Vx
    public double Ky = 0.2; // коэфициент умножения Vy
    public double tang_max = 30;
    private double windX =0;
    private double windY =0;
    public double alpha = 0;
    public double tang = 0;
    public double a_t = 0;
    public double wz = 0;
    public double Tet = 0;
    public double H = 4000;
    public double Lx = 0;
    public GameObject cam;

    double DBu = 0;
    double dVg = 0;// величина изменения скорости изза гравитации
    double ddVg = 1;// скорость изменения скорости из-за гравитации

    double DBspeed = 1;// изменение скорости отклонения руля высоты
    // управление
    Gyroscope gyroscope;
    public double[] DBmax=new double[] { -30, 30 };
    // Start is called before the first frame update
    void Start()
    {
        Prodol.WakeBasa();
        Initiate();
    }
    void Initiate()
    {
        m = Prodol.m0[PlaneId];
        Jz = Prodol.Jz0[PlaneId];
        ThisPlaneTr = GetComponent<Transform>();
        // берем гироскоп устройства
        Gyroscope gyroscope = Input.gyro;
    }

    void FixedUpdate()
    {   
        Quaternion phoneOrientation;
        UserInput();

        // ветер, управление через клавиатуру
        if (Input.GetKey(KeyCode.RightArrow))
        { windY = 10; }
        else if (Input.GetKey(KeyCode.LeftArrow))
        { windY = -10; }
        else { windY = 0; }

        gyroscope =  Input.gyro;
        double Bx = 2 * (double)(gyroscope.attitude.x * gyroscope.attitude.w); // возможно это не то  угол что нужен

        
        dt = Time.fixedDeltaTime;
        MehSim();
        GravityDeltaV();
        ChangeVengine();

        // моделирование (добавь обратную связь по wz и может a как входной сигнал в функцию sim)
        Prodol.sim(dt, m, Jz, DB, V+dVg, Kx, Ky, PlaneId, windX, windY,
        ref alpha, ref tang, ref a_t, ref wz, ref Tet, ref H, ref Lx);
        
        // задание ротации
        ThisPlaneTr.position = new Vector3((float)Lx, (float)(H - 4000), 0);
        ThisPlaneTr.rotation = Quaternion.Euler(0, 180, (float)-tang); //new Quaternion(0, (float)-tang + 180, 0, 0);
        DBtransform.localRotation = Quaternion.Euler(0, 0, (float)-DB);
        // задание позиции камеры
        cam.transform.position = new Vector3(ThisPlaneTr.position.x, ThisPlaneTr.position.y,-10);
    }
    void MehSim()
    {
        DB = DB + (DBu - DB) * dt * DBspeed; // привод РУЛЯ ВЫСОТЫ
    }
    void GravityDeltaV()
    {
        dVg = dVg + (-Math.Sin(tang / 180 * Math.PI)*30 - dVg) * dt*ddVg;
    }
    void ChangeVengine()
    {
        double changeV;
        changeV = (Vi - V) > 0 ? Math.Min((Vi - V), maxdV) : Math.Max((Vi - V), -maxdV);
        V = V + changeV * dt;
    }
    void OtkazSU()
    {
        DBspeed = 0.5;
    }
    void UserInput()
    {
        // РВ управление через клавиатуры
        gTime = Time.fixedTime;
        if (Input.GetKey(KeyCode.DownArrow))
        { DBu = 50; }
        else if (Input.GetKey(KeyCode.UpArrow))
        { DBu = -50; }
        else { DBu = 0; }
    }

    void OnCollisionEnter2D(Collision2D col)
    {
        if (col.gameObject.CompareTag("ground"))
        {
            double contactX, contactY, A, B, C;
            double Xmin, Ymin, Xmin90, Ymin90; // координаты точки на линии, ближайшей к ц.м. (минимальное расстояние) 

            ContactPoint2D[] contactPoints2D = new ContactPoint2D[20];
            this.GetComponent<Collider2D>().GetContacts(contactPoints2D);

            foreach (ContactPoint2D c in contactPoints2D)
            {
                if (c.collider is null) { break; }
                else
                {
                    double rad_stuck;// радиус точки касания для определения скорости вращения самолета
                                     
                    contactX = c.point.x; contactY = c.point.y; // ОПРЕДЕЛЕНИЕ ЗНАКА ВРАЩЕНИЯ ОТ ТОЧКИ КАСАНИЯ С ОПОРОЙ
                    double[] cont_in_SV = new double[2] { contactX - ThisPlaneTr.position.x, contactY - ThisPlaneTr.position.y };
                    double[] norm_in_SV = new double[2] { cont_in_SV[0] + c.normal.x * 2, cont_in_SV[1] + c.normal.y * 2 };// точка контакта смещенная на удвоенную норму
                    FlMath.ToSvSK(cont_in_SV, this.GetComponent<Transform>().rotation.eulerAngles.z);
                    FlMath.ToSvSK(norm_in_SV, this.GetComponent<Transform>().rotation.eulerAngles.z);

                    // центр масс в связанной с.к. = 0, 0
                    rad_stuck = Math.Sqrt(cont_in_SV[1] * cont_in_SV[1] + cont_in_SV[0] * cont_in_SV[0]);// радиус касания относительно центра масс 
                    // коэффициенты прямой нормали(из точки столкновения) в связанной с.к.
                    A = cont_in_SV[1] - norm_in_SV[1];
                    B = norm_in_SV[0] - cont_in_SV[0];
                    C = cont_in_SV[0] * norm_in_SV[1] - cont_in_SV[1] * norm_in_SV[0];
                    // расстояние от центра масс до прямой по координатам
                    Xmin = -A * C / (A * A + B * B);
                    Ymin = -B * C / (A * A + B * B);
                    // (поворачиваем перпендикуляр на 90 градусов)так как центр масс в начале координат  точки Xmin, Ymin - вектор перпендикуляра
                    Xmin90 = Ymin * (-1);
                    Ymin90 = Xmin;
                    //Вектор нормали в связанной с.к.
                    norm_in_SV[0] = norm_in_SV[0] - cont_in_SV[0];
                    norm_in_SV[1] = norm_in_SV[1] - cont_in_SV[1];

                    //глобальные скорости
                    double new_teta;
                    double teta = (tang - alpha) * Math.PI / 180;
                    double Vx = V * Math.Cos(teta);
                    double Vy = V * Math.Sin(teta);
                    //
                    double kinrad = rad_stuck * Math.Abs(V/10);
                    if (Math.Sign(Math.Sign(Xmin90 * norm_in_SV[0]) + Math.Sign(Ymin90 * norm_in_SV[1])) > 0)
                    { //кобрирование
                        if (wz < 0) { wz = Math.Max(-wz / 2, kinrad); }
                        else if (wz > 0) { wz = Math.Max(1.5 * wz, kinrad); }
                        else { wz = rad_stuck * Math.Abs(V); }
                        ChangesOnCollision(c.normal.x, c.normal.y, ref Vx, ref Vy);
                    }
                    else if (Math.Sign(Math.Sign(Xmin90 * norm_in_SV[0]) + Math.Sign(Ymin90 * norm_in_SV[1])) < 0)
                    { //пикирование
                        if (wz > 0) { wz = Math.Min(-wz / 2, -kinrad); }
                        else if (wz < 0) { wz = Math.Min(1.5 * wz, -kinrad); }
                        else { wz = -rad_stuck * Math.Abs(V); }
                        ChangesOnCollision(c.normal.x, c.normal.y, ref Vx, ref Vy);
                    }
                    else
                    { //нет момента только сила
                        ChangesOnCollision(c.normal.x, c.normal.y, ref Vx, ref Vy);
                    }

                    //  ПЕРЕВОД ИЗ ОБЫЧНЫХ КООРДИНАТ В УГЛОВЫЕ
                    if (Vx > 0 && Vy >= 0) { new_teta = Math.Atan(Vy / Vx) * 180 / Math.PI; }
                    else if (Vx > 0 && Vy < 0) { new_teta = (Math.Atan(Vy / Vx) + 2 * Math.PI) * 180 / Math.PI; }
                    else if (Vx < 0) { new_teta = (Math.Atan(Vy / Vx) + Math.PI) * 180 / Math.PI; }
                    else if (Vx == 0 && Vy > 0) { new_teta = 90; }
                    else if (Vx == 0 && Vy < 0) { new_teta = 270; }
                    else { new_teta = 0; }//if (Vx == 0 && Vy == 0) 

                    if (new_teta > 180) { new_teta = new_teta-360; }
                    alpha = tang - new_teta;
                    Tet = new_teta;
                    V = Math.Sqrt(Vx * Vx + Vy * Vy);
                    return;
                }
            }
        }
    }

    void ChangesOnCollision(double NormX, double NormY, ref double Vx, ref double Vy, double ky = 1d, double kx = 1d)// метод изменяющий силы, вектора и угловые скорости
    {
        double NormX90, NormY90;
        double A11, A12, A21, A22;
        double M11, M12, M21, M22;
        double Ad11, Ad12, Ad21, Ad22;
        double AdT11, AdT12, AdT21, AdT22;
        double Arev11, Arev12, Arev21, Arev22;
        double Vx_, Vy_, detA;

        // поварачиваем нормаль на 90 градусов против часовой
        NormX90 = NormY; NormY90 = -NormX;
        // создаем базисную матрицу 
        A11 = NormX90; A12 = NormX; A21 = NormY90; A22 = NormY;
        // рассчитываем определитель матрицы
        detA = A11 * A22 - A12 * A21;

        if (detA != 0)
        {
            M11 = A22; M12 = A21; M21 = A12; M22 = A11;

            Ad11 = M11; Ad12 = -M12; Ad21 = -M21; Ad22 = M22;

            AdT11 = Ad11; AdT12 = Ad21; AdT21 = Ad12; AdT22 = Ad22;
            // рассчитываем коэффициенты для перевода векторов в с.к. столкновения
            Arev11 = AdT11 / detA; Arev12 = AdT12 / detA;
            Arev21 = AdT21 / detA; Arev22 = AdT22 / detA;
            // переводим скорости в с.к. столкновения
            Vx_ = Arev11 * Vx + Arev12 * Vy;
            Vy_ = Arev21 * Vx + Arev22 * Vy;
            // ихменяем скорости делая скорость в плоскости столкновения всегда положительной
            Vy_ = Math.Abs(Vy_) * ky;
            Vx_ = Vx_ * kx;
            // переводим скорости обратно в земную с.к.
            Vx = A11 * Vx_ + A12 * Vy_;
            Vy = A21 * Vx_ + A22 * Vy_;
        }
    }

}
