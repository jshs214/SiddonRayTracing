#include "Siddon.h"

Siddon::Siddon(int matrix_x, int matrix_y, int matrix_z, double dx,
             double dy, double dz, double s_X, double s_Y, double s_Z, double d_X, double d_Y, double d_Z)
{
    m_x = matrix_x;         //set matrix(x*y*z)
    m_y = matrix_y;
    m_z = matrix_z;

    m_nx = m_x + 1;         //set nx, ny, nz
    m_ny = m_y + 1;
    m_nz = m_z + 1;

    m_dx = dx;              //set dx, dy, dz
    m_dy = dy;
    m_dz = dz;

    m_sourceX = s_X;        //set source(x,y,z)
    m_sourceY = s_Y;
    m_sourceZ = s_Z;

    m_detectorX = d_X;      //set detector(x,y,z)
    m_detectorY = d_Y;
    m_detectorZ = d_Z;

    m_xPlane1 = -dx * (matrix_x / 2.0);     // x,y,z의 plane(1)
    m_yPlane1 = -dy * (matrix_y / 2.0);
    m_zPlane1 = -dz * (matrix_z / 2.0);

    Eq_3();
    Eq_4();
    Eq_5();
    Eq_6();
    Eq_7();
    Eq_8();
    Eq_9();
    Eq_1011();
    Eq_1213();      

    //PrintInfo();
    //Eq_14(rho);
}

void Siddon::Eq_3()
{
    for (int i = 1; i <= m_nx; i++) {      //set xPlane(i)
        xPlane[i] = m_xPlane1 + (i - 1) * m_dx;
    }

    for (int i = 1; i <= m_ny; i++) {      //set yPlane(i)
        yPlane[i] = m_yPlane1 + (i - 1) * m_dy;
    }
    
    for (int i = 1; i <= m_nz; i++) {      //set zPlane(i)
        zPlane[i] = m_zPlane1 + (i - 1) * m_dz;
    }
    
}

void Siddon::Eq_4()
{
    if (m_detectorX - m_sourceX != 0) {                   //set alphaX[1], [Nx]
        alphaX[1] = (xPlane[1] - m_sourceX) / (m_detectorX - m_sourceX);
        alphaX[m_nx] = (xPlane[m_nx] - m_sourceX) / (m_detectorX - m_sourceX);
    }
    if (m_detectorY - m_sourceY != 0) {                  // set alphaY[1], [Ny]
        alphaY[1] = (yPlane[1] - m_sourceY) / (m_detectorY - m_sourceY);
        alphaY[m_ny] = (yPlane[m_ny] - m_sourceY) / (m_detectorY - m_sourceY);
    }
    if (m_detectorZ - m_sourceZ != 0) {                  // set alphaZ[1], [Nz]
        alphaZ[1] = (zPlane[1] - m_sourceZ) / (m_detectorZ - m_sourceZ);
        alphaZ[m_nz] = (zPlane[m_nz] - m_sourceZ) / (m_detectorZ - m_sourceZ);
    }
    
}
void Siddon::Eq_5() //Eq(5)
{
    double aX_Min, aX_Max, aY_Min, aY_Max, aZ_Min, aZ_Max;
    
    // set max[aX(1),aX(Nx)], min[aX(1),aX(Nx)]
    if(alphaX[1] < alphaX[m_nx]) {  
        aX_Min = alphaX[1];
        aX_Max = alphaX[m_nx];
    }
    else{
        aX_Min = alphaX[m_nx];
        aX_Max = alphaX[1];
    }
    // set max[aY(1),aY(Ny)], min[aY(1),aY(Ny)]
    if(alphaY[1] < alphaY[m_ny]){
        aY_Min = alphaY[1];
        aY_Max = alphaY[m_ny];
    }
    else{
        aY_Min = alphaY[m_ny];
        aY_Max = alphaY[1];    
    }
    // set max[aZ(1),aZ(Nz)], min[aZ(1),aZ(Nz)]
    if(alphaZ[1] < alphaZ[m_nz]){
        aZ_Min = alphaZ[1];
        aZ_Max = alphaZ[m_nz];
    }
    else{
        aZ_Min = alphaZ[m_nz];
        aZ_Max = alphaZ[1];    
    }

    //alphaMin, alphaMax 를 구하기 위한 Function
    alphaMin = maximum4(0, aX_Min, aY_Min, aZ_Min);
    alphaMax = minimum4(1, aX_Max, aY_Max, aZ_Max);
}
void Siddon::Eq_6() //Eq(6)         올림 내림 연산 (a=min axmin 과 같은 case) 추가 수정 필요
{
    // set iMin, iMax
    if(m_detectorX-m_sourceX > 0){
        iMin = m_nx - (xPlane[m_nx] - alphaMin*(m_detectorX-m_sourceX) - m_sourceX) / m_dx;
        iMax = 1 + (m_sourceX + alphaMax*(m_detectorX-m_sourceX) - xPlane[1]) / m_dx;
    }
    else if(m_detectorX-m_sourceX < 0){
        iMin = m_nx - (xPlane[m_nx] - alphaMax*(m_detectorX-m_sourceX) - m_sourceX ) / m_dx;
        iMax = 1 + (m_sourceX + alphaMin*(m_detectorX-m_sourceX) - xPlane[1]) / m_dx;
    }
    // set jMin, jMax
    if(m_detectorY-m_sourceY > 0){
        jMin = m_ny - (yPlane[m_ny] - alphaMin*(m_detectorY-m_sourceY) - m_sourceY) / m_dy;
        jMax = 1 + (m_sourceY + alphaMax*(m_detectorY-m_sourceY) - yPlane[1]) / m_dy;
    }
    else if(m_detectorY-m_sourceY < 0){
        jMin = m_ny - (yPlane[m_ny] - alphaMax*(m_detectorY-m_sourceY) - m_sourceY ) / m_dy;
        jMax = 1 + (m_sourceY + alphaMin*(m_detectorY-m_sourceY) - yPlane[1]) / m_dy;
    }
    // set kMin, kMax
    if(m_detectorZ-m_sourceZ > 0){
        kMin = m_nz - (zPlane[m_nz] - alphaMin*(m_detectorZ-m_sourceZ) - m_sourceZ) / m_dz;
        kMax = 1 + (m_sourceZ + alphaMax*(m_detectorZ-m_sourceZ) - zPlane[1]) / m_dz;
    }
    else if(m_detectorZ-m_sourceZ < 0){
        kMin = m_nz - (zPlane[m_nz] - alphaMax*(m_detectorZ-m_sourceZ) - m_sourceZ ) / m_dz;
        kMax = 1 + (m_sourceZ + alphaMin*(m_detectorZ-m_sourceZ) - zPlane[1]) / m_dz;
    }

    iMin = m_ceil(iMin);
    iMax = m_floor(iMax);
    jMin = m_ceil(jMin);
    jMax = m_floor(jMax);
    kMin = m_ceil(kMin);
    kMax = m_floor(kMax);
}
void Siddon::Eq_7() //Eq(7) alphaX(Min~Max 나 Max~Min)범위에 대한 alpha 집합 생성 과정은 Eq(8)에서.
                    // if(X2-X1) < 0 에서 수정사항 있을 수도 있음.
{
    int temp;

    //set alphaX(iMin) ~ alphaX(iMax)
    for(int i = 1; i <= m_nx; i ++){      //alphaX[i]값을 모두 계산
        if(i != 1 && i != m_nx){
            alphaX[i] = alphaX[i-1] + (m_dx/(m_detectorX-m_sourceX) );
        }
    }
    //set alphaX(iMax) ~ alphaX(iMin)
    if(m_detectorX-m_sourceX < 0){      //set iMin~iMax를 swap하는 과정 
        int count =0;
        for (int i = iMin; i <= iMax/2; i++) {
            if(i != 1 && i != m_nx){
                temp = alphaX[i];
                alphaX[i] = alphaX[iMax - count];
                alphaX[iMax - count] = temp;
                count++;
            }
        }
    }

    //set alphaY(jMin) ~ alphaY(jMax)
    for(int i = 1; i <= m_ny; i ++){      //alphaY[i]값을 모두 계산
            if(i != 1 && i != m_ny){
                alphaY[i] = alphaY[i-1] + (m_dy/(m_detectorY-m_sourceY) );
            }
    }
    //set alphaY(jMax) ~ alphaY(jMin)
    if(m_detectorY-m_sourceY < 0){      //set jMin~jMax를 swap하는 과정   
        int count =0;
        for (int i = jMin; i <= jMax/2; i++) {
            if(i != 1 && i != m_ny){
                temp = alphaY[i];
                alphaY[i] = alphaY[jMax - count];
                alphaY[jMax - count] = temp;
                count++;
            }
        }
    }
    
    //set alphaZ(kMin) ~ alphaZ(kMax)
    for(int i = 1; i <= m_nz; i ++){      //alphaZ[i]값을 모두 계산
            if(i != 1 && i != m_nz){
                alphaZ[i] = alphaZ[i-1] + (m_dz/(m_detectorZ-m_sourceZ) );
            }
    }
    //set alphaY(jMax) ~ alphaY(jMin)
    if(m_detectorZ-m_sourceZ < 0){      //set jMin~jMax를 swap하는 과정   
        int count =0;
        for (int i = kMin; i <= kMax/2; i++) {
            if(i != 1 && i != m_nz){
                temp = alphaZ[i];
                alphaZ[i] = alphaZ[kMax - count];
                alphaZ[kMax - count] = temp;
                count++;
            }
        }
    }
        
}
void Siddon::Eq_8() //Eq(8)
{
    //alpha vector 에 alphaX, alphaY, alphaZ, alphaMin, alphaMax 값 삽입 후 asc 정렬
    for (auto iter = alphaX.begin(); iter != alphaX.end(); ++iter) {
        if(iter->first >= iMin && iter->first<=iMax)       //alphaX의 Index 범위에 맞는값 alpha 집합에 삽입
            alpha.push_back(iter->second);
    }
    for (auto iter = alphaY.begin(); iter != alphaY.end(); ++iter) {
        if(iter->first >= jMin && iter->first<=jMax)       //alphaY의 Index 범위에 맞는값 alpha 집합에 삽입
            alpha.push_back(iter->second);
    }
    for (auto iter = alphaZ.begin(); iter != alphaZ.end(); ++iter) {
        if(iter->first >= kMin && iter->first<=kMax)       //alphaY의 Index 범위에 맞는값 alpha 집합에 삽입
            alpha.push_back(iter->second);
    }

    alpha.push_back(alphaMin);
    alpha.push_back(alphaMax);

    sort(alpha.begin(), alpha.end());
}
void Siddon::Eq_9() //Eq(9)
{
     n = (iMax - iMin + 1) + (jMax - jMin + 1) + (kMax - kMin + 1) + 1;
}
void Siddon::Eq_1011()    //Eq(10) 11과정 포함.
{
    int m;     //siddon 논문에서의 m, (m = 1, ... , n)

    d12 = m_sqrt( m_pow((m_detectorX-m_sourceX),2) + m_pow((m_detectorY-m_sourceY),2) + m_pow((m_detectorZ-m_sourceZ),2));
    
    for(auto iter = alpha.begin()+1; iter != alpha.end(); ++iter){        //alpha(1), ... , alpha(n)
        m = distance(alpha.begin(), iter);
        if(alpha[m] - alpha[m-1] < 0.00000000001){        //부동소수점 예외 처리
            lm[m] = 0;
        }
        else
            lm[m] = d12*(alpha[m] - alpha[m-1]);
    }
}
void Siddon::Eq_1213()    //Eq 13은 수식안에 바로 적용하였음.
{    
    int m;     //siddon 논문에서의 m, (m = 1, ... , n)

    for(auto iter = alpha.begin()+1; iter != alpha.end(); ++iter){        //alpha(1), ... , alpha(n)
        m = distance(alpha.begin(), iter);
        double alphaMid = (alpha[m]+alpha[m-1])/2;
        indexI[m] = m_floor( 1 + ((m_sourceX + alphaMid * (m_detectorX-m_sourceX) - m_xPlane1  ) / m_dx) );
        indexJ[m] = m_floor( 1 + ((m_sourceY + alphaMid * (m_detectorY-m_sourceY) - m_yPlane1  ) / m_dy) );
        indexK[m] = m_floor( 1 + ((m_sourceZ + alphaMid * (m_detectorZ-m_sourceZ) - m_zPlane1  ) / m_dz) );
    }
}
double Siddon::Eq_14(int ***rho)
{
    double sum=0.0;     //Eq(14) 첫번째 수식
    int count = 0;

    for(int m = 1; m <= n; m ++){               //Eq(14) 시그마 수식
        if(indexI[m] <= 0 || indexJ[m] <= 0 || indexK[m] <= 0 ||
             indexI[m] > m_x || indexJ[m] > m_y || indexK[m] > m_z || lm[m] ==0){
            continue;
        }
        else{
            int i = (int)indexI[m] - 1;     //index 는 1부터 시작하나,
            int j = (int)indexJ[m] - 1;     //rho의 좌표는 0부터 시작하므로 i, j, k 에 -1
            int k = (int)indexK[m] - 1;
            sum += lm[m] * rho[i][j][k];
            //cout << endl <<m <<", lm[m] : "<< lm[m] <<", rho"<<"[" <<i + 1 <<"]["<<j + 1 <<"][" << k +1 <<"] ="<<rho[i][j][k] 
            //        << ", l[m] * rho[i][j][k] : "<< lm[m] * rho[i][j][k] << " -> sum : "<<sum;
        }
    }
    //cout << " source : ("<<m_sourceX<<","<< m_sourceY <<") -> "<< "Detector : ("<<m_detectorX<<","<< m_detectorY <<") = " << sum <<endl;

    return sum;
}
void Siddon::PrintInfo()
{
    cout << "========parameter========" << endl;
    cout << "Matrix(x*y*z) : " << m_x << " " << m_y <<" "<<m_z << endl;
    cout << "dx, dy, dz : " << m_dx << " " << m_dy << " " << m_dz << endl;
    cout << "Nx, Ny, Nz : " << m_nx << " " << m_ny <<" " << m_nz << endl;
    cout << "source (x, y, z) : (" << m_sourceX << ", " << m_sourceY <<", "<< m_sourceZ << ")" << endl;
    cout << "detector (x, y, z) : (" << m_detectorX << ", " << m_detectorY << ", " << m_detectorZ <<")" << endl;

    //Eq(3)
    cout << "========Eq(3)========" << endl;
    for (auto iter = xPlane.begin(); iter != xPlane.end(); ++iter) {
        cout << "xPlane(" << iter->first << ") : " << iter->second << "   ";
    }
    cout << endl;
    for (auto iter = yPlane.begin(); iter != yPlane.end(); ++iter) {
        cout << "yPlane(" << iter->first << ") : " << iter->second << "   ";
    }
    cout << endl;
    for (auto iter = zPlane.begin(); iter != zPlane.end(); ++iter) {
        cout << "zPlane(" << iter->first << ") : " << iter->second << "   ";
    }
    cout << endl;

    //Eq(4)
    cout << "========Eq(4)========" << endl;
    cout << "Eq(6),(7) ...." << endl;
    cout << " alphaX[1] , alphaX[Nx] : " << alphaX[1] << ", " << alphaX[m_nx] << endl;
    cout << " alphaY[1] , alphaY[Ny] : " << alphaY[1] << ", " << alphaY[m_ny] << endl;
    cout << " alphaZ[1] , alphaZ[Nz] : " << alphaZ[1] << ", " << alphaZ[m_nz] << endl;
    
    //Eq(5)  
    cout << "========Eq(5)========" << endl;
    cout << "alpha Min, Max : " << alphaMin << ", " << alphaMax <<endl;

    //Eq(6)
    cout << "========Eq(6)========" << endl;
    cout << "iMin, iMax : " << iMin << ", " << iMax << endl;
    cout << "jMin, jMax : " << jMin << ", " << jMax << endl;
    cout << "kMin, kMax : " << kMin << ", " << kMax << endl;
    
    //Eq(7)     alphaX, Y값이 맞는지 print
    cout << "========Eq(7)========" << endl;
    for(auto iter = alphaX.begin(); iter != alphaX.end(); ++iter){
        cout <<"alphaX[" << iter->first <<"]:" << iter->second <<"  ";
    }
    cout <<endl;
    for(auto iter = alphaY.begin(); iter != alphaY.end(); ++iter){
        cout <<"alphaY[" << iter->first <<"]:" << iter->second <<"  ";
    }
    cout <<endl;
    for(auto iter = alphaZ.begin(); iter != alphaZ.end(); ++iter){
        cout <<"alphaZ[" << iter->first <<"]:" << iter->second <<"  ";
    }
    cout <<endl;

    //Eq(8)
    cout << "========Eq(8)========" << endl;
    int count =0;
    for(auto iter = alpha.begin(); iter != alpha.end(); ++iter){
        cout << "alpha["<<count<<"] : " << *iter << "   ";
        count ++;
    }

    //Eq(9)
    //alpha 집합의 요소 수
    cout <<endl << "========Eq(9)========" << endl;
    cout << " n : " << n << endl;
    cout << " alpha(0), ... , alpha("<<n <<")" << endl;

    //Eq(10, 11) 
    cout << "========Eq(1011)========" << endl;
    cout << "d12 : "<< d12 << endl;
    for(auto iter = lm.begin(); iter != lm.end(); ++iter){
        cout <<"lm["<< iter->first <<"]: "<<iter->second << "  |  ";
    }
    
    //Eq(1213)
    cout <<endl;
    cout << "========Eq(1213)========" << endl;
    for(auto iter = indexI.begin(); iter != indexI.end(); ++iter){
        cout <<"indexI["<< iter->first <<"]: "<<iter->second << "  |  ";
    }
    cout <<endl;
    for(auto iter = indexJ.begin(); iter != indexJ.end(); ++iter){
        cout <<"indexJ["<< iter->first <<"]: "<<iter->second << "  |  ";
    }
    cout <<endl;
    for(auto iter = indexK.begin(); iter != indexK.end(); ++iter){
        cout <<"indexK["<< iter->first <<"]: "<<iter->second << "  |  ";
    }
    cout <<endl;
}
double Siddon::maximum4(int num1,double num2, double num3, double num4){
    double maximum = num1;

    if (num2 >= maximum) {
        maximum = num2;
    }
    if (num3 >= maximum) {
        maximum = num3;
    }
    if (num4 >= maximum) {
        maximum = num4;
    }
    return maximum;
}

	//4개의 정수 중 최솟값 구하는 메소드
double Siddon::minimum4(int num1,double num2, double num3, double num4){
    double minimum = num1;

    if (num2 <= minimum) {
        minimum = num2;
    }
    if (num3 <= minimum) {
        minimum = num3;
    }
    if (num4 <= minimum) {
        minimum = num4;
    }
    return minimum;
}

double Siddon::m_ceil(double n)
{
    int r = (int) n;

    if(n < 0 || n == (double) r)
        return r;
    else 
        return r+1;
}
double Siddon::m_floor(double n)
{
    double r;
    
    r = (int) n - ((n>0)? 0 : 1);

    return r;
}
double Siddon::m_pow(double n, double m)
{
    double nn = n; 
    if (m == 0) 
      return 1.0;
    for (int i = 0; i < (m-1); i++)
        n *= nn;
    return n;
}
double Siddon::m_sqrt(double n)
{
    double s=0;
    double t=0;
 
    s = n / 2;
    for (;s != t;)
    {
        t = s;
        s = ( (n / t) + t) / 2;
    }
    return s;
}