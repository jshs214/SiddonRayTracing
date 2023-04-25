#pragma once
#ifndef __SIDDON_H__
#define __SIDDON_H__

#include <map>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

/**
* siddon raytraicing
*/
class Siddon {
public:
	Siddon(int matrix_x, int matrix_y, int matrix_z, double dx,
             double dy, double dz, double s_X, double s_Y, double s_Z, double d_X, double d_Y, double d_Z);
             
    void PrintInfo();
	void Eq_3();
    void Eq_4();
    void Eq_5();
    void Eq_6();
    void Eq_7();
    void Eq_8();
    void Eq_9();
    void Eq_1011();
    void Eq_1213();
    double Eq_14(int ***rho);

    double maximum4(int, double, double, double);
    double minimum4(int, double, double, double);
    double m_ceil(double);
    double m_floor(double);
    double m_pow(double, double);
    double m_sqrt(double);

private:
    int m_x, m_y, m_z;                       //Matrix size (x*y*z)
    double m_dx, m_dy, m_dz, m_nx, m_ny, m_nz;
    double m_sourceX, m_sourceY, m_sourceZ, m_detectorX, m_detectorY, m_detectorZ;

    double m_xPlane1, m_yPlane1, m_zPlane1;
    map<int, double> xPlane, yPlane, zPlane;
    map<int, double> alphaX, alphaY, alphaZ;

    double alphaMax = -0.1, alphaMin = 1.1;
    double iMin, iMax, jMin, jMax, kMin, kMax ;
    vector<double> alpha;   // alphaMin ~ alphaMax ASC ????? ????
    int n;
    double d12;
    map<int, double> lm, indexI, indexJ, indexK;
};
#endif