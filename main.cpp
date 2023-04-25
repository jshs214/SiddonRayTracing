#include "Siddon.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

using namespace std;

map<int, pair<double, double> > detector_map;		//detector의 좌표를 저장하는 detector_map <Key, <y, z> > y축은 고정.

/* detector_Size(Width, Height) 홀수, 짝수에 맞게 detector_map 생성.
@param detector_pixelDistance, detector_X, detector_Z, detector Width, detectorHegiht 
*/ 
void setDetectorSize(double detectorPixel, double detector_Y, double detector_Z, int detectorWidth, int detectorHeight)
{
	//홀수 * 홀수 case
	if( detectorWidth %2 != 0 && detectorHeight %2 != 0){
		double width = detectorWidth/2;
		double height = detectorHeight/2;
		int offset =0;

		for(double i = -height; i <= height; i++){					//detector 크기 만큼 detector_map에 detector (x, y)저장
			for(double j = -width; j <= width; j++){
				detector_map[offset] = {detector_Y + (detectorPixel*j), detector_Z + (detectorPixel*i) };	//detector_map에 기존 detector 좌표 + dx*offset의 좌표 값 저장
				offset++;
			}
		}
	}
	//짝수 * 짝수 case
	else if(detectorWidth %2 == 0 && detectorHeight %2 == 0){
		double width = detectorWidth/2 - 0.5;
		double height = detectorHeight/2 - 0.5;
		int offset =0;

		for(double i = -height; i <= height; i++){					//detector 크기 만큼 detector_map에 detector (x, y)저장
			for(double j = -width; j <= width; j++){
				detector_map[offset] = {detector_Y + (detectorPixel*j), detector_Z + (detectorPixel*i) };	//detector_map에 기존 detector 좌표 + dx*offset의 좌표 값 저장
				offset++;
			}
		}
	}
	//홀수 * 짝수 case
	else if(detectorWidth%2 != 0 && detectorHeight %2 == 0){
		double width = detectorWidth/2;
		double height = detectorHeight/2 - 0.5;
		int offset =0;

		for(double i = -height; i <= height; i++){					//detector 크기 만큼 detector_map에 detector (x, y)저장
			for(double j = -width; j <= width; j++){
				detector_map[offset] = {detector_Y + (detectorPixel*j), detector_Z + (detectorPixel*i) };	//detector_map에 기존 detector 좌표 + dx*offset의 좌표 값 저장
				offset++;
			}
		}
	}
	//짝수 * 홀수 case
	else{
		double width = detectorWidth/2 - 0.5;
		double height = detectorHeight/2;
		int offset =0;

		for(double i = -height; i <= height; i++){					//detector 크기 만큼 detector_map에 detector (x, y)저장
			for(double j = -width; j <= width; j++){
				detector_map[offset] = {detector_Y + (detectorPixel*j), detector_Z + (detectorPixel*i) };	//detector_map에 기존 detector 좌표 + dx*offset의 좌표 값 저장
				offset++;
			}
		}
	}
}

// x축, y축 동시 회전 행렬 구현
void rotateXY(double& x, double& y, double& z, double theta) {
    double cosTheta = cos(theta);
    double sinTheta = sin(theta);
    double newX = x * cosTheta - y * sinTheta;
    double newY = x * sinTheta + y * cosTheta;
    x = newX;
    y = newY;
}

int main() 
{
	int matrix_x, matrix_y, matrix_z;										//Matrix x*y*z 크기 설정
	double dx, dy, dz, source_X, source_Y, source_Z, detector_X, detector_Y, detector_Z;	//dx, dy, source 좌표(x, y, z), detector 좌표(x, y)

    cout << "Input | Object Matrix (x,y,z): ";		// Matrix 크기
    cin >> matrix_x >> matrix_y >> matrix_z;
    cout << "Input | dx, dy, dz : ";					// x,y,z pixel_distance (0 > 실수)
    cin >> dx >> dy >> dz;
    cout << "Input | source (x, y, z) : ";				// source (x,y,z)
    cin >> source_X >> source_Y >> source_Z;
	cout << "Input | detector : (x, y, z) : ";			// detector (x,y,z)
    cin >> detector_X >> detector_Y >> detector_Z;    

	double projValue;		//projValue를 저장하기 위한 임시변수
	double s_X, s_Y, s_Z;	//회전하기위한 값 복사 source, detector (x,y,z 좌표)
	double d_X, d_Y, d_Z;

	double thetaXY = M_PI / 180.0;  // x축, y축 동시 회전 각도: 1도를 라디안으로 변환;	//x축, y축 동시 회전 각도

	int detectorWidth = 128, detectorHeight = 128;		//detector Size (x축, z축) y축 = 0(임시)
	double detectorPixel = 1;

    setDetectorSize(detectorPixel, detector_Y, detector_Z, detectorWidth, detectorHeight);		//detector Size에따라 detector_map에 detector 좌표 할당하는 function.

	FILE *fp;		//fwrite 하기위한 fp


	float *projData;							//detector_map key값의 ray와 object의 교차 거리
	projData = (float*) malloc(sizeof(float) * detectorWidth*detectorHeight);

	// Matrix의 3차원 동적 배열 rho 생성 및 메모리 할당 
    int*** rho = new int**[matrix_x];
    for (int i = 0; i < matrix_x; ++i) {
        rho[i] = new int*[matrix_y];
        for (int j = 0; j < matrix_y; ++j) {
            rho[i][j] = new int[matrix_z];
        }
    }

    // 3차원 배열에 값 저장
    for (int i = 0; i < matrix_x; ++i) {
        for (int j = 0; j < matrix_y; ++j) {
            for (int k = 0; k < matrix_z; ++k) {
                rho[i][j][k] = 1;
            }
        }
    }

	cout << endl <<"===================Start siddon's Ray-Traicing Algorithm===================" ;
	cout << endl <<"detector Size : ("<<detectorWidth <<"*"<<detectorHeight<<"), detectorPixelSize : 1" << endl;
	
	//n도씩 회전
	for(int n = 0; n <=360; n+=45){

		//detector size만큼 siddon algorithm 수행.
		for(auto iter = detector_map.begin(); iter != detector_map.end(); ++iter){		//detector 좌표를 저장하는 map 개수 만큼 siddon ray tracing 수행
			d_X = detector_X;	//detector_map key값의 detector(x, Z) 값
			d_Y = iter->second.first;
			d_Z = iter->second.second;
			s_X = source_X;
			s_Y = source_Y;
			s_Z = source_Z;

			string fileName = "./projectionData/" + to_string(n) + "_projectionData.raw";		//file path + n도 projData

			// file Write and Close
			if ((fp = fopen(fileName.c_str(), "wb")) == NULL) {
				cout <<"No such file" <<endl;
				return -1;
			}
			
			if(n != 0){	//회전 각도 설정
				rotateXY(d_X, d_Y, d_Z, n * thetaXY);		// source, detector 좌표 n * 1 각도 회전
				rotateXY(s_X, s_Y, s_Z, n * thetaXY);
			}

			Siddon siddon(matrix_x, matrix_y, matrix_z, dx, dy, dz, s_X, s_Y, s_Z, d_X, d_Y, d_Z);		//siddon algorithm projData를 얻기 위한 function 수행
			projValue = siddon.Eq_14(rho);																				//projData를 구하는 수식인 Eq_14 function Call.
			projData[iter->first] = (float)projValue;														//double->float
		
			fwrite(projData, sizeof(float)*detectorWidth*detectorHeight, 1, fp);
			fclose(fp);
		}
	}
		
	// //detector_map의 projection Value값 출력 확인
	// cout <<endl ;
	// for(auto iter = detector_map.begin(); iter != detector_map.end(); ++iter){
	// 	int index = iter->first;
    //     if(index % detectorWidth ==0)		//Width 마다 개행
    //         cout<<endl;
	// 	//detector_map은 y,z축으로만 확장되있으므로 detector_map(key, y, z) 임.
	// 	cout << "(" << detector_X <<", "<< iter->second.first <<", " <<iter->second.second << ") : " <<   projData[index] <<" | ";
	// }
	

	free(projData);

	// 3차원 동적 배열 rho 메모리 해제
    for (int i = 0; i < matrix_x; ++i) {
        for (int j = 0; j < matrix_y; ++j) {
            delete[] rho[i][j];
        }
        delete[] rho[i];
    }
    delete[] rho;

	cout << endl << endl << "Main Function Exit success" << endl;

	return 0;
}