#pragma once
#ifndef BSPLINEVOLUME_H
#define BSPLINEVOLUME_H
#include "Volume.h"
#include "FitPoint.hpp"

const double PI = 3.141592653;

class BsplineVolume :public Volume{
public:
	Point3d eval(double u, double v, double w) override;

	int readVolume(string) override;
	int saveVolume(string) override;

	void get_isoparam_surface(NURBSSurface& surface, double t, char dir);

	// 以追加的方式打开
	static void writeFileAdd(char *filename, std::string s);

	// 点乘
	static double dot(Point3d v1, Point3d v2);

	// 叉乘
	static Point3d cross(Point3d v1, Point3d v2);

	//计算出两个圆锥之间的角度
	//@ 负值的绝对值越来越大的话表示两个圆锥越来越远离 ，负值的绝对值越来越小的话 代表section_angle-两个圆心角越来越小，从而两个圆锥越来越靠近
	//@ 总体而言就是越来越小圆锥越来越远离
	static double degreeofTwoCone(std::vector<Point3d> &cone1, std::vector<Point3d> &cone2);
	/*
	*	判断圆锥是不是相交的
	*  @返回值为true代表的是相交 ,false 代表的是不相交
	*/
	static bool coneIntersect(std::vector<Point3d> &cone1, std::vector<Point3d> &cone2);

	//利用另外一个计算圆锥的方法来计算圆锥的角度
	static double degreeofTwoConeWithAnotherFunc(std::vector<Point3d> &cone1, std::vector<Point3d> &cone2);

	//得到包含区域的最小的圆锥
	static void getconescale(Point3d &v1, Point3d &v2, std::vector<Point3d> &cone1);

	//矩阵乘法
	static Point3d matmul(double mat[][3], Point3d v);

	//二维空间上三角形的外心
	static Point3d getCircumcenter(Point3d v1, Point3d v2, Point3d v3);

	//利用投影变换计算出三维空间上的外心位置
	//需要保证v1,v2, v3 不为空间上的共线的点
	static Point3d getCircumcenterOnThreeDim(Point3d v1, Point3d v2, Point3d v3);

	//判断一个顶点和一个
	static bool isContainAllLine(Point3d v1, Point3d v2, std::vector<Point3d> &cone1);

	//得到包含区域的最小的圆锥
	static void getconescaleWithAnotherFunc(Point3d &v1, Point3d &v2, std::vector<Point3d> &cone1);


	//////////////////////////////////////////////////////////////////////////
	//判断圆锥两两之间之间不相交
	//相交返回true
	//不相交则返回false
	static bool coneIntersect3D(std::vector<Point3d> &conex, std::vector<Point3d> &coney, std::vector<Point3d> &conez);


	//通过control_grid来构建三个圆锥
	//cone 为三个方向上线段构成的圆锥的集合
	void constructsolidcone(int seg[], std::vector<Point3d> &conex, vector<Point3d> &coney,
		vector<Point3d> &conez);

	double FindBestratioOnSoildBysearch(const vector<vector<vector<Point3d>>>& diff_vector, vector<Point3d> &conex, vector<Point3d> &coney,
		vector<Point3d> &conez, int x_points, int y_points, int z_points);


	//求点(x ,y ,z)处Ep3的值
	double GetApartFuncvalueOnSoild(int x_points, int y_points, int z_points);

	Point3d getDiffofApartFuncOnSoild(int x_points, int y_points, int z_points, double delta, int x, int y, int z);

	void constructKnotVector(int x_points, int y_points, int z_points);

	//计算一个参数点返回Bi(u)Bj(v)Bk(w)
	//para为计算点的参数
	//kont_vector节点向量，大小为3，表示xyz三个轴
	//result表示4*4*4的结果矩阵，每一维表示Bi*Bj*Bk
	//Bi_start_index表示返回值Bi的向量的基函数的起始节点向量的下标
	void calculateBaseFunctionOfSolidBSpline(Point3d para,
		vector<vector<vector<double>>>& result, int &Bi_start_index, int &Bj_start_index, int &Bk_start_index);

	// 使用优化方法，生成一个拟合指定点，且雅克比值很好的B样条体
	void fitBsplineSolid(vector<FitPoint3D>& fit_points, int x_points, int y_points, int z_points, double alpha, double delta);

	//体拟合迭代的误差
	double GetSoildFiterror(vector<FitPoint3D>& fit_points,
		int x_points, int y_points, int z_points, double alpha, double delta);

	static double getConeAnglerror(std::vector<Point3d> &cones);

	static double getConeAngleandAngle(std::vector<Point3d> &conex, std::vector<Point3d> &coney);

	void lspia(vector<FitPoint3D>& fit_points, int x_points, int y_points, int z_points, int max_iter_num = 100, double eps = 1e-5);
private:
	void drawTmesh() override;
	void drawControlpolygon() override;

public:
	vector<vector<vector<Point3d>>> control_grid;
	vector<Eigen::VectorXd> knot_vector;
	vector<vector<vector<vector<double>>>> matri;
	vector<int> Bi_start_indexs;
	vector<int> Bj_start_indexs;
	vector<int> Bk_start_indexs;
};
#endif // !BSPLINEVOLUME_H

