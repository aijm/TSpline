#include "VolumePiaMethod.h"

void VolumePiaMethod::calculate()
{
	parameterize();
	init();
	insert();
	//sample_fitPoints_2();
	sample_fitPoints_bvolume();
	fitPoints = surface_points;
	fitPoints.insert(fitPoints.end(), inter_points.begin(), inter_points.end());
	cal_basis_cache(); // 计算基函数在fitPoints参数位置的值， 这样在fit()和pia()时不需要再重新计算
	cout << "finished cal_basis_cache() " << endl;
	fit();
	cout << "finished fit() " << endl;
	pia();
	cout << "finished pia() " << endl;
	update();
	cout << "finished update() " << endl;

	for (int i = 0; i < 3; i++) {
		fit();
		cout << "finished fit() " << endl;
		pia();
		cout << "finished pia() " << endl;
		update();
		cout << "finished update() " << endl;
	}

	//for (int j = 0; j < 0; j++) {
	//	Point3d low, high;
	//	param_helper_points(low, high); // 辅助点参数化

	//	//fitPoints = curve_points;
	//	fitPoints.clear();
	//	fitPoints.insert(fitPoints.end(), helper_points.begin(), helper_points.end());
	//	for (auto point : inter_points) {
	//		if (!point.inRectangle(low, high)) {
	//			fitPoints.push_back(point);
	//		}
	//	}
	//	for (auto point : surface_points) {
	//		if (!point.inRectangle(low, high)) {
	//			fitPoints.push_back(point);
	//		}
	//	}
	//	fit();
	//	pia();
	//}
	//update();
}

void VolumePiaMethod::param_helper_points(Point3d & low, Point3d & high)
{

}


void VolumePiaMethod::sample_fitPoints_2()
{
	const int sampleNum =20;
	for (int i = 1; i < surfaces_num - 1; i++) {

		for (int j = 0; j <= sampleNum; j++) {
			for (int k = 0; k <= sampleNum; k++) {
				FitPoint3D point;
				point.param[0] = 1.0*j / sampleNum;
				point.param[1] = 1.0*k / sampleNum;
				point.param[2] = this->w_params(i);
				point.origin = surfaces[i].eval(point.param[0], point.param[1]);
				surface_points.push_back(point);
				
			}
		}

	}


	// 纵向采样，拟合出一个B样条曲线
	const int v_sample_num = 5;
	const int u_sample_num = 5;
	const int w_sample_num = 30;

	VectorXd params = w_params;
	params(0) = 0; params(params.size() - 1) = 1;
	cout << "params: \n" << params << endl;
	VectorXd knots(params.size() + 6);
	knots(0) = 0; knots(1) = 0, knots(2) = 0;
	knots(knots.size() - 1) = 1; knots(knots.size() - 2) = 1, knots(knots.size() - 3) = 1;
	knots.block(3, 0, params.size(), 1) = params;
	cout << "knots for interpolate: \n" << knots << endl;

	vector<vector<NURBSCurve>> sample_curves(u_sample_num + 1, vector<NURBSCurve>(v_sample_num + 1));
	
	for (int i = 0; i <= u_sample_num; i++) {
		for (int j = 0; j <= v_sample_num; j++) {
			double u = 1.0*i / u_sample_num;
			double v = 1.0*j / v_sample_num;
			MatrixXd points(surfaces_num, 3);
			MatrixXd tangent(surfaces_num, 3);
			for (int k = 0; k < surfaces_num; k++) {
				points.row(k) = surfaces[k].eval(u, v).toVectorXd();
				tangent.row(k) = surfaces[k].normal(u, v).toVectorXd();
				MatrixXd normal = MatrixXd::Zero(1, 3);
				// venus --> 0.01
				// tooth --> 5
				double factor = 5;
				normal.row(0) = tangent.row(k) * factor;
				MatrixXd start = MatrixXd::Zero(1, 3);
				start.row(0) = points.row(k);
				//viewer->data().add_points(start, red);
				//viewer->data().add_edges(start, start + normal, green);
			}
			
			sample_curves[i][j].interpolate(points, knots);
			//sample_curves[i][j].interpolate_optimize(points, tangent, params, 0.005);
			//sample_curves[i][j].interpolate_optimize1(points, tangent, params, 0.005);
			//sample_curves[i][j].interpolate_tangent(points, tangent, params);
			//sample_curves[i][j].interpolate_tangent_improve(points, tangent, params);
			sample_curves[i][j].draw(*viewer, false, true, 0.001);
			cout << "finished " <<  i << ", " << j << endl;
		}
	}

	for (int i = 0; i <= w_sample_num; i++) {
		bool valid = true;
		for (int k = 0; k < w_params.size(); k++) {
			if (abs(w_params(k) - 1.0*i / w_sample_num) < 0.01) {
				valid = false;
				break;
			}
		}
		if (!valid) {
			continue;
		}
		for (int jj = 0; jj <= u_sample_num; jj++) {
			for (int kk = 0; kk <= v_sample_num; kk++) {
				FitPoint3D point;
				point.param[0] = 1.0*jj / u_sample_num;
				point.param[1] = 1.0*kk / v_sample_num;
				point.param[2] = 1.0*i / w_sample_num;
				point.origin.fromVectorXd(sample_curves[jj][kk].eval(point.param[2]));
				/*cout << "param: " << point.param[1] << ", " << point.param[1] << ", " << point.param[2] << endl;
				cout << "sample origin: ";
				point.origin.output(cout);
				cout << endl;*/
				inter_points.push_back(point);
			}
		}	
	}
}

void VolumePiaMethod::sample_fitPoints_bvolume()
{
	const int sampleNum = 20;
	for (int i = 1; i < surfaces_num - 1; i++) {

		for (int j = 0; j <= sampleNum; j++) {
			for (int k = 0; k <= sampleNum; k++) {
				FitPoint3D point;
				point.param[0] = 1.0*j / sampleNum;
				point.param[1] = 1.0*k / sampleNum;
				point.param[2] = this->w_params(i);
				point.origin = surfaces[i].eval(point.param[0], point.param[1]);
				surface_points.push_back(point);

			}
		}

	}


	// 纵向采样，拟合出一个B样条曲线

	const int x_points = 21;
	const int y_points = 21;
	//const int z_points = surfaces_num + 2;
	const int z_points = 11;

	const int v_sample_num = 30;
	const int u_sample_num = 30;
	const int w_sample_num = 30;

	VectorXd params = w_params;
	params(0) = 0; params(params.size() - 1) = 1;
	cout << "params: \n" << params << endl;
	VectorXd knots(params.size() + 6);
	knots(0) = 0; knots(1) = 0, knots(2) = 0;
	knots(knots.size() - 1) = 1; knots(knots.size() - 2) = 1, knots(knots.size() - 3) = 1;
	knots.block(3, 0, params.size(), 1) = params;
	cout << "knots for interpolate: \n" << knots << endl;

	vector<vector<NURBSCurve>> sample_curves(x_points, vector<NURBSCurve>(y_points));

	//vector<FitPoint3D> fit_points;

	//for (int i = 0; i <= x_points - 1; i++) {
	//	for (int j = 0; j <= y_points - 1; j++) {
	//		double u = 1.0*i / (x_points - 1);
	//		double v = 1.0*j / (y_points - 1);
	//		MatrixXd points(surfaces_num, 3);
	//		for (int k = 0; k < surfaces_num; k++) {
	//			points.row(k) = surfaces[k].eval(u, v).toVectorXd();
	//			FitPoint3D fit_point;
	//			fit_point.origin.fromVectorXd(points.row(k).transpose());
	//			fit_point.param[0] = u;
	//			fit_point.param[1] = v;
	//			fit_point.param[2] = 1.0 * k / (surfaces_num - 1);
	//			fit_points.push_back(fit_point);
	//		}

	//		sample_curves[i][j].interpolate(points, knots);
	//		cout << "finished " << i << ", " << j << endl;
	//	}
	//}

	//BsplineVolume bvolume;
	//bvolume.control_grid = vector<vector<vector<Point3d>>>(x_points, vector<vector<Point3d>>(y_points, vector<Point3d>(z_points)));
	//bvolume.constructKnotVector(x_points, y_points, z_points);
	//for (int i = 0; i < x_points; i++) {
	//	for (int j = 0; j < y_points; j++) {
	//		for (int k = 0; k < z_points; k++) {
	//			bvolume.control_grid[i][j][k].fromVectorXd(sample_curves[i][j].controlPw.row(k).transpose());
	//		}
	//	}
	//}
	//// tooth --> alpha = 0.3 delta = 0.3   iter_num = 20
	//// venus --> alpha = 0.45 delta = 0.45  iter_num = 20
	//// head --> alpha = 0.495 delta = 0.495  iter_num = 1000
	//double alpha = 0.495;
	//double delta = 0.495;

	//// tooth --> iter_num = 20
	//for (int i = 0; i < 100; i++) {
	//	double error = bvolume.GetSoildFiterror(fit_points, x_points, y_points, z_points, alpha, delta);
	//	cout << "iter: " << i << ", error: " << error << endl;
	//	bvolume.fitBsplineSolid(fit_points, x_points, y_points, z_points, alpha, delta);
	//}
	//// moai --> true
	//bvolume.setReverse(true);
	//bvolume.saveAsHex("../out/volume/tooth_fitbspline", 0.01);


	vector<FitPoint3D> fit_points;

	for (int i = 0; i <= x_points - 1; i++) {
		for (int j = 0; j <= y_points - 1; j++) {
			double u = 1.0*i / (x_points - 1);
			double v = 1.0*j / (y_points - 1);
			MatrixXd points(surfaces_num, 3);
			for (int k = 0; k < surfaces_num; k++) {
				points.row(k) = surfaces[k].eval(u, v).toVectorXd();
				/*FitPoint3D fit_point;
				fit_point.origin.fromVectorXd(points.row(k).transpose());
				fit_point.param[0] = u;
				fit_point.param[1] = v;
				fit_point.param[2] = 1.0 * k / sample_num;
				fit_points.push_back(fit_point);*/
			}

			sample_curves[i][j].interpolate(points, knots);
			cout << "finished " << i << ", " << j << endl;
		}
	}
	const int u_num = 20;
	const int v_num = 20;
	for (int i = 0; i < surfaces_num; i++) {
		for (int j = 0; j <= u_num; j++) {
			for (int k = 0; k <= v_num; k++) {
				FitPoint3D fit_point;
				double u = 1.0 * j / u_num;
				double v = 1.0 * k / v_num;
				fit_point.param[0] = u;
				fit_point.param[1] = v;
				fit_point.param[2] = 1.0 * i / (surfaces_num - 1);
				fit_point.origin = surfaces[i].eval(u, v);
				fit_points.push_back(fit_point);
			}
		}
	}


	BsplineVolume bvolume;
	bvolume.control_grid = vector<vector<vector<Point3d>>>(x_points, vector<vector<Point3d>>(y_points, vector<Point3d>(z_points)));
	bvolume.constructKnotVector(x_points, y_points, z_points);
	for (int i = 0; i < x_points; i++) {
		for (int j = 0; j < y_points; j++) {
			for (int k = 0; k < z_points; k++) {
				//bvolume.control_grid[i][j][k].fromVectorXd(sample_curves[i][j].controlPw.row(k).transpose());
				double w = 1.0 * k / (z_points - 1);
				bvolume.control_grid[i][j][k].fromVectorXd(sample_curves[i][j].eval(w));
			}
		}
	}


	// using lspia
	bvolume.lspia(fit_points, x_points, y_points, z_points, 10, 1e-12);


	double alpha = 0.495;
	double delta = 0.495;
	
	for (int i = 0; i < 20; i++) {
		double error = bvolume.GetSoildFiterror(fit_points, x_points, y_points, z_points, alpha, delta);
		cout << "iter: " << i << ", error: " << error << endl;
		bvolume.fitBsplineSolid(fit_points, x_points, y_points, z_points, alpha, delta);
	}
	bvolume.setReverse(true);
	string modelname = "moai_fitbspline";
	//bvolume.saveVolume("../out/volume/" + modelname + "_fitbspline");
	bvolume.saveAsHex("../out/volume/" + modelname + "_fitbspline", 0.01);





	for (int i = 0; i <= w_sample_num; i++) {
		/*bool valid = true;
		for (int k = 0; k < w_params.size(); k++) {
			if (abs(w_params(k) - 1.0*i / w_sample_num) < 0.01) {
				valid = false;
				break;
			}
		}
		if (!valid) {
			continue;
		}*/
		for (int jj = 0; jj <= u_sample_num; jj++) {
			for (int kk = 0; kk <= v_sample_num; kk++) {
				/*cout << "-----------------------" << endl;*/
				FitPoint3D point;
				point.param[0] = 1.0*jj / u_sample_num;
				point.param[1] = 1.0*kk / v_sample_num;
				point.param[2] = 1.0*i / w_sample_num;
				point.origin = bvolume.eval(point.param[0], point.param[1], point.param[2]);
				/*cout << "param: " << point.param[1] << ", " << point.param[1] << ", " << point.param[2] << endl;
				cout << "sample origin: ";
				point.origin.output(cout);
				cout << endl;*/
				inter_points.push_back(point);
			}
		}
	}
}

void VolumePiaMethod::sample_fitPoints()
{
	const int sampleNum = 10;
	for (int i = 1; i < surfaces_num - 1; i++) {

		for (int j = 0; j <= sampleNum; j++) {
			for (int k = 0; k <= sampleNum; k++) {
				FitPoint3D point;
				point.param[0] = 1.0*j / sampleNum;
				point.param[1] = 1.0*k / sampleNum;
				point.param[2] = this->w_params(i);
				point.origin = surfaces[i].eval(point.param[0], point.param[1]);
				surface_points.push_back(point);
				MatrixXd P;
				array2matrixd(point.origin, P);
				(*viewer).data().add_points(P, red);
			}
		}

	}

	// calculate sample points by linear interpolate
	for (int i = 0; i <= surfaces_num - 2; i++) {
		double w_now = w_params(i);
		double w_next = w_params(i + 1);
		double w_inter1 = 2.0 / 3 * w_now + 1.0 / 3 * w_next;
		double w_inter2 = 1.0 / 3 * w_now + 2.0 / 3 * w_next;

		for (int j = 0; j <= sampleNum; j++) {
			for (int k = 0; k <= sampleNum; k++) {
				FitPoint3D point1, point2;
				point1.param[0] = 1.0*j / sampleNum;
				point1.param[1] = 1.0*k / sampleNum;
				point2.param[0] = point1.param[0];
				point2.param[1] = point1.param[1];
				point1.param[2] = w_inter1;
				point2.param[2] = w_inter2;

				Point3d now_coor = surfaces[i].eval(point1.param[0], point1.param[1]);
				Point3d next_coor = surfaces[i + 1].eval(point1.param[0], point1.param[1]);
				point1.origin = 2.0 / 3 * now_coor + 1.0 / 3 * next_coor;
				point2.origin = 2.0 / 3 * next_coor + 1.0 / 3 * now_coor;

				inter_points.push_back(point1);
				inter_points.push_back(point2);

				MatrixXd P;
				array2matrixd(point1.origin, P);
				(*viewer).data().add_points(P, green);
				array2matrixd(point2.origin, P);
				(*viewer).data().add_points(P, green);
			}
		}
	}
	cout << "number of points: " << surface_points.size() + inter_points.size() << endl;
}

void VolumePiaMethod::fit()
{
	error = 0.0;
	for (int i = 0; i < fitPoints.size();i++) {
		auto& point = fitPoints[i];

		Point3d val;
		int count = 0;
		for (auto entry : volume.w_map) {
			for (auto node : entry.second->nodes) {
				val.add(node->data * basis_cache[count][i]);
				count++;
			}
			
		}
		point.eval = val;
	
		point.error = point.geterror();
		error += point.error;
	}
	error /= fitPoints.size();
}

void VolumePiaMethod::pia()
{
	for (int i = 0; i < maxIterNum; i++) {
		// 计算差向量并更新曲面控制点

		int count = 0;
		for (auto &entry : volume.w_map) {
			// 一层
			if (entry.first <= 0.0000 || entry.first >= 1.0) {
				count += entry.second->get_num();
				continue;
			}

			for (auto node : entry.second->nodes) {
				/*if (node->s[2] <= 0.0001 || node->s[2] >= 0.9999) {
					continue;
				}*/
				
				Point3d sum2;
				for (int j = 0; j < fitPoints.size();j++) {
					double blend = basis_cache[count][j];
					
					Point3d delta = fitPoints[j].origin - fitPoints[j].eval;
					delta.scale(blend);
					sum2.add(delta);
				}
				double sum1 = basis_cache_sum[count];
				double factor = 0.0;
				if (abs(sum1) > 0.0001) {
					factor = 1.0 / sum1;
				}
				sum2.scale(factor); // 差向量
				node->data.add(sum2); // 更新坐标	

				count++;
			}
			
		}
		

		fit();
		cout << "iter: " << i + 1 << ", error: " << error << endl;
		if (error < eps) {
			break;
		}

	}
}

void VolumePiaMethod::cal_basis_cache()
{
	basis_cache.resize(volume.get_num());
	basis_cache_sum.resize(basis_cache.size(), 0);
	// 计算每个节点处的基函数在要拟合点参数处的值 B_i(t_j)
	int layer = 0; // 用于计算w向基函数
	int count = 0;
	for (auto entry : volume.w_map) {
		for (auto node : entry.second->nodes) {
			basis_cache[count].resize(fitPoints.size());
			for (int i = 0; i < fitPoints.size(); i++) {
				auto point = fitPoints[i];
				double blend = node->basis(point.param[0], point.param[1]) * Basis(volume.w_knots, point.param[2], layer);
				basis_cache[count][i] = blend;
				basis_cache_sum[count] += blend;
			}
			count++;
		}
		layer++;
	}
}
