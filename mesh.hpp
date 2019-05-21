#include"utility.h"
#include<cassert>
#include<functional>
#include<algorithm>
#include<map>
#include<set>
#include<list>
#include<cmath>
#include<string>
#include<sstream>
#include<fstream>
#include"draw.h"
#include"NURBSCurve.h"
#include"NURBSSurface.h"
//#include<boost/thread/thread.hpp>
//#include<boost/bind.hpp>
//#include<boost/progress.hpp>
//#include<boost/timer.hpp>
//#include<boost/unordered_map.hpp>

namespace t_mesh{
    using namespace std;
	using namespace igl::opengl::glfw;
    template<class T>
    class Mesh{
        public:
            int loadMesh(string);
            int saveMesh(string);
			T eval(double s,double t);
			void drawTmesh(Viewer &viewer);
			void drawControlpolygon(Viewer &viewer);
			void drawSurface(Viewer &view, double resolution = 0.01);
			void draw(Viewer &viewer, bool tmesh, bool polygon, bool surface,double resolution = 0.01);
			int get_num() const { return nodes.size(); }
			VectorXd skinning_parameterize(const vector<NURBSCurve> &curves);
			void skinning_init(const VectorXd &s_knots, const vector<NURBSCurve> &curves);
			void skinning_insert(const VectorXd &s_knots, const vector<NURBSCurve> &curves);
			void skinning_insert1(const VectorXd &s_knots, const vector<NURBSCurve> &curves);
			void skinning_intermediate(const VectorXd &s_knots, const vector<NURBSCurve> &curves,Viewer &viewer);
			void intermediate_init(const VectorXd &s_knots, const vector<NURBSCurve>& curves, map<double,map<double,T>> &initial_cpts,Viewer &viewer, int n_sample = 100);
			double intermediate_update(const VectorXd &s_knots, const vector<NURBSCurve>& curves, map<double, map<double, T>> &initial_cpts, Viewer &viewer, int n_sample = 100, bool showiteration = false);
			void skinning_update_cross(const VectorXd &s_knots, const vector<NURBSCurve> &curves);
			void skinning(const vector<NURBSCurve> &curves,Viewer &viewer);
			
			void pia_fit(const VectorXd &s_params, const vector<VectorXd>& ttparams, const vector<NURBSCurve> &curves, Viewer &viewer, int max_iterations = 100, double eps = 1e-5);
			void pia_skinning(const vector<NURBSCurve> &curves, Viewer &viewer, int max_iterations = 100, double eps = 1e-5);
			double maxoffset(double s, int n, const NURBSCurve& curve, double& max_off);

          
            void insert(double s,double t);
            

            Node<T>*    new_node();
            Node<T>*    get_node(int num);
            Node<T>*    get_node(double s,double t);
            Node<T>     get_knot(double x,double y); // get the knot vector of (x,y)

			Eigen::MatrixXd mesh_V;
			Eigen::MatrixXi mesh_F;
        private:
			void basis_split(const map<double, Node<T>*> &fewer_map, const map<double, Node<T>*> &more_map, map<double,T> &coeff);
            int insert_helper(double s,double t,bool changedata=true);
			void adjust(Node<T>* n, bool changedata = true);
            void merge_all();
            bool check_valid();
			void clear();

            vector<Node<T>*>    nodes;
            
            // organizing node in a good data structure 
            map<double,map<double,Node<T>*> > s_map; // s_map[s][t]
            map<double,map<double,Node<T>*> > t_map; // t_map[t][s]
            list<Node<T> >              pool;

			double width = 1.0;
			double height = 1.0;
           
            string          iter_str;
			
    };

	template<class T>
	inline void t_mesh::Mesh<T>::clear()
	{
		this->pool.clear();
		this->nodes.clear();
		this->s_map.clear();
		this->t_map.clear();
	}
	
	 template<class T>
		 T Mesh<T>::eval(double s, double t) {
			 T result;
			 for (int i = 0; i < nodes.size(); i++) {
				 if (nodes[i]->is_ok(s, t)) {
					 double blend = Basis((nodes[i]->s).toVectorXd(), s)*Basis((nodes[i]->t).toVectorXd(), t);
					 T temp(nodes[i]->data);
					 //temp.output(cout);
					 temp.scale(blend);
					 result.add(temp);
				 }
			 }

			 /*for (auto it = s_map.begin(); it != s_map.end(); ++it) {
				 for (auto it1 = (it->second).begin(); it1 != (it->second).end(); ++it1) {
					 if (it1->second->t[4] < t)
						 continue;
					 if (it1->second->t[0] > t)
						 break;
					 if (it1->second->is_ok(s, t)) {
						 double blend = Basis((it1->second->s).toVectorXd(), s)*Basis((it1->second->t).toVectorXd(), t);
						 T temp = it1->second->data;
						 temp.scale(blend);
						 result.add(temp);
					 }
				 }
			 }*/
			 return result;
		 }

	 template<class T>
		 void Mesh<T>::drawTmesh(Viewer &viewer){
			 Eigen::MatrixXd P1(1,2), P2(1,2);
			 Eigen::MatrixXd nodes_st(nodes.size(), 2);
			 typedef typename map<double, map<double, Node<T>*> >::iterator   map_t;
			 typedef typename map<double, Node<T>*>::iterator             map2_t;
			 viewer.data().add_label(Vector3d(0, 0, 0), "haha");
			 //cout << "1" << endl;
			 for (int i = 0; i < nodes.size(); i++) {
				 nodes_st.row(i) << nodes[i]->s[2], nodes[i]->t[2];
				 std::stringstream label;
				 label << nodes[i]->get_order()<<":"<<nodes_st(i, 0) << ", " << nodes_st(i, 1);
				
				 viewer.data().add_label(nodes_st.row(i), label.str());
				 viewer.data().add_points(nodes_st.row(i), red);
				 //cout << "position: " << nodes_st.row(i) << ", label: " << label.str() << endl;
			 }
			 //cout << "2" << endl;
			 for (auto iter = s_map.begin(); iter != s_map.end(); ++iter) {
				 for (auto iter1 = (iter->second).begin(); iter1 != (iter->second).end(); ++iter1) {
					 if ((iter1->second)->adj[2]) {
						 P1(0, 0) = iter->first; 
						 P1(0, 1) = iter1->first;
						 P2(0, 0) = iter->first;
						 P2(0, 1) = (iter1->second)->adj[2]->t[2];
						 viewer.data().add_edges(P1, P2, white);
					 }
				 }
			 }

			 for (auto iter = t_map.begin(); iter != t_map.end(); ++iter) {
				 for (auto iter1 = (iter->second).begin(); iter1 != (iter->second).end(); ++iter1) {
					 if ((iter1->second)->adj[1]) {
						 P1(0, 0) = iter1->first;
						 P1(0, 1) = iter->first;
						 P2(0, 0) = (iter1->second)->adj[1]->s[2];
						 P2(0, 1) = iter->first;
						 viewer.data().add_edges(P1, P2,white);
					 }
				 }
			 }
			 viewer.core.align_camera_center(nodes_st); // center
		 }

	  template<class T>
	      void Mesh<T>::drawControlpolygon(Viewer &viewer) {
			  Eigen::MatrixXd P1, P2;
			  array2matrixd(nodes[0]->data, P1);
			  Eigen::MatrixXd nodes_point(nodes.size(), P1.cols());
			  
			  for (int i = 0; i < nodes.size(); i++) {
				  array2matrixd(nodes[i]->data, P1);
				  nodes_point.row(i) = P1;
			  }
			  
			  for (auto iter = s_map.begin(); iter != s_map.end(); ++iter) {
				  for (auto iter1 = (iter->second).begin(); iter1 != (iter->second).end(); ++iter1) {
					  if ((iter1->second)->adj[2]) { 
						  array2matrixd(iter1->second->data, P1);
						  array2matrixd((iter1->second->adj[2])->data, P2);
						  viewer.data().add_edges(P1, P2, green);
					  }
				  }
			  }

			  for (auto iter = t_map.begin(); iter != t_map.end(); ++iter) {
				  for (auto iter1 = (iter->second).begin(); iter1 != (iter->second).end(); ++iter1) {
					  if ((iter1->second)->adj[1]) {
						  array2matrixd(iter1->second->data, P1);
						  array2matrixd((iter1->second->adj[1])->data, P2);
						  viewer.data().add_edges(P1, P2, blue);
					  }
				  }
			  }
			  viewer.data().add_points(nodes_point, red);
			  viewer.core.align_camera_center(nodes_point); // center
			  
		  }

	  template<class T>
	      void Mesh<T>::drawSurface(Viewer &viewer, double resolution) {
			  // cut apart the parameter domain
			  double u_low = (++s_map.begin())->first;
			  
			  double u_high = (++s_map.rbegin())->first;
			  cout << "u_low: " << u_low <<", u_high: "<<u_high<< endl;
			  const int uspan = (u_high - u_low) / resolution;
			  double u_resolution = (u_high - u_low) / uspan;
			  
			  double v_low = (++t_map.begin())->first;
			  //v_low = 1.0;
			  double v_high = (++t_map.rbegin())->first;
			  //v_high = 2.0;
			  const int vspan = (v_high - v_low) / resolution;
			  double v_resolution = (v_high - v_low) / vspan;
			  cout << "v_low: " << v_low << ", v_high: " << v_high << endl;
			  mesh_V = Eigen::MatrixXd((uspan + 1)*(vspan + 1), 3);
			  mesh_F = Eigen::MatrixXi(2 * uspan*vspan, 3);
			  // discretize T-Spline Surface into triangular mesh(V,F) in libigl mesh structure
			  // calculate 
			  
			  for (int j = 0; j <= vspan; j++)
				  for (int i = 0; i <= uspan; i++){
					  Eigen::MatrixXd curvePoint;
					  double u = u_low + i*u_resolution;
					  double v = v_low + j*v_resolution;
					  array2matrixd(eval(u_low + i*u_resolution, v_low + j*v_resolution), curvePoint);
					  //cout << "u, v, point: \n" <<u<<" "<<v<<" "<< curvePoint << endl;
					  mesh_V.row(j*(uspan + 1) + i) = curvePoint;
				  }

			  for (int j = 0; j<vspan; j++)
				  for (int i = 0; i < uspan; i++){
					  int V_index = j*(uspan + 1) + i;
					  int F_index = 2 * j*uspan + 2 * i;
					  mesh_F.row(F_index) << V_index, V_index + 1, V_index + uspan + 1;
					  mesh_F.row(F_index + 1) << V_index + uspan + 1, V_index + 1, V_index + uspan + 2;
				  }
			  
			  viewer.data().set_mesh(mesh_V, mesh_F);
		  }

		  template<class T>
		  void Mesh<T>::draw(Viewer & viewer,bool tmesh, bool polygon, bool surface, double resolution){
			  if (tmesh) {
				  drawTmesh(viewer);
				  return;
			  }
			  if (polygon) {
				  drawControlpolygon(viewer);
			  }
			  if (surface) {
				  drawSurface(viewer, resolution);
			  }
		  }

		  template<class T>
		  inline VectorXd Mesh<T>::skinning_parameterize(const vector<NURBSCurve>& curves)
		  {
			  // 1. compute s-knot for curves
			  assert(curves.size() > 0);
			  const int curves_num = curves.size();
			  const int dimension = curves[0].controlPw.cols();
			  
			  VectorXd s_knots(curves_num);
			  MatrixXd u_cpts(curves_num, dimension);
			  for (int i = 0; i < u_cpts.rows(); i++) {
				  u_cpts.row(i) = curves[i].controlPw.row(0);
			  }
			  s_knots = NURBSCurve::parameterize(u_cpts);
			  s_knots(0) = 0.0001; s_knots(s_knots.size() - 1) = 0.9999;

			  cout << "s_knots: " << s_knots.transpose() << endl;
			  return s_knots;
		  }

		  template<class T>
		  inline void Mesh<T>::skinning_init(const VectorXd & s_knots, const vector<NURBSCurve>& curves)
		  {
			  // 2. construct basis T-mesh 
			  //add 0 and 1
			  const int curves_num = curves.size();
			  for (int j = 0; j <= curves[0].n; j++) {
				  double vknot = curves[0].knots(j + 2);
				  if (j == 1) { vknot = 0.0001; }
				  if (j == curves[0].n-1) { vknot = 0.9999; }
				  
				  insert_helper(0.0, vknot, false);
				  Node<T>* node = get_node(0.0, vknot);
				  (node->data).fromVectorXd(curves[0].controlPw.row(j));
				  /*(node->data).output(cout);
				  cout << endl;*/
			  }

			  for (int i = 0; i < curves_num; i++) {
				  for (int j = 0; j <= curves[i].n; j++) {
					  double vknot = curves[i].knots(j + 2);
					  if (j == 1) { vknot = 0.0001; }
					  if (j == curves[i].n - 1) { vknot = 0.9999; }
					  insert_helper(s_knots(i), vknot, false);
					  Node<T>* node = get_node(s_knots(i), vknot);
					  (node->data).fromVectorXd(curves[i].controlPw.row(j));
					  /*(node->data).output(cout);
					  cout << endl;*/
					  //merge_all();
				  }
			  }
			  for (int j = 0; j <= curves[curves_num - 1].n; j++) {
				  double vknot = curves[curves_num-1].knots(j + 2);
				  if (j == 1) { vknot = 0.0001; }
				  if (j == curves[curves_num-1].n - 1) { vknot = 0.9999; }
				  insert_helper(1.0, vknot, false);
				  Node<T>* node = get_node(1.0, vknot);
				  (node->data).fromVectorXd(curves[curves_num - 1].controlPw.row(j));
				  /*(node->data).output(cout);
				  cout << endl;*/
			  }

			  cout << "pool size:" << pool.size() << endl;
			  pool.clear();

			  if (!check_valid()) {
				  cout << "skinning: invalid T-mesh!" << endl;
				  return;
			  }
		  }

		  template<class T>
		  inline void Mesh<T>::skinning_insert(const VectorXd & s_knots, const vector<NURBSCurve>& curves)
		  {
			  const int curves_num = curves.size();
			  // 3. insert intermediate vertices
			  // the coordinate of vertices is the midpoint of the corresponding points in C_r and C_(r+1)
			  assert(curves_num >= 3);
			  for (int i = 0; i <= curves_num - 2; i++) {
				  double s_now = s_knots(i);
				  auto s_nodes = s_map[s_now];
				  for (auto it = s_nodes.begin(); it != s_nodes.end(); ++it) {
					  if (it->second->adj[1]) {
						  double s_mid = (s_now + s_knots(i + 1)) / 2;
						  insert_helper(s_mid, it->first, false);
						  Node<T>* node = get_node(s_mid, it->first);
						  //(node->data).add(it->second->data);
						  //(node->data).add(it->second->adj[1]->data);
						  // bug, it->second->adj[1] has changed to current node
						  (node->data).add(node->adj[1]->data);
						  (node->data).add(node->adj[3]->data);
						  (node->data).scale(0.5);
					  }
				  }
			  }
			  pool.clear();
			  if (!check_valid()) {
				  cout << "skinning: invalid T-mesh!" << endl;
				  return;
			  }
		  }

		  template<class T>
		  inline void Mesh<T>::skinning_insert1(const VectorXd & s_knots, const vector<NURBSCurve>& curves)
		  {
			  const int curves_num = curves.size();
			  // 3. insert intermediate vertices
			  // the coordinate of vertices is the midpoint of the corresponding points in C_r and C_(r+1)
			  assert(curves_num >= 3);
			  for (int i = 0; i <= curves_num-1; i++) {
				  double s_now = s_knots(i);
				  auto s_nodes = s_map[s_now];
				  if (i == 0) {
					  for (auto it = s_nodes.begin(); it != s_nodes.end(); ++it) {
						  double s_insert = 2.0 / 3 * s_now + 1.0 / 3 * s_knots(i + 1);
						  insert_helper(s_insert, it->first, false);
						  MatrixXd position = 2.0 / 3 * curves[i].eval(it->first) + 1.0 / 3 * curves[i + 1].eval(it->first);
						  (s_map[s_insert][it->first]->data).fromVectorXd(position.row(0).transpose());

					  }
				  }
				  else if (i == curves_num - 1) {
					  for (auto it = s_nodes.begin(); it != s_nodes.end(); ++it) {
						  double s_insert = 2.0 / 3 * s_now + 1.0 / 3 * s_knots(i - 1);
						  insert_helper(s_insert, it->first, false);
						  MatrixXd position = 2.0 / 3 * curves[i].eval(it->first) + 1.0 / 3 * curves[i - 1].eval(it->first);
						  (s_map[s_insert][it->first]->data).fromVectorXd(position.row(0).transpose());
					  }
				  }
				  else {
					  for (auto it = s_nodes.begin(); it != s_nodes.end(); ++it) {
						  double s_left = 2.0 / 3 * s_now + 1.0 / 3 * s_knots(i - 1);
						  double s_right = 2.0 / 3 * s_now + 1.0 / 3 * s_knots(i + 1);
						  insert_helper(s_left, it->first, false);
						  MatrixXd position = 2.0 / 3 * curves[i].eval(it->first) + 1.0 / 3 * curves[i - 1].eval(it->first);
						  (s_map[s_left][it->first]->data).fromVectorXd(position.row(0).transpose());

						  insert_helper(s_right, it->first, false);
						  position = 2.0 / 3 * curves[i].eval(it->first) + 1.0 / 3 * curves[i + 1].eval(it->first);
						  (s_map[s_right][it->first]->data).fromVectorXd(position.row(0).transpose());
					  }
				  }
				  
			  }
			  pool.clear();
			  if (!check_valid()) {
				  cout << "skinning: invalid T-mesh!" << endl;
				  return;
			  }
		  }

		  

		  template<class T>
		  inline void Mesh<T>::intermediate_init(
			  const VectorXd &s_knots, 
			  const vector<NURBSCurve>& curves, 
			  map<double, map<double,T>>& initial_cpts,
			  Viewer &viewer,
			  int n_sample)
		  {
			  const int curves_num = curves.size();
			  const int dimension = curves[0].controlPw.cols();
			  assert(curves_num >= 2);
			  for (int i = 0; i <= curves_num - 2; i++) {
				  double s_now = s_knots(i);
				  double s_next = s_knots(i + 1);
				  auto node_now = (s_map[s_now].begin())->second;
				  auto node_next = (s_map[s_next].begin())->second;
				  double s_inter1 = node_now->adj[1]->s[2];
				  double s_inter2 = node_next->adj[3]->s[2];
				  MatrixXd sample_inter1(n_sample + 1, dimension);
				  MatrixXd sample_inter2(n_sample + 1, dimension);
				  // sample points and linear interpolate
				  for (int j = 0; j <= n_sample; j++) {
					  RowVectorXd now_coor = curves[i].eval(1.0*j / n_sample);
					  RowVectorXd next_coor = curves[i + 1].eval(1.0*j / n_sample);
					  sample_inter1.row(j) = 2.0 / 3 * now_coor + 1.0 / 3 * next_coor;
					  sample_inter2.row(j) = 2.0 / 3 * next_coor + 1.0 / 3 * now_coor;
				  }
				  /*if (true) {
					  viewer.data().add_points(sample_inter1, blue);
					  viewer.data().add_points(sample_inter2, blue);
				  }*/
				  
				  /*if (i == 2) {
					  savepoints("sample_inter1" + std::to_string(i) + ".p", sample_inter1);
					  savepoints("sample_inter2" + std::to_string(i) + ".p", sample_inter2);
				  }*/
				  
				  // fit the sample points by LSPIA method with appointed knot vector
				  NURBSCurve inter1, inter2;
				  inter1.lspiafit(sample_inter1, curves[i].n + 1, curves[i].knots, 1000);
				  
				  inter2.lspiafit(sample_inter2, curves[i + 1].n + 1, curves[i + 1].knots, 1000);
				  /*if (true) {
					  inter1.draw(viewer, true, true, 0.001);
					  inter2.draw(viewer, true, true, 0.001);
					  cout << "inter1.knots: " << inter1.knots.transpose() << endl;
					  cout << "inter2.knots: " << inter2.knots.transpose() << endl;
				  }*/
				  

				  inter1.knots(3) = 0.0001; inter1.knots(inter1.n + 1) = 0.9999;
				  inter2.knots(3) = 0.0001; inter2.knots(inter2.n + 1) = 0.9999;

				  // update coordinate of intermediate control points and save as inital_cpts
				  for (int j = 0; j <= inter1.n; j++) {
					  T temp;
					  
					  
					  temp.fromVectorXd(inter1.controlPw.row(j));
					 
					  initial_cpts[s_inter1][inter1.knots(j + 2)] = temp;
					  s_map[s_inter1][inter1.knots(j + 2)]->data = temp;
					  
				  }
				  for (int j = 0; j <= inter2.n; j++) {
					  T temp;
					  temp.fromVectorXd(inter2.controlPw.row(j));

					  initial_cpts[s_inter2][inter2.knots(j + 2)] = temp;
					  s_map[s_inter2][inter2.knots(j + 2)]->data = temp;
				  }
			  }
			  cout << "finished intermediate_init()!" << endl;

			  
		  }

		  template<class T>
		  inline double Mesh<T>::intermediate_update(
			  const VectorXd & s_knots, 
			  const vector<NURBSCurve>& curves,
			  map<double, map<double, T>>& initial_cpts,
			  Viewer & viewer, int n_sample,bool showiteration)
		  {
			  if (showiteration) {
				  draw(viewer, false, false, true);
			  }
			  
			  const int curves_num = curves.size();
			  const int dimension = curves[0].controlPw.cols();
			  double error = 0.0;
			  assert(curves_num >= 2);
			  map<double, map<double, T>> T_cpts;

			  for (int i = 0; i <= curves_num - 2; i++) {
				  double s_now = s_knots(i);
				  double s_next = s_knots(i + 1);
				  auto node_now = (s_map[s_now].begin())->second;
				  auto node_next = (s_map[s_next].begin())->second;
				  double s_inter1 = node_now->adj[1]->s[2];
				  double s_inter2 = node_next->adj[3]->s[2];
				  MatrixXd T_inter1(n_sample + 1, dimension);
				  MatrixXd T_inter2(n_sample + 1, dimension);
				  // points on T-spline surface
				  for (int j = 0; j <= n_sample; j++) {
					  if (j == 0) {
						  T_inter1.row(j) = (node_now->adj[1]->data).toVectorXd();
						  T_inter2.row(j) = (node_next->adj[3]->data).toVectorXd();
						  
					  }
					  else if (j == n_sample) {
						  T_inter1.row(j) = ((s_map[s_now].rbegin())->second->adj[1]->data).toVectorXd();
						  T_inter2.row(j) = ((s_map[s_next].rbegin())->second->adj[3]->data).toVectorXd();
					  }
					  else {
						  T_inter1.row(j) = eval(s_inter1, 1.0*j / n_sample).toVectorXd();
						  T_inter2.row(j) = eval(s_inter2, 1.0*j / n_sample).toVectorXd();
					  }
					  
				  }
				  if (showiteration) {
					  viewer.data().add_points(T_inter1, green);
					  viewer.data().add_points(T_inter2, green);
				  }
				  

				  // fit the sample points by LSPIA method with appointed knot vector
				  NURBSCurve inter1, inter2;
				  inter1.lspiafit(T_inter1, curves[i].n + 1, curves[i].knots,1000);

				  inter2.lspiafit(T_inter2, curves[i + 1].n + 1, curves[i + 1].knots,1000);
				  if (showiteration) {
					  inter1.draw(viewer, true, true);
					  inter2.draw(viewer, true, true);
				  }
				  /*cout << "inter1.knots: " << inter1.knots.transpose() << endl;
				  cout << "inter2.knots: " << inter2.knots.transpose() << endl;*/

				  inter1.knots(3) = 0.0001; inter1.knots(inter1.n + 1) = 0.9999;
				  inter2.knots(3) = 0.0001; inter2.knots(inter2.n + 1) = 0.9999;

				  // update coordinate of intermediate control points
				  for (int j = 0; j <= inter1.n; j++) {
					  T_cpts[s_inter1][inter1.knots(j + 2)].fromVectorXd(inter1.controlPw.row(j).transpose());
				  }

				  for (int j = 0; j <= inter2.n; j++) {
					  T_cpts[s_inter2][inter2.knots(j + 2)].fromVectorXd(inter2.controlPw.row(j).transpose());
				  }
			  }

			  // update by T_cpts and initial_cpts
			  for (auto it = T_cpts.begin(); it != T_cpts.end(); ++it) {
				  for (auto it1 = (it->second).begin(); it1 != (it->second).end(); ++it1) {
					  VectorXd delta = (initial_cpts[it->first][it1->first]).toVectorXd() - (it1->second).toVectorXd();
					  if (delta.norm() > error) {
						  error = delta.norm();
					  }
					  T temp;
					  temp.fromVectorXd(delta);
					  (s_map[it->first][it1->first]->data).add(temp);
				  }
			  }

			  cout << "finished intermediate_update()!" << endl;
			  
			  return error;
		  }

		  template<class T>
		  inline void Mesh<T>::skinning_intermediate(
			  const VectorXd & s_knots,
			  const vector<NURBSCurve>& curves,
			  Viewer &viewer)
		  {
			  // 1. calculate sample points by linear interpolate

			  // 2. fit sample points by B-spline using PIA method
			  //    (the control points of B-Spline is the initial X,Y;
			  //     and the B-Spline is the shape to be preserved.)
			  // 3. update V to W by X,Y using formula aX+bW+cY=V
			  // 4. calculate s-curve on T-spline surface, forcing the s-curve to approximate the shape by iteration

			  const int curves_num = curves.size();
			  const int dimension = curves[0].controlPw.cols();
			  map<double, map<double, T>> initial_cpts;

			  intermediate_init(s_knots, curves, initial_cpts, viewer);
			  skinning_update_cross(s_knots, curves);

			  const int max_iter_num = 0;
			  const double eps = 1e-5;
			  int iter_num = 0;
			  double error = 1.0;
			  while (error > eps&&iter_num < max_iter_num) {
				  
				  if (iter_num == max_iter_num - 1) {
					  error = intermediate_update(s_knots, curves, initial_cpts, viewer,100,false);
					  iter_num++;
					  cout << "*******************iter " << iter_num << ": ,error " << error << "********************" << endl;
					  break;
				  }
				  else {
					  error = intermediate_update(s_knots, curves, initial_cpts, viewer);
				  }
				  //error = intermediate_update(s_knots, curves, initial_cpts, viewer);
				  skinning_update_cross(s_knots, curves);
				  iter_num++;
				  cout << "*******************iter " << iter_num << ": ,error " << error << "********************" << endl;
				  //skinning_update_cross(s_knots, curves);

			  }
		  }

		  template<class T>
		  inline void Mesh<T>::skinning_update_cross(const VectorXd & s_knots, const vector<NURBSCurve>& curves)
		  {
			  const int curves_num = curves.size();
			  // 4. update coordinates of control points by the formula from (nasri 2012)
			  // aX' + bW + cY' = V
			  map<double, T> coeff_X; // X'
			  map<double, T> coeff_Y; // Y'

			  for (int i = 1; i <= curves_num - 2; i++) {
				  const double s_now = s_knots(i);
				  auto node = s_map[s_now].begin()->second;
				  auto left_node = node->adj[3];
				  auto right_node = node->adj[1];
				  const double s_left = left_node->s[2];
				  const double s_right = right_node->s[2];

				  double a = Basis((left_node->s).toVectorXd(), s_now);
				  double b = Basis((node->s).toVectorXd(), s_now);
				  double c = Basis((right_node->s).toVectorXd(), s_now);

				  basis_split(s_map[s_left], s_map[s_now], coeff_X); // X'
				  basis_split(s_map[s_right], s_map[s_now], coeff_Y); // Y'

				  for (int j = 0; j <= curves[i].n; j++) {
					  double vknot = curves[i].knots(j + 2);
					  if (j == 1) {
						  vknot = 0.0001;
					  }
					  if (j == curves[i].n - 1) {
						  vknot = 0.9999;
					  }
					  auto temp_X = coeff_X[vknot];
					  auto temp_Y = coeff_Y[vknot];
					  temp_X.scale(-a); // -aX'
					  temp_Y.scale(-c); // -bY'
					  T V;
					  V.fromVectorXd(curves[i].controlPw.row(j).transpose());
					  s_map[s_now][vknot]->data = V.add(temp_X).add(temp_Y).scale(1.0 / b); // W=(V-aX'-bY')/b
				  }

			  }
		  }
		  // nasri 2012 local T-spline skinning
		  // knot vector of every curve should be [0,0,0,0,...,1,1,1,1]
		  // but I need to change multiple knots to single knots [0,0,0,0.0001,...,0.9999,1,1,1]
		  template<class T>
		  inline void Mesh<T>::skinning(const vector<NURBSCurve>& curves, Viewer & viewer)
		  {
			  assert(curves.size() > 0);
			  const int curves_num = curves.size();
			  //const int dimension = curves[0].controlPw.cols();

			  // 1. compute s-knot for curves
			  VectorXd s_knots = skinning_parameterize(curves);

			  // 2. construct basis T-mesh 
			  skinning_init(s_knots, curves);

			  // 3. insert intermediate vertices
			  // the coordinate of vertices is the midpoint of the corresponding points in C_r and C_(r+1)
			  //skinning_insert(s_knots, curves);
			  skinning_insert1(s_knots, curves);
			  skinning_intermediate(s_knots, curves,viewer);

			  // 4. update coordinates of control points by the formula from (nasri 2012)
			  // aX' + bW + cY' = V
			  //skinning_update_cross(s_knots, curves);
		  }

		  

		  template<class T>
		  inline void Mesh<T>::pia_fit(const VectorXd & s_params, const vector<VectorXd>& ttparams, const vector<NURBSCurve>& curves, Viewer &viewer, int max_iterations, double eps)
		  {
			  // 每个网格点对应: Q --> 要逼近的点， params --> 对应的参数值
			  map<double, map<double, T>> Q;
			  map<double, map<double, pair<double, double>>> params;

			  // 初始化T样条的T-preimage,并将控制点坐标设置为要拟合的点坐标Q
			  VectorXd s_nodes(s_params.size());
			  s_nodes(0) = 0.0; s_nodes(1) = 0.0001;
			  s_nodes(s_nodes.size() - 2) = 0.9999; s_nodes(s_nodes.size() - 1) = 1.0;
			  s_nodes.segment(2, s_nodes.size() - 4) = s_params.segment(2, s_nodes.size() - 4);
			  //cout << "s_nodes:" << s_nodes << endl;
			  for (int i = 0; i < s_nodes.size(); i++) {
				  double s = s_nodes(i);
				  double mapped_s = s_params(i);
				  // 取出曲线节点向量的参数部分
				  /*VectorXd t_params = curves[i].knots.segment(3, curves[i].n - 1);
				  t_params[0] = 0.0001; t_params[t_params.size() - 1] = 0.9999;*/
				  VectorXd t_params = ttparams[i];
				  // 计算相应的T样条节点
				  VectorXd t_nodes(t_params.size());
				  t_nodes(0) = 0.0; t_nodes(1) = 0.0001;
				  t_nodes(t_nodes.size() - 2) = 0.9999; t_nodes(t_nodes.size() - 1) = 1.0;
				  t_nodes.segment(2, t_nodes.size() - 4) = t_params.segment(2, t_nodes.size() - 4);

				  for (int j = 0; j < t_params.size(); j++) {
					  double t = t_nodes(j);
					  params[s][t] = make_pair(mapped_s, t_params(j));
					  //cout << "params[" << s << ", " << t << "] = " << mapped_s << ", " << t_params(j) << endl;
					  T point;
					  point.fromVectorXd(curves[i].eval(t_params(j)).row(0).transpose());
					  Q[s][t] = point;
					  MatrixXd P;
					  array2matrixd(point, P);
					  viewer.data().add_points(P, green);
					  insert_helper(s, t, false);
					  Node<T>* node = get_node(s, t);
					  (node->data) = point;
				  }
			  }

			  double error = 100.0;
			  int iter_num = 0;
			  // 开始迭代
			  while (error > eps && iter_num < max_iterations) {
				  error = 0.0;
				  for (auto s_it = s_map.begin(); s_it != s_map.end(); s_it++) {
					  for (auto t_it = (s_it->second).begin(); t_it != (s_it->second).end(); t_it++) {
						  double s = s_it->first;
						  double t = t_it->first;
						  T cur = eval(params[s][t].first, params[s][t].second);
						  T delta = Q[s][t] - cur;
						  /*MatrixXd P1, P2;
						  array2matrixd(Q[s][t], P1);
						  array2matrixd(cur, P2);

						  viewer.data().add_edges(P1, P2, white);*/

						  (s_map[s][t]->data).add(delta);
						  error += delta.toVectorXd().norm();
					  }
				  }

				  error /= get_num();
				  iter_num++;
				  //cout << "iter: " << iter_num << ", error: " << error << endl;
			  }
		  
		  }

		  template<class T>
		  inline double Mesh<T>::maxoffset(double s, int n, const NURBSCurve& curve, double& max_off)
		  {
			  int index = 1;
			  max_off = 0.0;
			  for (int i = 1; i < n; i++) {
				  double t = 1.0*i / n;
				  double offset = (curve.eval(t).row(0).transpose() - eval(s, t).toVectorXd()).norm();
				  if (offset > max_off) {
					  max_off = offset;
					  index = i;
				  }
			  }
			  cout << "maxoffset: " << max_off << endl;
			  return 1.0 * index / n;
		  }

		  template<class T>
		  inline void t_mesh::Mesh<T>::pia_skinning(const vector<NURBSCurve>& curves, Viewer &viewer, int max_iterations, double eps)
		  {
			  assert(curves.size() > 0);
			  const int curves_num = curves.size();
			  //const int dimension = curves[0].controlPw.cols();

			  // 1. compute s-knot for curves
			  // 比如 要拟合点的参数: 0, 0.4, 0.6, 0.8, 1.0
			  // 则对应的样条节点 0, 0,           0, 0, 0.4, 0.6, 0.8, 1.0, 1.0,          1.0, 1.0
			  VectorXd s_params = skinning_parameterize(curves);
			  s_params[0] = 0.0001; s_params[s_params.size() - 1] = 0.9999;


			  vector<VectorXd> ttparams(curves_num);
			  for (int i = 0; i < ttparams.size(); i++) {
				  // 取出曲线节点向量的参数部分
				  VectorXd t_params = curves[i].knots.segment(3, curves[i].n - 1);
				  t_params[0] = 0.0001; t_params[t_params.size() - 1] = 0.9999;
				  ttparams[i] = t_params;
			  }

			  /*for (int i = 0; i < ttparams.size(); i++) {
				  const int n = 10;
				  VectorXd t_params(n + 1);
				  for (int i = 0; i <= n; i++) {
					  t_params(i) = 1.0*i / n;
				  }
				  t_params(0) = 0.0001;
				  t_params(n) = 0.9999;
				  ttparams[i] = t_params;
			  }*/

			  for (int i = 0; i < 10; i++) {
				  this->clear();
				  cout << i << ": *****************************" << endl;
				  /*for (int j = 0; j < ttparams.size(); j++) {
					  cout << "curve " << j << ":    " << ttparams[j].transpose() << endl;
				  }*/
				  pia_fit(s_params, ttparams, curves, viewer, max_iterations, eps);
				  cout << "num of nodes: " << get_num() << endl;
				  const int n = 50; // [0,1]区间分为n份
				  for (int i = 0; i < curves_num; i++) {
					  double s = s_params(i);
					  double max_off = 0.0;
					  double t = maxoffset(s, n, curves[i], max_off);
					  //t_mesh::insert(ttparams[i], t);

					  if (max_off > eps) {
						  t_mesh::insert(ttparams[i], t);
					  }  
				  }
			  }
		  }

		  
	  

    template<class T>
        int Mesh<T>::loadMesh(string name){
            ifstream in(name.c_str());
            if(!in)
                return -1;
            int node_num;
            in>>node_num;
			cout << "node_num: " << node_num << endl;
            for(int i=0;i<node_num;++i){
                new_node();
            }
            for(int i=0;i<node_num;++i){
                nodes[i]->load(in,*this);
                s_map[nodes[i]->s[2]][nodes[i]->t[2]]=nodes[i];
                t_map[nodes[i]->t[2]][nodes[i]->s[2]]=nodes[i];
            }
            return 0;
        }

   

    template<class T>
        int Mesh<T>::saveMesh(string name){
            ofstream out((name+iter_str+".cfg").c_str());
            if(!out)
                return -1;
            int node_num=nodes.size();
            out<<node_num<<endl;
            for(int i=0;i<node_num;++i){
                nodes[i]->save(out);
            }
            return 0;
        }

    template<class T>
        Node<T>* Mesh<T>::new_node(){
            Node<T>* t=new Node<T>(nodes.size()+1);
            nodes.push_back(t);
            return t;
        }

    template<class T>
        Node<T>* Mesh<T>::get_node(int num){
            if(num>(int)nodes.size()||num<=0)
                return 0;
            return nodes[num-1];
        }
    template<class T>
        Node<T>* Mesh<T>::get_node(double s,double t){
            if(s_map.find(s)!=s_map.end()){
                if(s_map[s].find(t)!=s_map[s].end()){
                    return s_map[s][t];
                }
            }
            return 0;
        }

    template<class T>
        Node<T> Mesh<T>::get_knot(double s2,double t2){
            // calculate knot vector in s&t direction
            Node<T> node(0);
            node.s[2]=s2;
            node.t[2]=t2;

            typedef typename map<double,map<double,Node<T>*> >::iterator   map_t;
            typedef typename map<double,Node<T>*>::iterator             map2_t;
            typedef typename map<double,map<double,Node<T>*> >::reverse_iterator   rmap_t;

            // calculate s3,s4 by judging whether [s2+a,t2] (a>0)intersects with s-edge
			int offset = 2;
            for(map_t iter=s_map.begin();iter!=s_map.end();++iter){
                if(iter->first<=s2)
                    continue;
                for(map2_t iter1=(iter->second).begin();iter1!=(iter->second).end();++iter1){
                    if(t2<iter1->second->t[2])
                        break;
                    if(iter1->second->t[2]==t2){
                        node.s[++offset]=iter->first;
                        break;
                    }
                    if(iter1->second->adj[2]!=0){
                        if(iter1->second->t[2]<=t2&&t2<=iter1->second->adj[2]->t[2]){
                            node.s[++offset]=iter->first;
                            break;
                        }
                    }
                }
                if(offset>=4)
                    break;
            }
            while(offset<4)
                node.s[++offset]=width;

            // calculate s1,s0 by judging whether [s2-a,t2] (a>0)intersects with s-edge
            offset=2;
            for(rmap_t iter=s_map.rbegin();iter!=s_map.rend();++iter){
                if(iter->first>=s2)
                    continue;
                for(map2_t iter1=(iter->second).begin();iter1!=(iter->second).end();++iter1){
                    if(t2<iter1->second->t[2])
                        break;
                    if(iter1->second->t[2]==t2){
                        node.s[--offset]=iter->first;
                        break;
                    }
                    if(iter1->second->adj[2]!=0){
                        if(iter1->second->t[2]<=t2&&t2<=iter1->second->adj[2]->t[2]){
                            node.s[--offset]=iter->first;
                            break;
                        }
                    }
                }
                if(offset<=0)
                    break;
            }
            while(offset>0)
                node.s[--offset]=0.0;

            // calculate t3,t4 by judging whether [s2,t2+a] (a>0)intersects with t-edge
            offset=2;
            for(map_t iter=t_map.begin();iter!=t_map.end();++iter){
                if(iter->first<=t2)
                    continue;
                for(map2_t iter1=(iter->second).begin();iter1!=(iter->second).end();++iter1){
                    if(s2<iter1->second->s[2])
                        break;
                    if(iter1->second->s[2]==s2){
                        node.t[++offset]=iter->first;
                        break;
                    }
                    if(iter1->second->adj[1]!=0){
                        if(iter1->second->s[2]<=s2&&s2<=iter1->second->adj[1]->s[2]){
                            node.t[++offset]=iter->first;
                            break;
                        }
                    }
                }
                if(offset>=4)
                    break;
            }
            while(offset<4)
                node.t[++offset]=height;

            // calculate t1,t0 by judging whether [s2,t2-a] (a>0)intersects with t-edge
            offset=2;
            for(rmap_t iter=t_map.rbegin();iter!=t_map.rend();++iter){
                if(iter->first>=t2)
                    continue;
                for(map2_t iter1=(iter->second).begin();iter1!=(iter->second).end();++iter1){
                    if(s2<iter1->second->s[2])
                        break;
                    if(iter1->second->s[2]==s2){
                        node.t[--offset]=iter->first;
                        break;
                    }
                    if(iter1->second->adj[1]!=0){
                        if(iter1->second->s[2]<=s2&&s2<=iter1->second->adj[1]->s[2]){
                            node.t[--offset]=iter->first;
                            break;
                        }
                    }
                }
                if(offset<=0)
                    break;
            }
            while(offset>0)
                node.t[--offset]=0.0;

            return node;
        }

   

   

    template<class T>
        void Mesh<T>::insert(double s,double t){
            if(get_node(s,t)!=0)
                return;

            Node<T> tmp=get_knot(s,t);
            int count=0;
            if(get_node(tmp.s[1],tmp.t[2])!=0)
                ++count;
            if(get_node(tmp.s[3],tmp.t[2])!=0)
                ++count;
            if(get_node(tmp.s[2],tmp.t[1])!=0)
                ++count;
            if(get_node(tmp.s[2],tmp.t[3])!=0)
                ++count;
            if(count==0){
                insert_helper(tmp.s[1],tmp.t[2]);
                merge_all();
                insert_helper(tmp.s[2],tmp.t[1]);
                merge_all();
            }else if(count==1){
                if(get_node(tmp.s[1],tmp.t[2])==0)
                    insert_helper(tmp.s[1],tmp.t[2]);
                else
                    insert_helper(tmp.s[2],tmp.t[1]);
                merge_all();
            }
            insert_helper(tmp.s[2],tmp.t[2]);
            merge_all();
        }
		template<class T>
		inline void Mesh<T>::basis_split(const map<double, Node<T>*>& fewer_map, const map<double, Node<T>*>& more_map, map<double, T>& coeff)
		{
			vector<Node<T>> split_pool;
			// 1. copy node of fewer_map to new map
			//vector<Node<T>> fewer_nodes;
			map<double, Node<T>*> new_map;
			for (auto it = fewer_map.begin(); it != fewer_map.end(); ++it) {
				Node<T>* node = new Node<T>();
				*node = *(it->second);
				new_map[it->first] = node;
			}
			for (auto it = more_map.begin(); it != more_map.end(); ++it) {
				if (new_map.find(it->first) == new_map.end()) {
					// split basis function by knot (it->first)
					for (auto it1 = new_map.begin(); it1 != new_map.end(); ++it1) {
						Node<T> temp;
						if (it1->second->split(0, it->first, &temp, true)) {
							split_pool.push_back(temp);
						}
					}
					// merge same node from split_pool to new_map
					for (int i = 0; i < split_pool.size(); i++) {
						double t = split_pool[i].t[2];
						if (new_map.find(t) != new_map.end()) {
							(new_map[t]->data).add(split_pool[i].data);
						}
						else {
							Node<T>* node = new Node<T>();
							*node = split_pool[i];
							new_map[t] = node;
						}
					}
					//bug: split_pool should be cleared
					split_pool.clear();
					cout << "split_pool size: " << split_pool.size() << endl;
				}
			}

			// return data after spliting by reference and delete node
			coeff.clear();
			for (auto it = new_map.begin(); it != new_map.end(); ++it) {
				coeff[it->first] = it->second->data;
				delete it->second;
			}

		}
		template<class T>
        int Mesh<T>::insert_helper(double s,double t,bool changedata){
            if(get_node(s,t)!=0)
                return 0;
            Node<T>* node=new_node();
            *node=get_knot(s,t);
            node->adj[0]=get_node(node->s[2],node->t[1]);
            node->adj[1]=get_node(node->s[3],node->t[2]);
            node->adj[2]=get_node(node->s[2],node->t[3]);
            node->adj[3]=get_node(node->s[1],node->t[2]);

            if(node->adj[0]){
                node->adj[0]->adj[2]=node;
            }
            if(node->adj[1]){
                node->adj[1]->adj[3]=node;
            }
            if(node->adj[2]){
                node->adj[2]->adj[0]=node;
            }
            if(node->adj[3]){
                node->adj[3]->adj[1]=node;
            }

            s_map[s][t]=node;
            t_map[t][s]=node;

			adjust(node, changedata);
            /*node->s.output(cout);
            node->t.output(cout);
            cout<<endl;*/
            return 1;
        }
     template<class T>
         bool Mesh<T>::check_valid(){
             for(size_t i=0;i<nodes.size();++i){
                 if(nodes[i]->s[2]==0||nodes[i]->t[2]==0||nodes[i]->s[2]>=width||nodes[i]->t[2]>=height)
                     continue;
                 Node<T> tmp=get_knot(nodes[i]->s[2],nodes[i]->t[2]);
                 if(tmp.s!=nodes[i]->s){
                     tmp.save(cout);
                     nodes[i]->save(cout);
                     //outMesh("error");
					 cout << "error: invalid T-mesh *********************!" << endl;
                     return false;
                 }
                 if(tmp.t!=nodes[i]->t){
                     tmp.save(cout);
                     nodes[i]->save(cout);
                     //outMesh("error");
					 cout << "error: invalid T-mesh *********************!" << endl;
                     return false;
                 }
                 if(nodes[i]->adj[0]!=get_node(nodes[i]->s[2],nodes[i]->t[1])){
                     nodes[i]->save(cout);
					 cout << "error: invalid T-mesh *********************!" << endl;
                     return false;
					 
                 }
                 if(nodes[i]->adj[1]!=get_node(nodes[i]->s[3],nodes[i]->t[2])){
                     nodes[i]->save(cout); 
					 cout << "error: invalid T-mesh *********************!" << endl;
                     return false;
                 }
                 if(nodes[i]->adj[2]!=get_node(nodes[i]->s[2],nodes[i]->t[3])){
                     nodes[i]->save(cout); 
					 cout << "error: invalid T-mesh *********************!" << endl;
                     return false;
                 }
                 if(nodes[i]->adj[3]!=get_node(nodes[i]->s[1],nodes[i]->t[2])){
                     nodes[i]->save(cout); 
					 cout << "error: invalid T-mesh *********************!" << endl;
                     return false;
                 }
             }
             return true;
         }

    template<class T>
        void Mesh<T>::merge_all(){
            while(!pool.empty()){
                Node<T> tmp=pool.back();
                pool.pop_back();

				// ��ȡpool�д�������ʱnode��Ӧ����ʵnode
                Node<T>* m_node=get_node(tmp.s[2],tmp.t[2]);
				// �����ʱnode�����ڣ���������node
                if(!m_node){
                    pool.push_back(tmp);
                    insert_helper(tmp.s[2],tmp.t[2]);
                    continue;
                }
				// ����ʱnode����ʵnode ��ȫһ��,��ϲ�����ʵnode
                if(m_node->s==tmp.s&&m_node->t==tmp.t){
                    m_node->data.add(tmp.data);
                    continue;
                }
				// ����ʵnodeΪ [s0,s1,s2,s3,s4]
				// ����ʱnodeΪ[a,b,c,d,e], ������[s0,s1,s2,s3,s4]���½ڵ㴦�������µ�node
				// ͬʱ�µ�node��knot vector֮����ܻ�䣬��Ҫ�ȼ���pool
                bool check=0;
                for(int i=0;i<5;++i){
                    if(i==2)
                        continue;
                    if(tmp.s[i]<=m_node->s[0]||tmp.s[i]>=m_node->s[4])
                        continue;
					
                    if(!m_node->s.have(tmp.s[i])){
                        insert_helper(tmp.s[i],tmp.t[2]);
                        pool.push_back(tmp);
                        check=1;
                        break;
                    }
                }
                if(check)
                    continue; // ��ͣ�����µ�node
				// t����ͬ��
                for(int i=0;i<5;++i){
                    if(i==2)
                        continue;
                    if(tmp.t[i]<=m_node->t[0]||tmp.t[i]>=m_node->t[4])
                        continue;

                    if(!m_node->t.have(tmp.t[i])){
                        insert_helper(tmp.s[2],tmp.t[i]);
                        pool.push_back(tmp);
                        check=1;
                        break;
                    }
                }
                if(check)
                    continue;

				// ��û�в����κε㣬����ʵnodeΪ [s0,s1,s2,s3,s4]
				// ����ʱnodeΪ���� [s0,s2,s4,a,b],[s0,s1,s2,s4,a] ��������ʽ��
				// ��Ҫ����ʵnode�е�knot��ϸ
                if(m_node->merge(&tmp,pool)!=0){
                    pool.push_back(tmp);
                }
            }
        }

    template<class T>
	void Mesh<T>::adjust(Node<T>* n, bool changedata) {
		if (!n)
			return;
		double knots[4] = { n->t[2],n->s[2],n->t[2],n->s[2] };
		// ��4�������8�����ϣ����� blending function refinement
		// for example: N[s0,s1,s2,s3,s4](s) = c1*N[k,s1,s2,s3,s4](s)+c2*N[s0,k,s1,s2,s3](s)
		// B(s2,t2) = N[s0,s1,s2,s3,s4](s)*N[t0,t1,t2,t3,t4](t)
		//          =    c1*N[k,s1,s2,s3,s4](s)*N[t0,t1,t2,t3,t4](t)
		//    		   + c2*N[s0,k,s1,s2,s3](s)*N[t0,t1,t2,t3,t4](t)
		//			= c1*B(s2,t2) + c2*N[s0,k,s1,s2,s3](s)*N[t0,t1,t2,t3,t4](t)
		// ����[s2,t2]����blending function ��ϸ�������� [s2,t2] �� [s1,t2]
		// ��[s1,t2] ����kont vector ��[t0,t1,t2,t3,t4] ��ͬ��
		// ��Ҫ���� violation test
		// violation1: ��[t0,t1,t2,t3,t4] ������ [s1,t2] ����kont vector,�������ϸ
		// violation2: ��[t0,t1,t2,t3,t4] �д��� [s1,t2] �� knot vectorû�еĽڵ㣬������µĵ�
		for (int i = 0; i < 4; ++i) {
			if (n->adj[i]) {
				Node<T> tmp2;
				if (n->adj[i]->split(i, knots[i], &tmp2, changedata)) {
					pool.push_back(tmp2);
				}
				if (n->adj[i]->adj[i]) {
					Node<T> tmp2;
					if (n->adj[i]->adj[i]->split(i, knots[i], &tmp2, changedata)) {
						pool.push_back(tmp2);
					}
				}
			}
		}

		set<Node<T>*>   node_set;
		typedef typename map<double, map<double, Node<T>*> >::iterator   map_t;
		typedef typename map<double, Node<T>*>::iterator             map2_t;
		typedef typename set<Node<T>*>::iterator                 set_t;
		// �ҳ�s����������knot vector���ܲ���ȷ��node,����node_set
		for (map_t iter = s_map.begin(); iter != s_map.end(); ++iter) {
			if (iter->first < n->s[1])
				continue;
			if (iter->first > n->s[3])
				break;
			for (map2_t iter1 = (iter->second).begin(); iter1 != (iter->second).end(); ++iter1) {
				if (n->t[2] > iter1->second->t[4])
					continue;
				if (n->t[2] < iter1->second->t[0])
					break;
				node_set.insert(iter1->second);
			}
		}
		// �ҳ�s����������knot vector���ܲ���ȷ��node, ����node_set
		for (map_t iter = t_map.begin(); iter != t_map.end(); ++iter) {
			if (iter->first < n->t[1])
				continue;
			if (iter->first > n->t[3])
				break;
			for (map2_t iter1 = (iter->second).begin(); iter1 != (iter->second).end(); ++iter1) {
				if (n->s[2] > iter1->second->s[4])
					continue;
				if (n->s[2] < iter1->second->s[0])
					break;
				node_set.insert(iter1->second);
			}
		}
		//��ÿ��blending function��Ҫ��ϸ��node,  ���м�ϸ 
		for (set_t iter = node_set.begin(); iter != node_set.end(); ++iter) {
			Node<T> tmp = get_knot((*iter)->s[2], (*iter)->t[2]);
			if (tmp.s != (*iter)->s) {
				(*iter)->valid = false;
				for (int j = 0; j < 5; ++j) {
					if (j == 2)
						continue;
					// ����ȷ��get_knot�Ľڵ��ϣ�����ϸ 
					if (!(*iter)->s.have(tmp.s[j])) {
						Node<T> tmp2;
						int dir = j > 2 ? 1 : 3;
						if ((*iter)->split(dir, tmp.s[j], &tmp2, changedata)) {
							pool.push_back(tmp2);
						}
					}
				}
			}
			if (tmp.t != (*iter)->t) {
				(*iter)->valid = false;
				for (int j = 0; j < 5; ++j) {
					if (j == 2)
						continue;
					if (!(*iter)->t.have(tmp.t[j])) {
						Node<T> tmp2;
						int dir = j > 2 ? 0 : 2;
						if ((*iter)->split(dir, tmp.t[j], &tmp2, changedata)) {
							pool.push_back(tmp2);
						}
					}
				}
			}
		}
	}

};