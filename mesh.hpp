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

namespace t_mesh{
    using namespace std;
	using namespace igl::opengl::glfw;
    template<class T>
    class Mesh{
		private:
			Mesh& operator=(const Mesh&){}
        public:
			Mesh():width(1.0),height(1.0),id(-1){}
			Mesh(const Mesh& other);// deep copy
			~Mesh();
            int loadMesh(string);
			istream& loadMesh(istream&);
            int saveMesh(string);
			ostream& saveMesh(ostream&);
			void saveAsObj(string, double resolution=0.01);
			T eval(double s,double t);
			vector<std::tuple<double, double, double, double>> region(double u, double v);
			void draw(bool tmesh, bool polygon, bool surface,double resolution = 0.01);
			void setViewer(Viewer* viewer) { this->viewer = viewer; }
			void piafit(const map<double, map<double, T>>& targetPoints, int maxIterNum=10, double eps=1e-5);
			int get_num() const { return nodes.size(); }

			T du(double u, double v);
			T dv(double u, double v);
			T normal(double u, double v);


			double du2(Node<T>* node, double u, double v);
			double dv2(Node<T>* node, double u, double v);
			double duv(Node<T>* node, double u, double v);
          
            void insert(double s,double t);
			int insert_helper(double s, double t, bool changedata = true);

            Node<T>*    new_node();
            Node<T>*    get_node(int num);
            Node<T>*    get_node(double s,double t);
            Node<T>     get_knot(double x,double y); // get the knot vector of (x,y)
			bool check_valid();
			void merge_all();
        private:
			void drawTmesh();
			void drawControlpolygon();
			void drawSurface(double resolution = 0.01);

			void adjust(Node<T>* n, bool changedata = true);
            
            
			void clear();
			
		public:
			// organizing node in a good data structure 
			map<double, map<double, Node<T>*> > s_map; // s_map[s][t]
			map<double, map<double, Node<T>*> > t_map; // t_map[t][s]
			list<Node<T> >              pool;
			vector<Node<T>*>    nodes;

		private:
			double width;
			double height;
           
            string          iter_str;

			Eigen::MatrixXd mesh_V;
			Eigen::MatrixXi mesh_F;
			Viewer* viewer;
			int id;
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
	inline double Mesh<T>::du2(Node<T>* node, double u, double v) {
		if (u == 0.0) {
			u = 0.0001;
		}
		if (u == 1.0) {
			u = 0.9999;
		}
		return DersBasis(node->s.toVectorXd(), u)(2)*Basis(node->t.toVectorXd(), v);
	}
	template<class T>
	inline double Mesh<T>::dv2(Node<T>* node, double u, double v) {
		if (v == 0.0) {
			v = 0.0001;
		}
		if (v == 1.0) {
			v = 0.9999;
		}
		return DersBasis(node->t.toVectorXd(), v)(2)*Basis(node->s.toVectorXd(), u);
	}
	template<class T>
	inline double Mesh<T>::duv(Node<T>* node, double u, double v) {
		if (u == 0.0) {
			u = 0.0001;
		}
		if (u == 1.0) {
			u = 0.9999;
		}
		if (v == 0.0) {
			v = 0.0001;
		}
		if (v == 1.0) {
			v = 0.9999;
		}
		return DersBasis(node->s.toVectorXd(), u)(1)*DersBasis(node->t.toVectorXd(), v)(1);
	}
	template<class T>
		 T Mesh<T>::eval(double s, double t) {
			 T result;
			 for (int i = 0; i < nodes.size(); i++) {
				 double blend = nodes[i]->basis(s, t);
				 T temp(nodes[i]->data);
				 //temp.output(cout);
				 temp.scale(blend);
				 result.add(temp);
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
		 inline vector<std::tuple<double, double, double, double>> Mesh<T>::region(double u, double v)
		 {
			 vector<tuple<double, double, double, double>> res;
			 // find the rectangle region of param u,v
			 if (get_node(u, v) != 0) {
				 return res;
			 }
			 Node<T> temp = get_knot(u, v);
			 if (get_node(temp.s[2], temp.t[1]) != 0 && get_node(temp.s[2], temp.t[3]) != 0) {
				 // on the u-edge of tspline

				 // 左边的矩形
				 double ulow = temp.s[1];
				 double u = temp.s[2];
				 double vlow = 0;
				 double vhigh = 0;
				 int count = 0;
				 for (auto it = s_map[u].begin(); it != s_map[u].end(); it++) {
					 if (u == ulow) {
						 break;
					 }
					 if (it->first >= temp.t[2]) {
						 break;
					 }
					 if (it->second->adj[3] != NULL && it->second->adj[3]->s[2] == ulow) {
						 vlow = it->first;
						 count++;
					 }
				 }
				 if (count != 0) {
					 count = 0;
					 for (auto it = s_map[u].rbegin(); it != s_map[u].rend(); it++) {
						 if (it->first <= temp.t[2]) {
							 break;
						 }
						 if (it->second->adj[3] != NULL && it->second->adj[3]->s[2] == ulow) {
							 vhigh = it->first;
							 count++;
						 }
					 }
					 if (count != 0) {
						 res.push_back(make_tuple(ulow, u, vlow, vhigh));
					 }
				 }


				 // 右边的矩形
				 double uhigh = temp.s[3];
				 count = 0;
				 for (auto it = s_map[u].begin(); it != s_map[u].end(); it++) {
					 if (u == uhigh) {
						 break;
					 }
					 if (it->first >= temp.t[2]) {
						 break;
					 }
					 if (it->second->adj[1] != NULL && it->second->adj[1]->s[2] == uhigh) {
						 vlow = it->first;
						 count++;
					 }
				 }
				 if (count != 0) {
					 count = 0;
					 for (auto it = s_map[u].rbegin(); it != s_map[u].rend(); it++) {
						 if (it->first <= temp.t[2]) {
							 break;
						 }
						 if (it->second->adj[1] != NULL && it->second->adj[1]->s[2] == uhigh) {
							 vhigh = it->first;
							 count++;
						 }
					 }
					 if (count != 0) {
						 res.push_back(make_tuple(u, uhigh, vlow, vhigh));
					 }
				 }
			 }
			 else if (get_node(temp.s[1], temp.t[2]) != 0 && get_node(temp.s[3], temp.t[2]) != 0) {
				 // on the v-edge of tspline
				 // 下面的矩形
				 double vlow = temp.t[1];
				 double v = temp.t[2];
				 double ulow = 0;
				 double uhigh = 0;
				 int count = 0;
				 for (auto it = t_map[v].begin(); it != t_map[v].end(); it++) {
					 if (v == vlow) {
						 break;
					 }
					 if (it->first >= temp.s[2]) {
						 break;
					 }
					 if (it->second->adj[0] != NULL && it->second->adj[0]->t[2] == vlow) {
						 ulow = it->first;
						 count++;
					 }
				 }
				 if (count != 0) {
					 count = 0;
					 for (auto it = t_map[v].rbegin(); it != t_map[v].rend(); it++) {
						 if (it->first <= temp.s[2]) {
							 break;
						 }
						 if (it->second->adj[0] != NULL && it->second->adj[0]->t[2] == vlow) {
							 uhigh = it->first;
							 count++;
						 }
					 }
					 if (count != 0) {
						 res.push_back(make_tuple(ulow, uhigh, vlow, v));
					 }
				 }


				 // 上面的矩形
				 double vhigh = temp.t[3];
				 count = 0;
				 for (auto it = t_map[v].begin(); it != t_map[v].end(); it++) {
					 if (v == vhigh) {
						 break;
					 }
					 if (it->first >= temp.s[2]) {
						 break;
					 }
					 if (it->second->adj[2] != NULL && it->second->adj[2]->t[2] == vhigh) {
						 ulow = it->first;
						 count++;
					 }
				 }
				 if (count != 0) {
					 count = 0;
					 for (auto it = t_map[v].rbegin(); it != t_map[v].rend(); it++) {
						 if (it->first <= temp.s[2]) {
							 break;
						 }
						 if (it->second->adj[2] != NULL && it->second->adj[2]->t[2] == vhigh) {
							 uhigh = it->first;
							 count++;
						 }
					 }
					 if (count != 0) {
						 res.push_back(make_tuple(ulow, uhigh, v, vhigh));
					 }
				 }
			 }
			 else {
				 // inside the rectangle
				 res.push_back(std::make_tuple(temp.s[1], temp.s[3], temp.t[1], temp.t[3]));
			 }
			 return res;
		 }

	template<class T>
		 void Mesh<T>::drawTmesh(){
			 assert(viewer != NULL); // use setViewer(Viewer* viewer)

			 Eigen::MatrixXd P1(1,2), P2(1,2);
			 Eigen::MatrixXd nodes_st(nodes.size(), 2);
			 
			 //(*viewer).data().add_label(Eigen::Vector3d(0, 0, 0), "haha");
			 //cout << "1" << endl;
			 for (int i = 0; i < nodes.size(); i++) {
				 nodes_st.row(i) << nodes[i]->s[2], nodes[i]->t[2];
				 std::stringstream label;
				 label <<nodes_st(i, 0) << ", " << nodes_st(i, 1);
				
				 (*viewer).data().add_label(nodes_st.row(i), label.str());
				 (*viewer).data().add_points(nodes_st.row(i), red);
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
						 (*viewer).data().add_edges(P1, P2, white);
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
						 (*viewer).data().add_edges(P1, P2,white);
					 }
				 }
			 }
			 (*viewer).core.align_camera_center(nodes_st); // center
		 }

	template<class T>
	      void Mesh<T>::drawControlpolygon() {
			  assert(viewer != NULL); // use setViewer(Viewer* viewer)

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
						  (*viewer).data().add_edges(P1, P2, green);
					  }
				  }
			  }

			  for (auto iter = t_map.begin(); iter != t_map.end(); ++iter) {
				  for (auto iter1 = (iter->second).begin(); iter1 != (iter->second).end(); ++iter1) {
					  if ((iter1->second)->adj[1]) {
						  array2matrixd(iter1->second->data, P1);
						  array2matrixd((iter1->second->adj[1])->data, P2);
						  (*viewer).data().add_edges(P1, P2, blue);
					  }
				  }
			  }
			  (*viewer).data().add_points(nodes_point, red);
			  (*viewer).core.align_camera_center(nodes_point); // center
			  
		  }

	template<class T>
	      void Mesh<T>::drawSurface(double resolution) {
			  assert(viewer != NULL); // use setViewer(Viewer* viewer)
		
			  // cut apart the parameter domain
			  //double u_low = (++s_map.begin())->first;
			  //double u_high = (++s_map.rbegin())->first;
			  
			  double u_low = s_map.begin()->first;
			  double u_high = s_map.rbegin()->first;

			  cout << "u_low: " << u_low <<", u_high: "<<u_high<< endl;
			  const int uspan = (u_high - u_low) / resolution;
			  double u_resolution = (u_high - u_low) / uspan;
			  
			  /*double v_low = (++t_map.begin())->first;
			  double v_high = (++t_map.rbegin())->first;*/
			  double v_low = t_map.begin()->first;
			  double v_high = t_map.rbegin()->first;
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
					  
					  /*Eigen::MatrixXd endPoint;
					  endPoint.setZero();

					  array2matrixd(normal(u, v), endPoint);
					  endPoint.row(0) = endPoint.row(0).normalized();
					  cout << "normal: \n" << endPoint << endl;
					  endPoint += curvePoint;
					  (*viewer).data().add_edges(curvePoint, endPoint, red);*/

					  /*array2matrixd(dv(u, v), endPoint);
					  endPoint.row(0) = endPoint.row(0).normalized();
					  cout << "normal: \n" << endPoint << endl;
					  endPoint += curvePoint;
					  (*viewer).data().add_edges(curvePoint, endPoint, blue);*/
				  }

			  for (int j = 0; j<vspan; j++)
				  for (int i = 0; i < uspan; i++){
					  int V_index = j*(uspan + 1) + i;
					  int F_index = 2 * j*uspan + 2 * i;
					  mesh_F.row(F_index) << V_index, V_index + 1, V_index + uspan + 1;
					  mesh_F.row(F_index + 1) << V_index + uspan + 1, V_index + 1, V_index + uspan + 2;
				  }

			  Eigen::MatrixXd HN;
			  Eigen::VectorXd H;
			  Eigen::SparseMatrix<double> L, M, Minv;
			  igl::cotmatrix(mesh_V, mesh_F, L);
			  igl::massmatrix(mesh_V, mesh_F, igl::MASSMATRIX_TYPE_VORONOI, M);
			  igl::invert_diag(M, Minv);
			  HN = -Minv*(L*mesh_V);
			  H = HN.rowwise().norm(); //up to sign

			  // compute curvatrue directions via quadric fitting
			  Eigen::MatrixXd PD1, PD2;
			  Eigen::VectorXd PV1, PV2;
			  igl::principal_curvature(mesh_V, mesh_F, PD1, PD2, PV1, PV2);
			  // mean curvature
			  H = 0.5*(PV1 + PV2);

			  (*viewer).data().set_mesh(mesh_V, mesh_F);

			  // Compute pseudocolor
			  Eigen::MatrixXd C;
			  igl::jet(H, true, C);
			  (*viewer).data().set_colors(C);
		  }

	template<class T>
		  void Mesh<T>::draw(bool tmesh, bool polygon, bool surface, double resolution){

			  assert(viewer != NULL); // use setViewer(Viewer* viewer)
			  if (id == -1) {
				  id = (*viewer).append_mesh();
			  }
			  (*viewer).selected_data_index = id;
			  if (tmesh) {
				  drawTmesh();
				  return;
			  }
			  if (polygon) {
				  drawControlpolygon();
			  }
			  if (surface) {
				  drawSurface(resolution);
			  }
		  }

		  template<class T>
		  inline void Mesh<T>::piafit(const map<double, map<double, T>>& targetPoints, int maxIterNum, double eps)
		  {
			  // 初始控制点坐标设为要逼近的目标点坐标
			  for (auto it = targetPoints.begin(); it != targetPoints.end(); it++) {
				  for (auto it1 = (it->second).begin(); it1 != (it->second).end(); it1++) {
					  double s = it->first;
					  double t = it1->first;
					  s_map[s][t]->data = targetPoints[s][t];
				  }
			  }

			  
			  // 迭代更新控制顶点
			  for (int i = 0; i < maxIterNum; i++) {
				  double error = 0.0;
				  for (auto it = targetPoints.begin(); it != targetPoints.end(); it++) {
					  for (auto it1 = (it->second).begin(); it1 != (it->second).end(); it1++) {
						  double s = it->first;
						  double t = it1->first;
						  T delta = targetPoints[s][t] - eval(s, t);
						  (s_map[s][t]->data).add(delta);
						  error = max(error, delta.toVectorXd().norm());
					  }
				  }
				  cout << "iter: " << i+1 << ", error: " << error << endl;
				  if (error < eps) {
					  break;
				  }
			  }
		  }

		  template<class T>
		  inline T Mesh<T>::du(double u, double v)
		  {
			  if (u == 0.0) {
				  u = 0.0001;
			  }
			  if (u == 1.0) {
				  u = 0.9999;
			  }
			  T res;
			  for (int i = 0; i < nodes.size(); i++) {
				  auto node = nodes[i];
				  double fac = DersBasis(node->s.toVectorXd(), u)(1) * Basis(node->t.toVectorXd(), v);
				  res.add(fac * node->data);
			  }
			  return res;
		  }

		  template<class T>
		  inline T Mesh<T>::dv(double u, double v)
		  {
			  if (v == 0.0) {
				  v = 0.0001;
			  }
			  if (v == 1.0) {
				  v = 0.9999;
			  }
			  T res;
			  for (int i = 0; i < nodes.size(); i++) {
				  auto node = nodes[i];
				  double fac = DersBasis(node->t.toVectorXd(), v)(1) * Basis(node->s.toVectorXd(), u);
				  res.add(fac * node->data);
			  }
			  return res;
		  }

		  template<class T>
		  inline T Mesh<T>::normal(double u, double v)
		  {
			  Eigen::Vector3d a = du(u, v).toVectorXd();
			  Eigen::Vector3d b = dv(u, v).toVectorXd();
			  Eigen::Vector3d n = a.cross(b);
			  n.normalize();
			  T res;
			  res.fromVectorXd(n);
			  return res;
		  }

		  template<class T>
		  inline Mesh<T>::Mesh(const Mesh & other)
		  {
			  this->nodes = other.nodes;
			  for (int i = 0; i < other.nodes.size(); i++) {
				  Node<T>* node = new Node<T>(*(other.nodes[i]));
				  this->nodes[i] = node;
				  this->s_map[node->s[2]][node->t[2]] = node;
				  this->t_map[node->t[2]][node->s[2]] = node;
			  }
			  // 更改邻接关系
			  for (int i = 0; i < this->nodes.size(); i++) {
				  Node<T>* node = this->nodes[i];
				  for (int j = 0; j < 4; j++) {
					  Node<T>* other_node = other.nodes[i]->adj[j];
					  if (other_node) {
						  node->adj[j] = this->get_node(other_node->order);
					  }
					  else {
						  node->adj[j] = NULL;
					  }
				  }
			  }
			  this->viewer = other.viewer;
			  this->mesh_V = other.mesh_V;
			  this->mesh_F = other.mesh_F;
			  this->pool = other.pool;
			  this->id = -1;
			  this->width = 1.0;
			  this->height = 1.0;
		  }

		  template<class T>
		  inline Mesh<T>::~Mesh()
		  {
			  // 析构所有new出的node
			  for (int i = 0; i < nodes.size(); i++) {
				  delete nodes[i];
				  nodes[i] = NULL;
			  }
		  }

		  template<class T>
		  int Mesh<T>::loadMesh(string name) {
			  ifstream in(name);
			  if (!in.is_open()) {
				  cout << "can't open file: " << name << endl;
				  return -1;
			  }
				  
			  loadMesh(in);
			  return 0;
		  }

		  template<class T>
		  inline istream & Mesh<T>::loadMesh(istream& in)
		  {
			  int node_num = 0;
			  in >> node_num;
			  cout << "node_num: " << node_num << endl;
			  for (int i = 0; i<node_num; ++i) {
				  new_node();
			  }
			  for (int i = 0; i<node_num; ++i) {
				  nodes[i]->load(in, *this);
				  s_map[nodes[i]->s[2]][nodes[i]->t[2]] = nodes[i];
				  t_map[nodes[i]->t[2]][nodes[i]->s[2]] = nodes[i];
			  }
			  return in;
		  }

   

    template<class T>
        int Mesh<T>::saveMesh(string name){
			ofstream out(name + ".cfg");
			if (!out.is_open()) {
				cout << "can't open file: " << name << ".cfg" << endl;
				return -1;
			}
			saveMesh(out);
            return 0;
        }

		template<class T>
		inline ostream & Mesh<T>::saveMesh(ostream& out)
		{
			int node_num = nodes.size();
			out << node_num << endl;
			for (int i = 0; i<node_num; ++i) {
				nodes[i]->save(out);
			}
			return out;
		}

		template<class T>
		inline void Mesh<T>::saveAsObj(string filename, double resolution)
		{
			filename += ".obj";
			ofstream out(filename);
			if (!out.is_open()) {
				cout << "error, can't open file: " << filename << endl;
				return;
			}

			double u_low = s_map.begin()->first;
			double u_high = s_map.rbegin()->first;

			cout << "u_low: " << u_low << ", u_high: " << u_high << endl;
			const int uspan = (u_high - u_low) / resolution;
			double u_resolution = (u_high - u_low) / uspan;
			double v_low = t_map.begin()->first;
			double v_high = t_map.rbegin()->first;
			const int vspan = (v_high - v_low) / resolution;
			double v_resolution = (v_high - v_low) / vspan;
			cout << "v_low: " << v_low << ", v_high: " << v_high << endl;
			mesh_V = Eigen::MatrixXd((uspan + 1)*(vspan + 1), 3);
			mesh_F = Eigen::MatrixXi(2 * uspan*vspan, 3);
			// discretize T-Spline Surface into triangular mesh(V,F) in libigl mesh structure
			// calculate 

			for (int j = 0; j <= vspan; j++)
				for (int i = 0; i <= uspan; i++) {
					Eigen::MatrixXd curvePoint;
					double u = u_low + i*u_resolution;
					double v = v_low + j*v_resolution;
					array2matrixd(eval(u_low + i*u_resolution, v_low + j*v_resolution), curvePoint);
					mesh_V.row(j*(uspan + 1) + i) = curvePoint;
				}

			for (int j = 0; j<vspan; j++)
				for (int i = 0; i < uspan; i++) {
					int V_index = j*(uspan + 1) + i;
					int F_index = 2 * j*uspan + 2 * i;
					mesh_F.row(F_index) << V_index, V_index + 1, V_index + uspan + 1;
					mesh_F.row(F_index + 1) << V_index + uspan + 1, V_index + 1, V_index + uspan + 2;
				}

			igl::writeOBJ(filename, mesh_V, mesh_F);
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