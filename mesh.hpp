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
    template<class T>
    class Mesh{
        public:
            int loadMesh(string);
            int saveMesh(string);
			T eval(double s,double t);
			void drawTmesh(igl::opengl::glfw::Viewer &viewer);
			void drawControlpolygon(igl::opengl::glfw::Viewer &viewer);
			void drawSurface(igl::opengl::glfw::Viewer &view, double resolution = 0.01);
			void draw(igl::opengl::glfw::Viewer &viewer, bool tmesh, bool polygon, bool surface,double resolution = 0.01);
			int get_num() const { return nodes.size(); }
			void skinning(const vector<NURBSCurve> &curves);
            /*int savePoints(string);
            int saveLog(string);
            int loadPoints(string);
            void outMesh(string);*/

            //int init();
            /*int get_iter_num();
            void fit_points(int thread_num);*/
            void insert(double s,double t);
            /*void update_mesh();
            int piafit(int thread_num);
            void solve();
            int solve2(int num,double e);
            void solve3();
            void solve4();
            double get_error();*/
            // void init_mesh();
            //void test(int threads);

            Node<T>*    new_node();
            Node<T>*    get_node(int num);
            Node<T>*    get_node(double s,double t);
            Node<T>     get_knot(double x,double y); // get the knot vector of (x,y)

			Eigen::MatrixXd mesh_V;
			Eigen::MatrixXi mesh_F;
        private:
            int insert_helper(double s,double t);
            void adjust(Node<T>* n);
            void merge_all();
            //bool check_valid();
            /*void update_iter();
            void update_log();
            void pia_thread(int c,int a);
            void fit_thread(int c,int a);
            double get_A(const vector<vector<pair<int,double> > >& A,int i,int j);*/

            vector<Node<T>*>    nodes;
            /*cv::Mat             origin;
            cv::Mat             data;*/

            // organizing node in a good data structure 
            map<double,map<double,Node<T>*> > s_map; // s_map[s][t]
            map<double,map<double,Node<T>*> > t_map; // t_map[t][s]
            list<Node<T> >              pool;

			double width = 3.0;
			double height = 3.0;
            // int             iter_num;
            string          iter_str;
            // int             width;
            // int             height;
            // Array<int,2>    offset;
            /*double          error;
            ostringstream   logger;*/
			
    };
	
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
		 void Mesh<T>::drawTmesh(igl::opengl::glfw::Viewer &viewer){
			 Eigen::MatrixXd P1(1,2), P2(1,2);
			 Eigen::MatrixXd nodes_st(nodes.size(), 2);
			 typedef typename map<double, map<double, Node<T>*> >::iterator   map_t;
			 typedef typename map<double, Node<T>*>::iterator             map2_t;
			
			 //cout << "1" << endl;
			 for (int i = 0; i < nodes.size(); i++) {
				 nodes_st.row(i) << nodes[i]->s[2], nodes[i]->t[2];
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
	      void Mesh<T>::drawControlpolygon(igl::opengl::glfw::Viewer &viewer) {
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
	      void Mesh<T>::drawSurface(igl::opengl::glfw::Viewer &viewer, double resolution) {
			  // cut apart the parameter domain
			  double u_low = (++s_map.begin())->first;
			  //u_low = 1.0;
			  double u_high = (++s_map.rbegin())->first;
			  //u_high = 2.0;
			  const int uspan = (u_high - u_low) / resolution;
			  double u_resolution = (u_high - u_low) / uspan;
			  
			  double v_low = (++t_map.begin())->first;
			  //v_low = 1.0;
			  double v_high = (++t_map.rbegin())->first;
			  //v_high = 2.0;
			  const int vspan = (v_high - v_low) / resolution;
			  double v_resolution = (v_high - v_low) / vspan;

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
		  void Mesh<T>::draw(igl::opengl::glfw::Viewer & viewer,bool tmesh, bool polygon, bool surface, double resolution){
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

		  // nasri 2012 local T-spline skinning
		  template<class T>
		  inline void Mesh<T>::skinning(const vector<NURBSCurve>& curves)
		  {
			  // 1. construct basis T-mesh 


			  // 2. insert intermediate vertices by knot insertion

			  // 3. update coordinates of cross-sectional NURBS curves by the formula from (nasri 2012)

		  }

		  
	  
//    template<class T>
//        void Mesh<T>::test(int threads){
//            vector<T>   backup;
//            for(int i=0;i<(int)nodes.size();++i){
//                backup.push_back(nodes[i]->data);
//            }
//            piafit(threads);
//            double e=get_error();
//            for(int i=0;i<(int)nodes.size();++i){
//                nodes[i]->data=backup[i];
//            }
//            fit_points(threads);
//            int num=solve2(0,e);
//            for(int i=0;i<(int)nodes.size();++i){
//                nodes[i]->data=backup[i];
//            }
//            fit_points(threads);
//            solve2(num,0);
//        }
//
//    template<class T>
//        void Mesh<T>::update_iter(){
//            stringstream ss;
//            ss<<iter_num;
//            ss>>iter_str;
//        }
//
//    template<class T>
//        void Mesh<T>::update_log(){
//            logger.str("");
//            logger<<"Iter:"<<iter_num<<endl;
//            logger<<"width:"<<width<<endl;
//            logger<<"height:"<<height<<endl;
//        }
//
//    template<class T>
//        int Mesh<T>::get_iter_num(){
//            return iter_num;
//        }
//    
//    template<class T>
//        void Mesh<T>::solve(){
//            boost::timer t;
//            int count=0;
//            cv::Mat A=cv::Mat::zeros(nodes.size(),nodes.size(),CV_64FC1);
//            cv::Mat B0=cv::Mat::zeros(nodes.size(),1,CV_64FC1);
//            cv::Mat B1=cv::Mat::zeros(nodes.size(),1,CV_64FC1);
//            cv::Mat B2=cv::Mat::zeros(nodes.size(),1,CV_64FC1);
//            cv::Mat X0=cv::Mat::zeros(nodes.size(),1,CV_64FC1);
//            cv::Mat X1=cv::Mat::zeros(nodes.size(),1,CV_64FC1);
//            cv::Mat X2=cv::Mat::zeros(nodes.size(),1,CV_64FC1);
//            for(int i=0;i<(int)nodes.size();++i){
//                for(int j=i;j<(int)nodes.size();++j){
//                    for(int y=max(nodes[i]->t[0],nodes[j]->t[0]);y<min(nodes[i]->t[4],nodes[j]->t[4]);++y){
//                        for(int x=max(nodes[i]->s[0],nodes[j]->s[0]);x<min(nodes[i]->s[4],nodes[j]->s[4]);++x){
//                            A.at<double>(i,j)+=nodes[i]->base(x+1,y+1)*nodes[j]->base(x+1,y+1);
//                        }
//                    }
//                    A.at<double>(j,i)=A.at<double>(i,j);
//                    if(A.at<double>(i,j)<1e-5){
//                        ++count;
//                        if(i!=j)
//                            ++count;
//                    }
//                    //cout<<A.at<double>(i,j)<<',';
//                }
////                cout<<endl;
//            }
//            cout<<"Total: "<<nodes.size()*nodes.size()<<" Zero:"<<count<<endl;
//            for(int i=0;i<(int)nodes.size();++i){
//                for(int y=nodes[i]->t[0];y<nodes[i]->t[4];++y){
//                    for(int x=nodes[i]->s[0];x<nodes[i]->s[4];++x){
//                        B0.at<double>(i,0)+=nodes[i]->base(x+1,y+1)*origin.ptr<unsigned char>(y)[(x)*(T::SIZE)+0];
//                        B1.at<double>(i,0)+=nodes[i]->base(x+1,y+1)*origin.ptr<unsigned char>(y)[(x)*(T::SIZE)+1];
//                        B2.at<double>(i,0)+=nodes[i]->base(x+1,y+1)*origin.ptr<unsigned char>(y)[(x)*(T::SIZE)+2];
//                    }
//                }
//            }
//            cout<<cv::solve(A,B0,X0,cv::DECOMP_CHOLESKY)<<endl;
//            cout<<cv::solve(A,B1,X1,cv::DECOMP_CHOLESKY)<<endl;
//            cout<<cv::solve(A,B2,X2,cv::DECOMP_CHOLESKY)<<endl;
//            for(int i=0;i<(int)nodes.size();++i){
//                nodes[i]->data[0]=X0.at<double>(i,0);
//                nodes[i]->data[1]=X1.at<double>(i,0);
//                nodes[i]->data[2]=X2.at<double>(i,0);
//                //    cout<<X0.at<double>(i,0)<<','<<X1.at<double>(i,0)<<','<<X2.at<double>(i,0)<<endl;
//            }
//            cout<<t.elapsed()<<endl;
//        }
//    
//    template<class T>
//        int Mesh<T>::solve2(int num,double e){
//            time_t start,end1,end2;
//            time(&start);
//            vector<boost::unordered_map<int,double> > A;
//
//            cout<<"initial matrix...."<<endl;
//                
//            vector<T> B;
//            B.reserve(nodes.size());
//            for(int i=0;i<(int)nodes.size();++i){
//                boost::unordered_map<int,double> mymap;
//                for(int j=i;j<(int)nodes.size();++j){
//                    for(int y=max(nodes[i]->t[0],nodes[j]->t[0]);y<min(nodes[i]->t[4],nodes[j]->t[4]);++y){
//                        for(int x=max(nodes[i]->s[0],nodes[j]->s[0]);x<min(nodes[i]->s[4],nodes[j]->s[4]);++x){
//                            mymap[j]+=nodes[i]->base(x+1,y+1)*nodes[j]->base(x+1,y+1);
//                       }
//                    }
//               }
//                A.push_back(mymap);
//            }
//
//            cout<<"initial A ok..."<<endl;
//            for(int i=0;i<(int)nodes.size();++i){
//                T tmp;
//                for(int y=nodes[i]->t[0];y<nodes[i]->t[4];++y){
//                    for(int x=nodes[i]->s[0];x<nodes[i]->s[4];++x){
//                        for(int k=0;k<T::SIZE;++k){
//                            tmp[k]+=nodes[i]->base(x+1,y+1)*origin.ptr<unsigned char>(y)[(x)*(T::SIZE)+k];
//                        }
//                    }
//                }
//                B.push_back(tmp);
//            }
//
//            cout<<"initial B ok..."<<endl;
//            
//            if(num==0){
//                int x=0;
//                do{
//                    for(int i=0;i<(int)nodes.size();++i){
//                        T last=nodes[i]->data;
//                        nodes[i]->data=B[i];
//                        for(int j=0;j<(int)nodes.size();++j){
//                            if(i==j)
//                                continue;
//                            double Aij=0;
//                            boost::unordered_map<int,double>::iterator iter;
//                            if(i>j){
//                                if((iter=A[j].find(i))!=A[j].end())
//                                    Aij=iter->second;
//                            }else{
//                                if((iter=A[i].find(j))!=A[i].end())
//                                    Aij=iter->second;
//                            }
//                            if(Aij>1e-3){
//                                T tmp=nodes[j]->data;
//                                nodes[i]->data.add(tmp.scale(-Aij));
//                            }
//                        }
//                        nodes[i]->data.scale(1/A[i][i]);
//                    }
//                    fit_points(8);
//                    cout<<"iter:"<<++x<<":"<<get_error()<<endl;
//                    if(get_error()<e)
//                        return x;
//                }while(1);
//            }else{
//                time(&end1);
//                logger<<"Gauss-seidel init time:"<<difftime(end1,start)<<endl;
//                int x=0;
//                do{
//                    double diff=0;
//                    for(int i=0;i<(int)nodes.size();++i){
//                        T last=nodes[i]->data;
//                        nodes[i]->data=B[i];
//                        for(int j=0;j<(int)nodes.size();++j){
//                            if(i==j)
//                                continue;
//                            double Aij=0;
//                            boost::unordered_map<int,double>::iterator iter;
//                            if(i>j){
//                                if((iter=A[j].find(i))!=A[j].end())
//                                    Aij=iter->second;
//                            }else{
//                                if((iter=A[i].find(j))!=A[i].end())
//                                    Aij=iter->second;
//                            }
//                            if(Aij>1e-3){
//                                T tmp=nodes[j]->data;
//                                nodes[i]->data.add(tmp.scale(-Aij));
//                            }
//                        }
//                        nodes[i]->data.scale(1/A[i][i]);
//                        for(int k=0;k<T::SIZE;++k){
//                            diff+=(last[k]-nodes[i]->data[k])*(last[k]-nodes[i]->data[k]);
//                        }
//                    }
//                    cout<<"iter "<<++x<<":"<<diff<<endl;
//                    if(x>=num)
//                        break;
//                }while(1);
//                time(&end2);
//                logger<<"Gauss-seidel iter time:"<<difftime(end2,end1)<<endl;
//                logger<<"Gauss-seidel steps:"<<x<<endl;
//                logger<<"Gauss-seidel mean iter:"<<difftime(end2,end1)/x<<endl;
//                logger<<"Gauss-seidel total time:"<<difftime(end2,start)<<endl;
//                fit_points(8);
//                logger<<"Gauss-seidel after error:"<<get_error()<<endl;
//            }
//            return 0;
//        }
//
//    template<class T>
//        double Mesh<T>::get_A(const vector<vector<pair<int,double> > >& A,int i,int j){
//            if(i>j){
//                int tmp=i;
//                i=j;
//                j=tmp;
//            }
//            int lp=0,rp=A[i].size()-1;
//            if(A[i][lp].first>j||A[i][rp].first<j)
//                return 0;
//            while(lp<=rp){
//                if(A[i][(lp+rp)/2].first==j)
//                    return A[i][(lp+rp)/2].second;
//                if(A[i][(lp+rp)/2].first<j){
//                    if(lp==(lp+rp)/2)
//                        return 0;
//                    lp=(lp+rp)/2;
//                }else{
//                    if(rp==(lp+rp)/2)
//                        return 0;
//                    rp=(lp+rp)/2;
//                }
//            }
//            return 0;
//        }
//
//    template<class T>
//        void Mesh<T>::solve3(){
//            boost::timer t;
//            vector<T> B;
//            B.reserve(nodes.size());
//
//            vector<vector<pair<int,double> > > A;
//            A.reserve(nodes.size());
//
//            for(int i=0;i<(int)nodes.size();++i){
//                vector<pair<int,double> > vec;
//                vec.reserve(nodes.size());
//                for(int j=i;j<(int)nodes.size();++j){
//                    double tmp=0;
//                    for(int y=max(nodes[i]->t[0],nodes[j]->t[0]);y<min(nodes[i]->t[4],nodes[j]->t[4]);++y){
//                        for(int x=max(nodes[i]->s[0],nodes[j]->s[0]);x<min(nodes[i]->s[4],nodes[j]->s[4]);++x){
//                            tmp+=nodes[i]->base(x+1,y+1)*nodes[j]->base(x+1,y+1);
//                        }
//                    }
//                    if(tmp>0){
//                        vec.push_back(pair<int,double>(j,tmp));
//                    }
//                }
//                A.push_back(vec);
//            }
//            for(int i=0;i<(int)nodes.size();++i){
//                T tmp;
//                for(int y=nodes[i]->t[0];y<nodes[i]->t[4];++y){
//                    for(int x=nodes[i]->s[0];x<nodes[i]->s[4];++x){
//                        for(int k=0;k<T::SIZE;++k){
//                            tmp[k]+=nodes[i]->base(x+1,y+1)*origin.ptr<unsigned char>(y)[(x)*(T::SIZE)+k];
//                        }
//                    }
//                }
//                B.push_back(tmp);
//            }
//
//            int x=0;
//            do{
//                double diff=0;
//                for(int i=0;i<(int)nodes.size();++i){
//                    T last=nodes[i]->data;
//                    nodes[i]->data=B[i];
//                    for(int j=0;j<(int)nodes.size();++j){
//                        if(i==j)
//                            continue;
//                        double b=get_A(A,i,j);
//                        if(b>0){
//                            T tmp=nodes[j]->data;
//                            nodes[i]->data.add(tmp.scale(-b));
//                        }
//                    }
//                    nodes[i]->data.scale(1/get_A(A,i,i));
//                    for(int k=0;k<T::SIZE;++k){
//                        diff+=(last[k]-nodes[i]->data[k])*(last[k]-nodes[i]->data[k]);
//                    }
//                }
//                cout<<"iter "<<++x<<":"<<diff<<endl;
//                if(diff<1)
//                    break;
//            }while(1||x<100);
//            cout<<t.elapsed()<<endl;
//        }
//    template<class T>
//        void Mesh<T>::solve4(){
//            boost::timer t;
//            int count=0;
//            cv::Mat A=cv::Mat::zeros(nodes.size(),nodes.size(),CV_64FC1);
//
//            vector<T> B;
//            B.reserve(nodes.size());
//            for(int i=0;i<(int)nodes.size();++i){
//                boost::unordered_map<int,double> mymap;
//                for(int j=i;j<(int)nodes.size();++j){
//                    for(int y=max(nodes[i]->t[0],nodes[j]->t[0]);y<min(nodes[i]->t[4],nodes[j]->t[4]);++y){
//                        for(int x=max(nodes[i]->s[0],nodes[j]->s[0]);x<min(nodes[i]->s[4],nodes[j]->s[4]);++x){
//                            A.at<double>(i,j)+=nodes[i]->base(x+1,y+1)*nodes[j]->base(x+1,y+1);
//                        }
//                    }
//
//                    A.at<double>(j,i)=A.at<double>(i,j);
//                }
//            }
//            cout<<"Total: "<<nodes.size()*nodes.size()<<" Zero:"<<count<<endl;
//            for(int i=0;i<(int)nodes.size();++i){
//                T tmp;
//                for(int y=nodes[i]->t[0];y<nodes[i]->t[4];++y){
//                    for(int x=nodes[i]->s[0];x<nodes[i]->s[4];++x){
//                        for(int k=0;k<T::SIZE;++k){
//                            tmp[k]+=nodes[i]->base(x+1,y+1)*origin.ptr<unsigned char>(y)[(x)*(T::SIZE)+k];
//                        }
//                    }
//                }
//                B.push_back(tmp);
//            }
//
//            int x=0;
//            do{
//                double diff=0;
//                for(int i=0;i<(int)nodes.size();++i){
//                    T last=nodes[i]->data;
//                    nodes[i]->data=B[i];
//                    for(int j=0;j<(int)nodes.size();++j){
//                        if(i==j)
//                            continue;
//                        T tmp=nodes[j]->data;
//                        nodes[i]->data.add(tmp.scale(-A.at<double>(i,j)));
//                    }
//                    nodes[i]->data.scale(1/A.at<double>(i,i));
//                    for(int k=0;k<T::SIZE;++k){
//                        diff+=(last[k]-nodes[i]->data[k])*(last[k]-nodes[i]->data[k]);
//                    }
//                }
//                cout<<"iter "<<++x<<":"<<diff<<endl;
//                if(diff<1)
//                    break;
//            }while(1);
//            cout<<t.elapsed()<<endl;
//        }

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

    /*template<class T>
        int Mesh<T>::loadPoints(string img){
            origin=cv::imread(img,CV_64FC3);
            width=origin.size().width;
            height=origin.size().height;

            return 0;
        }*/

    /*template<class T>
        int Mesh<T>::savePoints(string name){
            cv::imwrite(name+iter_str+".jpg",data);

            return 0;
        }*/

    /*template<class T>
        int Mesh<T>::saveLog(string name){
            ofstream out((name+iter_str+".log").c_str());
            if(!out)
                return -1;
            out<<logger.str();
            return 0;
        }*/

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
                node.s[++offset]=width+1;

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
                node.s[--offset]=-1.0;

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
                node.t[++offset]=height+1;

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
                node.t[--offset]=-1.0;

            return node;
        }

    /*template<class T>
        int Mesh<T>::init(){
            for(size_t i=0;i<nodes.size();++i){
                nodes[i]->valid=false;
                nodes[i]->update();
            }
            return 0;
        }*/

    // template<class T>
    //     void Mesh<T>::init_mesh(){
    //         for(int y=0;y<height;y+=100){
    //             for(int x=0;x<width;x+=100){
    //                 Node<T> tmp=get_knot(x+1,y+1);
    //                 insert_helper(tmp.s[1],tmp.t[2]);
    //                 merge_all();
    //                 insert_helper(tmp.s[2],tmp.t[1]);
    //                 merge_all();
    //                 insert_helper(tmp.s[3],tmp.t[2]);
    //                 merge_all();
    //                 insert_helper(tmp.s[2],tmp.t[3]);
    //                 merge_all();
    //                 insert_helper(tmp.s[2],tmp.t[2]);
    //                 merge_all();
    //             }
    //         }
    //         init();
    //     }

    //template<class T>
    //    void Mesh<T>::fit_thread(int current,int thread_num){
    //        for(size_t i=0;i<nodes.size();++i){
    //            for(int y=(nodes[i]->t[0]+1);y<=nodes[i]->t[4]-1;++y){
    //                if(y<height*current/thread_num){
    //                    y=height*current/thread_num;
    //                    continue;
    //                }
    //                if(y>height*(current+1)/thread_num)
    //                    break;
    //                for(int x=nodes[i]->s[0]+1;x<=nodes[i]->s[4]-1;++x){
    //                    T tmp=nodes[i]->get();
    //                    double b=nodes[i]->base(x,y);
    //                    for(int n=0;n<T::SIZE;++n){
    //                        data.ptr<double>(y-1)[(x-1)*(T::SIZE)+n]+=tmp[n]*b;
    //                    }
    //                }
    //            }
    //        }
    //    }

    //template<class T>
    //    void Mesh<T>::fit_points(int thread_num){
    //        using namespace boost;
    //        {

    //            data=cv::Mat::zeros(origin.size(),CV_64FC3);

    //            if(thread_num>1){
    //                thread_group threads;
    //                for(int i=0;i<thread_num;++i){
    //                    threads.create_thread(bind(&Mesh<T>::fit_thread,this,i,thread_num));
    //                }
    //                threads.join_all();
    //            }
    //            else{
    //                fit_thread(0,1);
    //            }
    //        }

    //        error=0;
    //        for(int y=0;y<height;++y){
    //            for(int x=0;x<width;++x){
    //                double e=0;
    //                for(int n=0;n<T::SIZE;++n){
    //                    e+=(origin.ptr<unsigned char>(y)[(x)*(T::SIZE)+n]-data.ptr<double>(y)[(x)*(T::SIZE)+n])*
    //                        (origin.ptr<unsigned char>(y)[(x)*(T::SIZE)+n]-data.ptr<double>(y)[(x)*(T::SIZE)+n]);
    //                }
    //                error+=e/width/height;
    //            }
    //        }
    //    }

    //template<class T>
    //    void Mesh<T>::update_mesh(){
    //        ++iter_num;
    //        update_iter();
    //        update_log();

    //        multimap<double,Array<int,2> > offsets;

    //        
    //        double rate=sqrt(iter_num);
    //        int a=rate*height/100;
    //        int b=rate*width/100;

    //        int last_node_num=nodes.size();
    //        logger<<"Before insert:"<<last_node_num<<endl;

    //        time_t start,end;
    //        time(&start);
    //        
    //        for(int i=0;i<a;++i){
    //            for(int j=0;j<b;++j){
    //                double max_error=0;
    //                double total_error=0;
    //                int count=0;
    //                for(int y=height*i/a;y<height*(i+1)/a;++y){
    //                    for(int x=width*j/b;x<width*(j+1)/b;++x){
    //                        double e=0;
    //                        for(int n=0;n<T::SIZE;++n){
    //                            e+=(origin.ptr<unsigned char>(y)[(x)*(T::SIZE)+n]-data.ptr<double>(y)[(x)*(T::SIZE)+n])*
    //                                (origin.ptr<unsigned char>(y)[(x)*(T::SIZE)+n]-data.ptr<double>(y)[(x)*(T::SIZE)+n]);
    //                        }
    //                        if(get_node(x+1,y+1)!=0)
    //                            continue;
    //                        if(e>max_error){
    //                            max_error=e;
    //                            offset[0]=x+1;
    //                            offset[1]=y+1;
    //                        }
    //                        total_error+=e;
    //                        ++count;
    //                    }
    //                }
    //                total_error/=count;
    //                offsets.insert(pair<double,Array<int,2> >(total_error,offset));
    //            }
    //        }

    //        multimap<double,Array<int,2> >::reverse_iterator iter;

    //        int x=0;
    //        for(iter=offsets.rbegin();iter!=offsets.rend();++iter){
    //            Node<T> tmp=get_knot((iter->second)[0],(iter->second)[1]);
    //            insert((tmp.s[1]+tmp.s[3])/2,(tmp.t[1]+tmp.t[3])/2);
    //            //insert(tmp.s[2],tmp.t[2]);
    //            if(++x>a*b*0.1)
    //                break;
    //        }

    //        time(&end);
    //        logger<<"Insert nodes:"<<nodes.size()-last_node_num<<endl;
    //        logger<<"After insert:"<<nodes.size()<<endl;
    //        logger<<"Insertion time:"<<difftime(end,start)<<endl;
    //        init();
    //    }

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
        int Mesh<T>::insert_helper(double s,double t){
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

            adjust(node);
            node->s.output(cout);
            node->t.output(cout);
            cout<<endl;
            return 1;
        }
    // template<class T>
    //     bool Mesh<T>::check_valid(){
    //         for(size_t i=0;i<nodes.size();++i){
    //             if(nodes[i]->s[2]==0||nodes[i]->t[2]==0||nodes[i]->s[2]>=width||nodes[i]->t[2]>=height)
    //                 continue;
    //             Node<T> tmp=get_knot(nodes[i]->s[2],nodes[i]->t[2]);
    //             if(tmp.s!=nodes[i]->s){
    //                 tmp.save(cout);
    //                 nodes[i]->save(cout);
    //                 outMesh("error");
    //                 return false;
    //             }
    //             if(tmp.t!=nodes[i]->t){
    //                 tmp.save(cout);
    //                 nodes[i]->save(cout);
    //                 outMesh("error");
    //                 return false;
    //             }
    //             if(nodes[i]->adj[0]!=get_node(nodes[i]->s[2],nodes[i]->t[1])){
    //                 nodes[i]->save(cout);                   
    //                 return false;
    //             }
    //             if(nodes[i]->adj[1]!=get_node(nodes[i]->s[3],nodes[i]->t[2])){
    //                 nodes[i]->save(cout);                   
    //                 return false;
    //             }
    //             if(nodes[i]->adj[2]!=get_node(nodes[i]->s[2],nodes[i]->t[3])){
    //                 nodes[i]->save(cout);                   
    //                 return false;
    //             }
    //             if(nodes[i]->adj[3]!=get_node(nodes[i]->s[1],nodes[i]->t[2])){
    //                 nodes[i]->save(cout);                   
    //                 return false;
    //             }
    //         }
    //         return true;
    //     }
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
        void Mesh<T>::adjust(Node<T>* n){
            if(!n)
                return;
            double knots[4]={n->t[2],n->s[2],n->t[2],n->s[2]};
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
            for(int i=0;i<4;++i){
                if(n->adj[i]){
                    Node<T> tmp2;
                    if(n->adj[i]->split(i,knots[i],&tmp2)){
                        pool.push_back(tmp2); 
                    }
                    if(n->adj[i]->adj[i]){
                        Node<T> tmp2;
                        if(n->adj[i]->adj[i]->split(i,knots[i],&tmp2)){
                            pool.push_back(tmp2); 
                        }
                    }
                }
            }

            set<Node<T>*>   node_set;
            typedef typename map<double,map<double,Node<T>*> >::iterator   map_t;
            typedef typename map<double,Node<T>*>::iterator             map2_t;
            typedef typename set<Node<T>*>::iterator                 set_t;
			// �ҳ�s����������knot vector���ܲ���ȷ��node,����node_set
            for(map_t iter=s_map.begin();iter!=s_map.end();++iter){
                if(iter->first<n->s[1])
                    continue;
                if(iter->first>n->s[3])
                    break;
                for(map2_t iter1=(iter->second).begin();iter1!=(iter->second).end();++iter1){
                    if(n->t[2]>iter1->second->t[4])
                        continue;
                    if(n->t[2]<iter1->second->t[0])
                        break;
                    node_set.insert(iter1->second);
                }
            }
			// �ҳ�s����������knot vector���ܲ���ȷ��node, ����node_set
            for(map_t iter=t_map.begin();iter!=t_map.end();++iter){
                if(iter->first<n->t[1])
                    continue;
                if(iter->first>n->t[3])
                    break;
                for(map2_t iter1=(iter->second).begin();iter1!=(iter->second).end();++iter1){
                    if(n->s[2]>iter1->second->s[4])
                        continue;
                    if(n->s[2]<iter1->second->s[0])
                        break;
                    node_set.insert(iter1->second);
                }
            }
			//��ÿ��blending function��Ҫ��ϸ��node,  ���м�ϸ 
            for(set_t iter=node_set.begin();iter!=node_set.end();++iter){
                Node<T> tmp=get_knot((*iter)->s[2],(*iter)->t[2]);
                if(tmp.s!=(*iter)->s){
                    (*iter)->valid=false;
                    for(int j=0;j<5;++j){
                        if(j==2)
                            continue;
						// ����ȷ��get_knot�Ľڵ��ϣ�����ϸ 
                        if(!(*iter)->s.have(tmp.s[j])){
                            Node<T> tmp2;
                            int dir=j>2?1:3;
                            if((*iter)->split(dir,tmp.s[j],&tmp2)){
                                pool.push_back(tmp2);
                            }
                        }
                    }
                }
                if(tmp.t!=(*iter)->t){
                    (*iter)->valid=false;
                    for(int j=0;j<5;++j){
                        if(j==2)
                            continue;
                        if(!(*iter)->t.have(tmp.t[j])){
                            Node<T> tmp2;
                            int dir=j>2?0:2;
                            if((*iter)->split(dir,tmp.t[j],&tmp2)){
                                pool.push_back(tmp2);
                            }
                        }
                    }
                }
            }
        }

    //template<class T>
    //    void Mesh<T>::pia_thread(int current,int thread_num){
    //        for(size_t i=nodes.size()*current/thread_num;i<nodes.size()*(current+1)/thread_num;++i){
    //            double sum1=0;
    //            T   sum2;
    //            for(int y=nodes[i]->t[0]+1;y<=nodes[i]->t[4]-1;++y){
    //                for(int x=nodes[i]->s[0]+1;x<=nodes[i]->s[4]-1;++x){
    //                    double b=nodes[i]->base(x,y);
    //                    T tmp;

    //                    sum1+=b;
    //                    for(int n=0;n<T::SIZE;++n){
    //                        tmp[n]=(origin.ptr<unsigned char>(y-1)[(x-1)*(T::SIZE)+n]-
    //                                data.ptr<double>(y-1)[(x-1)*(T::SIZE)+n])*b;
    //                    }
    //                    sum2.add(tmp);
    //                }
    //            }
    //            assert(sum1>0);
    //            //if(sum1>m)
    //            //m=sum1;
    //            //vec_sum2.push_back(sum2);
    //            sum2.scale(1/sum1);
    //            nodes[i]->data.add(sum2);
    //        }
    //    }

    //template<class T>
    //    int Mesh<T>::piafit(int thread_num){
    //        time_t start,end;
    //        time (&start);
    //        double init_error=get_error();

    //        double last_error;
    //        int x=0;
    //        for(size_t i=0;i<nodes.size();++i){
    //            nodes[i]->valid=false;
    //        }
    //        using namespace boost;

    //        do{
    //            last_error=get_error();
    //            if(thread_num>1){
    //                thread_group threads;
    //                for(int i=0;i<thread_num;++i){
    //                    threads.create_thread(bind(&Mesh<T>::pia_thread,this,i,thread_num));
    //                }
    //                threads.join_all();
    //            }
    //            else{
    //                pia_thread(0,1);
    //            }
    //            fit_points(thread_num);
    //            std::cout<<"pia:"<<++x<<' '<<get_error()<<endl;

    //            if(fabs(get_error()/last_error-1)<1e-3||x>50){
    //                time (&end);
    //                double totaltime=difftime(end,start);
    //                logger<<"Pia total time:"<<totaltime<<endl;
    //                logger<<"Pia steps:"<<x<<endl;
    //                logger<<"Pia mean time:"<<totaltime/x<<endl;
    //                logger<<"Pia before error:"<<init_error<<endl;
    //                logger<<"Pia after error:"<<get_error()<<endl;
    //                return x;
    //            }
    //        }while(1);
    //    }

    //template<class T>
    //    double Mesh<T>::get_error(){
    //        return error;
    //    }

    //template<class T>
    //    void Mesh<T>::outMesh(string name){
    //        IplImage *out=cvCreateImage(cvSize(width+2,height+2),8,1);
    //        cvZero(out);
    //        double a=1,b=1;
    //        typedef typename map<int,map<int,Node<T>*> >::iterator mm_t;
    //        typedef typename map<int,Node<T>*>::iterator m_t;
    //        for(mm_t iter=s_map.begin();iter!=s_map.end();++iter){
    //            for(m_t it=(iter->second).begin();it!=(iter->second).end();++it){
    //                if((it->second)->adj[2]){
    //                    cvLineAA(out,cvPoint((iter->first)/a+b,(it->second)->t[2]/a+b),cvPoint((iter->first)/a+b,(it->second)->adj[2]->t[2]/a+b),255);
    //                }
    //            }
    //        }
    //        for(mm_t iter=t_map.begin();iter!=t_map.end();++iter){
    //            for(m_t it=(iter->second).begin();it!=(iter->second).end();++it){
    //                if((it->second)->adj[1]){
    //                    cvLineAA(out,cvPoint((it->second)->s[2]/a+b,(iter->first)/a+b),cvPoint((it->second)->adj[1]->s[2]/a+b,(iter->first)/a+b),255);
    //                }
    //            }
    //        }
    //        /*
    //           for(size_t i=0;i<nodes.size();++i){
    //           CvPoint   centerpoint;
    //           centerpoint.x=nodes[i]->s[2]/a+b;
    //           centerpoint.y=nodes[i]->t[2]/a+b;
    //           cvCircle( out, centerpoint ,3 , CV_RGB(255,255,255),1, 8, 3);
    //           }
    //           */

    //        cvSaveImage((name+iter_str+".jpg").c_str(),out);
    //        cvReleaseImage(&out);
    //    }

};