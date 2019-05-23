#include"utility.h"
#include<list>

namespace t_mesh{
    using namespace std;
    
    template<class T>
    class Node{
            friend class Mesh<T>;
        public:
            Node(int num=0):order(num),data(),s(),t(),adj(),valid(true){}
            ~Node();
            const T& get() const{ return data; }
            void set(const T& v){ data=v; }
            Node<T>& operator=(const Node<T>& n);
            std::istream& load(istream& in,Mesh<T>& m);
            std::ostream& save(ostream& out) const;
            int get_order(){return order;}
            int split(int dir,double k,Node<T>* tmp,bool changedata=true);
            int merge(Node<T>* tmp,list<Node<T> >&pool);
			bool is_ok(double x, double y) {
				return (x >= s[0] && x <= s[4] && y >= t[0] && y <= t[4]);
			}
			double basis(double x, double y) {
				if (is_ok(x, y)) {
					return Basis(s.toVectorXd(), x) * Basis(t.toVectorXd(), y);
				}
				return 0.0;
			}

        public:
            int                 order;   // the order of node: 1,2,3,...
            T                   data;	 // point coordinate, 节点插入通过split能够正确更新点的坐标
            Array<double,5>        s;       // knot vector of s-direction
            Array<double,5>        t;       // knot vector of s-direction
            Array<Node<T>*,4>   adj;     // pointer of lower, right,upper,left node

            bool                valid;
    };

    template<class T>
    Node<T>::~Node(){

    }

    template<class T>
    Node<T>& Node<T>::operator=(const Node<T>& n)
    {
        s=n.s;
        t=n.t;
        data=n.data;
        return *this;
    }
    
    template<class T>
    std::istream& Node<T>::load(istream& in,Mesh<T>& m){
        s.input(in);
        t.input(in);
        for(int i=0;i<4;++i){
            int tmp;
            in>>tmp;
            adj[i]=m.get_node(tmp);
        }
        data.input(in);
        return in;
    }

    template<class T>
    std::ostream& Node<T>::save(ostream& out) const{
        s.output(out);
        t.output(out);
        for(int i=0;i<4;++i){
            if(adj[i])
                out<<adj[i]->get_order()<<' ';
            else
                out<<"0 ";
        }
        data.output(out);
        out<<endl;
        return out;
    }
    // insert k to [s0,s1,s2,s3,s4] or [t0,t1,t2,t3,t4], get temp node
    template<class T>
    int Node<T>::split(int dir,double k,Node<T>* tmp,bool changedata){
        // dir=1,3, s-direction; dir=0,2, t-direction
        //valid=false;
        Array<double,5>& this_st=dir%2?this->s:this->t;
        Array<double,5>& tmp_st=dir%2?tmp->s:tmp->t;

        if(this_st.have(k)||k<this_st[0]||k>this_st[4])
            return 0;
        
        tmp->s=this->s;
        tmp->t=this->t;
        tmp->data=this->data;

        // [s0,k,s1,s2,s3,s4]
        if(k<this_st[1]){
            double s1=1.0*(k-this_st[0])/(this_st[3]-this_st[0]);
			if (changedata) {
				tmp->data.scale(s1);
			}
            
            tmp_st[0]=this_st[0];
            tmp_st[1]=k;
            tmp_st[2]=this_st[1];
            tmp_st[3]=this_st[2];
            tmp_st[4]=this_st[3];
            this_st[0]=k;
        // [s0,s1,k,s2,s3,s4]
        }else if(k<this_st[2]){
            double s1=1.0*(k-this_st[0])/(this_st[3]-this_st[0]);
            double s2=1.0*(this_st[4]-k)/(this_st[4]-this_st[1]);
			if (changedata) {
				tmp->data.scale(s1);
				this->data.scale(s2);
			}
            
            tmp_st[0]=this_st[0];
            tmp_st[1]=this_st[1];
            tmp_st[2]=k;
            tmp_st[3]=this_st[2];
            tmp_st[4]=this_st[3];
            this_st[0]=this_st[1];
            this_st[1]=k;
        // [s0,s1,s2,k,s3,s4]
        }else if(k<this_st[3]){
            double s1=1.0*(this_st[4]-k)/(this_st[4]-this_st[1]);
            double s2=1.0*(k-this_st[0])/(this_st[3]-this_st[0]);
			if (changedata) {
				tmp->data.scale(s1);
				this->data.scale(s2);
			}
            
            tmp_st[0]=this_st[1];
            tmp_st[1]=this_st[2];
            tmp_st[2]=k;
            tmp_st[3]=this_st[3];
            tmp_st[4]=this_st[4];
            this_st[4]=this_st[3];
            this_st[3]=k;
        // [s0,s1,s2,s3,k,s4]
        }else if(k<this_st[4]){
            double s1=1.0*(this_st[4]-k)/(this_st[4]-this_st[1]);
			if (changedata) {
				tmp->data.scale(s1);
			}
            
            tmp_st[0]=this_st[1];
            tmp_st[1]=this_st[2];
            tmp_st[2]=this_st[3];
            tmp_st[3]=k;
            tmp_st[4]=this_st[4];
            this_st[4]=k;
        }
        return 1;
    }
    template<class T>
    int Node<T>::merge(Node<T>* tmp,list<Node<T> >&pool){
        if(tmp->s==s&&tmp->t==t){
            data.add(tmp->data);
            return 0;
        }
        //valid=false;
        for(int i=0;i<5;++i){
            if(!tmp->s.have(s[i])){
                Node<T> tmp2;
                int dir=i>2?1:3;
                if(tmp->split(dir,s[i],&tmp2)){
                    pool.push_back(tmp2);
                }
                return 1;
            }
        }
        for(int i=0;i<5;++i){
            if(!tmp->t.have(t[i])){
                Node<T> tmp2;
                int dir=i>2?2:0;
                if(tmp->split(dir,t[i],&tmp2)){
                    pool.push_back(tmp2);
                }
                return 1;
            }
        }
        return 1;
    }
};
