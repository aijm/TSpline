namespace t_mesh{
    template<class T>
    struct Point{
        T           data;
        T           origin;
        Parameter   parameter;
        double       get_error(){
            double x=0;
            for(int i=0;i<T::SIZE;++i){
                x+=(data[i]-origin[i])*(data[i]-origin[i]);
            }
            return x;
        }
        T           get_diff(){
            T tmp=origin;
            tmp.add(data.scale(-1));
            return tmp;
        }
   };

    template<class T>
    struct Linker{
        public:
            Linker(Point<T>*& p,Node<T>*& n):valid(false),point(p),node(n){
            }
            double get();
            void invalid();
            bool is_valid(){return valid;}
            void fit_point();
            // judging whether point is in the region of node 
            bool is_ok(){return point->parameter[0]>=node->s[0]&&
                                point->parameter[0]<=node->s[4]&&
                                point->parameter[1]>=node->t[0]&&
                                point->parameter[1]<=node->t[4];
            }
        private:
            void update();

        private:
            bool        valid;
            double       value;
        public:
            Point<T>*   point;
            Node<T>*    node;
    };

    template<class T>
    double Linker<T>::get(){
        if(!valid)
            update();
        return value;
    }
    template<class T>
    void Linker<T>::invalid(){
        valid=false;
    }

    template<class T>
    void Linker<T>::update(){
        value=node->base(point->parameter);
        valid=true;
    }

    template<class T>
    void Linker<T>::fit_point(){
        double tmp=this->get();
        const T& pt=node->get();
        point->data.add(pt.scale(tmp));
    }
};
