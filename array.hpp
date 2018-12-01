#include"utility.h"
#include<cstring>

namespace t_mesh{
    using namespace std;

    template<class T,int num>
    struct Array{
        public:
            Array():data(){}
			Array(const Array<T, num>& a) {
				for (int i = 0; i < num; i++) {
					data[i] = a[i];
				}
			}
            typedef T ValueType;
            enum{SIZE=num};
            T& operator[](int index){
                if(index<0||index>=num)
                    throw("out of range");
                return data[index];
            }
            const T& operator[](int index) const{
                if(index<0||index>=num)
                    throw("out of range");
                return data[index];
            }
            bool operator==(const Array<T,num>& other){
                for(int i=0;i<num;++i){
                    if(data[i]!=other.data[i])
                        return false;
                }
                return true;
            }
            bool operator!=(const Array<T,num>& other){
                return !(*this==other);
            }
            bool have(const T& t){
                for(int i=0;i<num;++i){
                    if(data[i]==t)
                        return true;
                }
                return false;
            }
            void output(ostream& out) const;
            void input(istream& in);
            Array<T,num>& scale(const T&);
            Array<T,num>& add(const Array<T,num>&);
            void clear();
        private:
            T data[num];
    };
    template<class T,int num>
    void Array<T,num>::output(ostream& out) const{
        for(int i=0;i<num;++i){
            out<<data[i]<<' ';
        }
    }

    template<class T,int num>
    void Array<T,num>::input(istream& in){
        for(int i=0;i<num;++i){
            in>>data[i];
        }
    }

    template<class T,int num>
    Array<T,num>& Array<T,num>::scale(const T& x) {
        for(int i=0;i<num;++i){
            data[i]*=x;
        }
        return *this;
    }
    template<class T,int num>
    Array<T,num>& Array<T,num>::add(const Array<T,num>& other){
        for(int i=0;i<num;++i){
            data[i]+=other.data[i];
        }
        return *this;
    }
    template<class T,int num>
    void Array<T,num>::clear(){
        memset(data,0,sizeof(data));
    }
};
