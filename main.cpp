#include <iostream>
#include <cmath>

enum HALF{
    LEFT = -1 ,
    RIGHT = 1
};

enum WAVETYPE{
    UDAR,
    RAZR
};


class Param{
public:
    Param();
    Param( float p,  float ro, float u, float X , float B,HALF type);
    friend float Uud(Param& left,Param &right);
    friend float Uraz(Param& left, Param& right);
    friend float Uvac(Param &left , Param& right);
    friend float deltU(Param &left, Param &right);
    friend float Pstart(Param &left, Param &right);
    friend float Func(float P, Param &param);
    friend float dFunc(float P , Param& param);
    friend void Speed(float P,Param& p);
    friend float SpeedContact(Param&left, Param&right);
    friend void Adiabad(float U, float P , Param& p);
 
    void Print();
    void setType(WAVETYPE type);
private:
    float u,p,ro,eps , X, B, c , a , D ,R = 0; 
    HALF course;
    WAVETYPE type;
};

Param::Param(){
        u=0;
        p=0;
        ro=0;
        eps=0;
        X = 0;
        B = 0;
        c= 0;
}




Param::Param(float p,  float ro, float u,  float X, float B, HALF cource){
        this->u=u;
        this->p=p;
        this->ro=ro;
        this->X = X;
        this->B = B;
        //this->eps = (p+X*B)/(X-1)/ro;
        this->c = sqrtf(X*(p+B)/ro);
        this->course = course;
}


float Uud(Param &left, Param &right){
    return (right.p-left.p)/sqrtf(left.ro*((left.X+1.f)/2.f*(right.p+left.B)+(left.X-1.f)/2.f*(left.p+left.B)));
}

float Uraz(Param &left, Param &right){
    return -2.f*right.c/(right.X-1.f)*(1-powf((left.p+right.B)/(right.p+left.B), (right.X - 1.f)/(2.f*right.X) ));
}

float Uvac(Param &left, Param &right){
    return -2.f*left.c/(left.X-1.f)-2.f*right.c*(right.X-1.f);
}

float deltU(Param &left, Param &right){
    return left.u-right.u;
}

float Pstart(Param &left, Param &right){
    return (left.p*right.ro*right.c + right.p*left.ro*left.c+deltU(left,right)*left.ro*left.c*right.ro*right.c)/(left.ro*left.c+right.ro*right.c);
}

float Func(float P, Param &param){
    float pik = (P+param.B)/(param.p+param.B);
    //float ck = sqrtf(param.X*(param.p + param.B)/param.ro);
    if (P >= param.p){
        return (P-param.p)/(param.ro*param.c*sqrtf((param.X+1.f)/2.f/param.X*pik+(param.X-1.f)/2.f/param.X));
    }else
    {
        return 2.f*param.c/(param.X-1.f)*(powf(pik,(param.X-1.f)/2.f/param.X)-1.f);
    }
}

float dFunc(float P , Param& param){
    float pik = (P+param.B)/(param.p+param.B);
    if (P>=param.p)
    {
        return ((param.X+1.f)*pik+3.f*param.X-1.f)/(4.f*param.X*param.ro*param.c*sqrtf(powf((param.X+1.f)/2.f/param.X*pik+(param.X-1.f)/2.f/param.X, 3)));
    }else
    {
        return param.c*powf(pik,(param.X-1.f)/2.f/param.X)/(param.X*(P+param.B));
    }
}

void Speed(float P, Param& p){
    switch (p.type )
    {
    case WAVETYPE::RAZR:
        p.a = (p.X - 1)/2.f/p.X*p.ro*p.c*(1.f-(P+p.B)/(p.p+p.B))/(1.f-powf((P - p.B)/(p.p+p.B),(p.X-1)/2.f/p.X ));
        if (HALF::LEFT == p.course)
        {
            p.D = p.u-p.c;
        }else
        {
            p.D = p.u+p.c;
        }
        break;
    case WAVETYPE::UDAR:
        p.a = sqrtf(p.ro*((p.X+1.f)/2.f*(P+p.B)+(p.X-1.f)/2.f*(p.p+p.B)));
        if (HALF::LEFT == p.course)
        {
            p.D = p.u-p.a/p.ro;
        }else
        {
            p.D = p.u+p.a/p.ro;
        }
        break;
    default:
        break;
    }
} 

float SpeedContact(Param&left, Param&right){
    return (left.a*left.u+right.a*right.u+left.p-right.p)/(left.a-right.a);
}

void Adiabad(float U,float P,Param& p){


}


void Param::Print(){
    if (course == HALF::LEFT)
    {
        std::cout<<"Left ";
    }else
    {
        std::cout<<"Right ";
    }
    if (type == WAVETYPE::UDAR)
    {
        std::cout<<"udar wave:\n";
    }
    if (type == WAVETYPE::RAZR)
    {
        std::cout<<"razr wave:\n";
    }
    std::cout<<"p = "<<p<< std::endl<<"ro = "<<ro<< std::endl<<"u = "<<u<< std::endl<<"X = "<< X << std::endl<<"B = "<<B<< std::endl<<"c = "<< c << std::endl<<"a = "<< a<<std::endl<<"D = "<< D<< std::endl<<"R = "<<R<<std::endl;
    
    
}


void Param::setType(WAVETYPE type){
    this->type=type;
}

int main(){
    
    Param Left(1.f,1.f,10.f,1.667f,0.f, HALF::LEFT);
    Param Right(480.f,8.f,10.f,1.667f,0.f, HALF::RIGHT);
    
    float uud, uraz, uvac;
    uud = Uud(Left,Right);
    uraz = Uraz(Left,Right);
    uvac = Uvac(Left,Right);
    if (uvac < deltU(Left,Right) && deltU(Left,Right) <uud)
    {
        std::cout.precision(8);
        std::cout<<"Left - udarnaya, Right - razryazheniya\n";
        float P;
        Left.setType(WAVETYPE::UDAR);
        Right.setType(WAVETYPE::RAZR);
        P = Pstart(Left, Right);
        std::cout<<"P = "<<P<<"\nStart:\n";
        float delta;
        do
        {
            delta = (Func(P,Left)+Func(P,Right)-deltU(Left,Right))/(dFunc(P, Left)+ dFunc(P,Right));
            P -=delta; 
            std::cout<<"delta = "<<delta<<std::endl;
        } while (fabs(delta) > 0.001f);
        std::cout<<"P = "<<P<<std::endl;
        
        Speed(P,Left);
        Speed(P,Right);
        float U = SpeedContact(Left,Right);
        std::cout<<"U = "<<U<<"\n";
        Left.Print();
        Right.Print();
    }
    return 0;
}