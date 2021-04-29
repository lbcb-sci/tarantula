#include <math.h>

#include <vector>
#include <memory>
#include <iostream>

namespace directedforce {

class MathVector{
  public: 
    MathVector(){}
    MathVector(double x_, double y_){
      x = x_;
      y = y_;  
    }
    ~MathVector(){

    }
    double x; 
    double y; 

    MathVector operator =(const int i){
      MathVector mv; 
      x=i; 
      y=i; 
      return mv; 
    }

    MathVector operator +(const MathVector op){
      MathVector mv; 
      mv.x = x + op.x; 
      mv.y = y + op.y; 
      return mv; 
    }
    MathVector operator -(const MathVector op){
      MathVector mv; 
      mv.x = x - op.x; 
      mv.y = y - op.y; 
      return mv; 
    }

    MathVector operator *(const MathVector op){
      MathVector mv; 
      mv.x = x*op.x; 
      mv.y = y*op.y; 
      return mv; 
    }

    MathVector operator *(const double d){
      MathVector mv; 
      mv.x = x*d; 
      mv.y = y*d; 
      return mv; 
    }
    
    MathVector operator /(const MathVector op){
      MathVector mv; 
      mv.x = x/op.x; 
      mv.y = y/op.y;
      return mv; 
    }
    
    MathVector operator /(const double d){
      MathVector mv; 
      mv.x = x/d; 
      mv.y = y/d; 
      return mv; 
    }

    void operator +=(const MathVector op){
      x += op.x; 
      y += op.y; 
    }

    void operator -=(const MathVector op){
      x -= op.x; 
      y -= op.y; 
    }

    void operator -=(const double d){
      x -= d; 
      y -= d; 
    }


    void operator *=(const double d){
      x *= d; 
      y *= d; 
    }

    void operator /=(const double d){
      x /= d; 
      y /= d; 
    }

    double abs(){
      return sqrt(x*x+y*y); 
    }

    MathVector min(double num){
      MathVector mv; 
      if (num<x)
        mv.x = num; 
      else
        mv.x = x; 
      
      if (num<y)
        mv.y = num; 
      else 
        mv.y = y; 
      
      return mv; 
    }

    MathVector max(double num){
      MathVector mv; 
      if (num>x)
        mv.x = num; 
      else 
        mv.x = x; 
      
      if (num>y)
        mv.y = y; 
      else 
        mv.y = y; 
      
      return mv; 
    }
}; 

struct Vertex{
  MathVector pos; 
  MathVector disp; 
}; 

class Box{
public:
	MathVector c1;
	MathVector c2;
	MathVector c3;
	MathVector c4;

	bool in(MathVector pos){
		if (pos.x >= c1.x && pos.x <= c2.x)
		{
			if (pos.y >= c1.y && pos.y <= c4.y)
			{
				return true;
			}
		}
		return false;
	}
};

class Node{
public:
	std::shared_ptr<Vertex> n; 
	std::shared_ptr<Node> first; 
	std::shared_ptr<Node> second; 
	std::shared_ptr<Node> third ; 
	std::shared_ptr<Node> fourth; 

	Box box;
	MathVector centreOfMass;
	double mass = 1;

	Node(
		std::shared_ptr<Vertex> n,
		std::shared_ptr<Node> first, 
		std::shared_ptr<Node> second, 
		std::shared_ptr<Node> third, 
		std::shared_ptr<Node> fourth, 
		Box box) : 
		n(n), 
		first(first), 
		second(second),
		third(third),
		fourth(fourth), 
		box(box) {}

	Node() {}
	bool noParticles(){
		if (first == nullptr && second == nullptr && third == nullptr && fourth == nullptr && n == nullptr){
			return true;
		}
		return false;
	}

	int numChild(){
		int sum = 0;
		if (first != nullptr){
			sum++;
		}
		if (second != nullptr){
			sum++;
		}
		if (third != nullptr){
			sum++;
		}
		if (fourth != nullptr){
			sum++;
		}
		return sum;
	}

	std::shared_ptr<Node>& getOnlyChild(){
		if (first != nullptr){
			return first;
		}
		else if (second != nullptr){
			return second;
		}
		else if (third != nullptr){
			return third;
		}
		else if (fourth != nullptr){
			return fourth;
		}
	}

	std::shared_ptr<Node>& getQuadrant(MathVector pos)
	{
		double xMidPoint = (box.c2.x - box.c1.x) / 2 + box.c1.x;
		double yMidPoint = (box.c4.y - box.c1.y) / 2 + box.c1.y;
		if (pos.x <= xMidPoint){
			if (pos.y <= yMidPoint){
				if (first == nullptr){
					Box box = {{box.c1.x, box.c1.y}, {xMidPoint, box.c1.y}, {xMidPoint, yMidPoint}, {box.c1.x, yMidPoint}};
					first = std::make_shared<Node>(nullptr, nullptr, nullptr, nullptr, nullptr, box); 
				}
				return first;
			}
			else{
				if (fourth == nullptr){
					Box box = {{box.c1.x, yMidPoint}, {xMidPoint, yMidPoint}, {xMidPoint, box.c4.y}, {box.c4.x, box.c4.y}};
					fourth = std::make_shared<Node>(nullptr, nullptr, nullptr, nullptr, nullptr, box); 
				}
				return fourth;
			}
		}
		else{
			if (pos.y <= yMidPoint) {
				if (second == nullptr) {
					Box box = {{xMidPoint, box.c2.y}, {box.c2.x, box.c2.y}, {box.c2.x, yMidPoint}, {xMidPoint, yMidPoint}};
					second = std::make_shared<Node>(nullptr, nullptr, nullptr, nullptr, nullptr, box); 
				}
				return second;
			} else{ 
				if (third == nullptr){
					//cerr << "new " << endl; 
					Box box = {{xMidPoint, yMidPoint}, {box.c3.x, yMidPoint}, {box.c3.x, box.c3.y}, {xMidPoint, box.c3.y}};
					third = std::make_shared<Node>(nullptr, nullptr, nullptr, nullptr, nullptr, box); 
					//cerr << "after new" << endl;
				}
				return third;
			}
		}
	}
};

}