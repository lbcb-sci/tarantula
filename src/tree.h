// Copyright (c) 2021 Cecilia Lee, Robert Vaser

#ifndef TARANTULA_TREE_HPP_
#define TARANTULA_TREE_HPP_


#include <vector>
#include <memory>

#include "vertex.h"

#include <iostream>

using namespace std;

namespace directedforce {

class Box{
public:
	MathVector c1;
	MathVector c2;
	MathVector c3;
	MathVector c4;

	bool in(MathVector pos) {
		if (pos.x >= c1.x && pos.x <= c2.x) {
			if (pos.y >= c1.y && pos.y <= c4.y) {
				return true;
			}
		}
		return false;
	}
};

class Node {
public:

	shared_ptr<Vertex> n; 
	shared_ptr<Node> first; 
	shared_ptr<Node> second; 
	shared_ptr<Node> third ; 
	shared_ptr<Node> fourth; 

	Box box;
	MathVector centreOfMass;
	double mass;

	bool noParticles() {
		if (first == nullptr && second == nullptr && third == nullptr && fourth == nullptr && n == nullptr){
			return true;
		}
		return false;
	}

	int numChild() {
		int sum = 0;
		if (first != nullptr) {
			sum++;
		}
		if (second != nullptr) {
			sum++;
		}
		if (third != nullptr) {
			sum++;
		}
		if (fourth != nullptr) {
			sum++;
		}
		return sum;
	}

	shared_ptr<Node>& getOnlyChild() {
		if (first != nullptr) {
			return first;
		} else if (second != nullptr) {
			return second;
		} else if (third != nullptr) {
			return third;
		} else if (fourth != nullptr) {
			return fourth;
		}
		//return nullptr;
	}

	shared_ptr<Node>& getQuadrant(MathVector pos) {
		//cerr << "before get quadrant" << endl ; 
		double xMidPoint = (box.c2.x - box.c1.x) / 2 + box.c1.x;
		double yMidPoint = (box.c4.y - box.c1.y) / 2 + box.c1.y;
		if (pos.x <= xMidPoint) {
			if (pos.y <= yMidPoint) {
				if (first == nullptr) {
					first = make_shared<Node>(); 
					*first = {nullptr, nullptr, nullptr, nullptr, nullptr, {{box.c1.x, box.c1.y}, {xMidPoint, box.c1.y}, {xMidPoint, yMidPoint}, {box.c1.x, yMidPoint}}};
				}
				return first;
			} else {
				if (fourth == nullptr) {
					fourth = make_shared<Node>(); 
					*fourth = {nullptr, nullptr, nullptr, nullptr, nullptr, {{box.c1.x, yMidPoint}, {xMidPoint, yMidPoint}, {xMidPoint, box.c4.y}, {box.c4.x, box.c4.y}}};
				}
				return fourth;
			}
		} else {
			if (pos.y <= yMidPoint) {
				if (second == nullptr) {
					second = make_shared<Node>(); 
					*second = {nullptr, nullptr, nullptr, nullptr, nullptr, {{xMidPoint, box.c2.y}, {box.c2.x, box.c2.y}, {box.c2.x, yMidPoint}, {xMidPoint, yMidPoint}}};
				}
				return second;
			} else {
				if (third == nullptr) {
					third = make_shared<Node>(); 
					*third = {nullptr, nullptr, nullptr, nullptr, nullptr, {{xMidPoint, yMidPoint}, {box.c3.x, yMidPoint}, {box.c3.x, box.c3.y}, {xMidPoint, box.c3.y}}};
				}
				return third;
			}
		}
	}
};

}  // namespace directedforce

#endif  // TARANTULA_TREE_HPP_
