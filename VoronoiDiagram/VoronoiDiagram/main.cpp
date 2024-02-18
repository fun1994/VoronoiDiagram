#include "VoronoiDiagram.h"

void generatePoints(std::vector<Point>& P) {
	P.resize(11);
	P[0].x = 213; P[0].y = -90;
	P[1].x = 164; P[1].y = -165;
	P[2].x = 230; P[2].y = -134;
	P[3].x = 356; P[3].y = -100;
	P[4].x = 244; P[4].y = -247;
	P[5].x = 405; P[5].y = -186;
	P[6].x = 343; P[6].y = -144;
	P[7].x = 222; P[7].y = -325;
	P[8].x = 308; P[8].y = -189;
	P[9].x = 212; P[9].y = -223;
	P[10].x = 351; P[10].y = -248;
}

void test() {
	VoronoiDiagram VD;
	std::vector<DCEL> D(2);
	std::vector<Point> P;
	generatePoints(P);
	VD.incrementalConstruction(D[0], P);
	VD.divideAndConquer(D[1], P);
}

int main() {
	test();
	return 0;
}
