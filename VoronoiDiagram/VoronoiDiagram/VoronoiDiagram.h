#pragma once
#include <vector>
#include "Point.h"
#include "DCEL.h"
#define M 1e3
#define EPS 1e-9

class VoronoiDiagram {
	void initialize(DCEL& D) {
		D.indexV = 0;
		Vertex v;
		v.x = -M; v.y = -M; v.incidentEdge = 0; D.V[D.indexV++] = v;
		v.x = M; v.y = -M; v.incidentEdge = 2; D.V[D.indexV++] = v;
		v.x = M; v.y = M; v.incidentEdge = 4; D.V[D.indexV++] = v;
		v.x = -M; v.y = M; v.incidentEdge = 6; D.V[D.indexV++] = v;
		D.indexE = 0;
		HalfEdge e;
		e.origin = 0; e.twin = 1; e.prev = 6; e.next = 2; e.incidentFace = 1; D.E[D.indexE++] = e;
		e.origin = 1; e.twin = 0; e.prev = 3; e.next = 7; e.incidentFace = 0; D.E[D.indexE++] = e;
		e.origin = 1; e.twin = 3; e.prev = 0; e.next = 4; e.incidentFace = 1; D.E[D.indexE++] = e;
		e.origin = 2; e.twin = 2; e.prev = 5; e.next = 1; e.incidentFace = 0; D.E[D.indexE++] = e;
		e.origin = 2; e.twin = 5; e.prev = 2; e.next = 6; e.incidentFace = 1; D.E[D.indexE++] = e;
		e.origin = 3; e.twin = 4; e.prev = 7; e.next = 3; e.incidentFace = 0; D.E[D.indexE++] = e;
		e.origin = 3; e.twin = 7; e.prev = 4; e.next = 0; e.incidentFace = 1; D.E[D.indexE++] = e;
		e.origin = 0; e.twin = 6; e.prev = 1; e.next = 5; e.incidentFace = 0; D.E[D.indexE++] = e;
		D.indexF = 0;
		Face f;
		f.outerComponent = -1; f.innerComponents.insert(1); f.site = -1; D.F[D.indexF++] = f;
		f.outerComponent = 0; f.innerComponents.clear(); f.site = 0; D.F[D.indexF++] = f;
	}
	void insert(DCEL& D, std::vector<Point>& P, int i) {
		int i0 = closest(D, P, i);
		std::vector<int> E;
		double a, b, c;
		bisector(P[i], P[D.F[i0].site], a, b, c);
		int j, k;
		int l = D.F[i0].outerComponent;
		do {
			bool flag1 = toLeft(D.V[D.E[l].origin], a, b, c);
			bool flag2 = toLeft(D.V[D.E[D.E[l].twin].origin], a, b, c);
			if (flag1 && !flag2) {
				j = l;
			}
			else if (!flag1 && flag2) {
				k = l;
			}
			l = D.E[l].next;
		} while (l != D.F[i0].outerComponent);
		int i1 = D.E[D.E[k].twin].incidentFace;
		int i2 = D.E[D.E[j].twin].incidentFace;
		int u0 = intersect(D, j, a, b, c);
		int v0 = intersect(D, k, a, b, c);
		addEdge(D, E, u0, v0);
		int u = v0;
		int v;
		bool bounded = true;
		while (true) {
			if (D.F[i1].site < 0) {
				bounded = false;
				break;
			}
			D.V[u].incidentEdge = D.E[D.E[D.V[u].incidentEdge].prev].twin;
			if (i1 == i2) {
				D.V[u0].incidentEdge = D.E[D.E[D.V[u0].incidentEdge].prev].twin;
				addEdge(D, E, u, u0);
				break;
			}
			bisector(P[i], P[D.F[i1].site], a, b, c);
			l = D.F[i1].outerComponent;
			do {
				bool flag1 = toLeft(D.V[D.E[l].origin], a, b, c);
				bool flag2 = toLeft(D.V[D.E[D.E[l].twin].origin], a, b, c);
				if (!flag1 && flag2) {
					k = l;
				}
				l = D.E[l].next;
			} while (l != D.F[i1].outerComponent);
			i1 = D.E[D.E[k].twin].incidentFace;
			v = intersect(D, k, a, b, c);
			addEdge(D, E, u, v);
			u = v;
		}
		if (!bounded) {
			std::vector<int> E1;
			v = u0;
			while (true) {
				if (D.F[i2].site < 0) {
					break;
				}
				D.V[v].incidentEdge = D.E[D.E[D.V[v].incidentEdge].prev].twin;
				bisector(P[i], P[D.F[i2].site], a, b, c);
				l = D.F[i2].outerComponent;
				do {
					bool flag1 = toLeft(D.V[D.E[l].origin], a, b, c);
					bool flag2 = toLeft(D.V[D.E[D.E[l].twin].origin], a, b, c);
					if (flag1 && !flag2) {
						j = l;
					}
					l = D.E[l].next;
				} while (l != D.F[i2].outerComponent);
				i2 = D.E[D.E[j].twin].incidentFace;
				u = intersect(D, j, a, b, c);
				addEdge(D, E1, u, v);
				v = u;
			}
			E.insert(E.begin(), E1.rbegin(), E1.rend());
		}
		mergeFaces(D, E, i);
	}
	int closest(DCEL& D, std::vector<Point>& P, int i) {
		int j = -1;
		for (std::unordered_map<int, Face>::iterator it = D.F.begin(); it != D.F.end(); it++) {
			if (0 <= it->second.site) {
				if (j < 0 || distance(P[i], P[it->second.site]) < distance(P[i], P[D.F[j].site])) {
					j = it->first;
				}
			}
		}
		return j;
	}
	double distance(Point& p, Point& q) {
		return sqrt(pow(p.x - q.x, 2) + pow(p.y - q.y, 2));
	}
	void bisector(Point& p, Point& q, double& a, double& b, double& c) {
		a = 2 * (p.x - q.x);
		b = 2 * (p.y - q.y);
		c = -pow(p.x, 2) + pow(q.x, 2) - pow(p.y, 2) + pow(q.y, 2);
	}
	bool toLeft(Vertex& v, double a, double b, double c) {
		return a * v.x + b * v.y + c > 0;
	}
	int intersect(DCEL& D, int i, double a, double b, double c) {
		Vertex p = D.V[D.E[i].origin];
		Vertex q = D.V[D.E[D.E[i].twin].origin];
		double t = -(a * q.x + b * q.y + c) / (a * (p.x - q.x) + b * (p.y - q.y));
		Vertex v;
		v.x = t * p.x + (1 - t) * q.x;
		v.y = t * p.y + (1 - t) * q.y;
		v.incidentEdge = D.indexE;
		D.V[D.indexV++] = v;
		HalfEdge e1, e2;
		e1.origin = D.indexV - 1;
		e2.origin = D.indexV - 1;
		D.E[D.E[i].twin].twin = D.indexE + 1;
		e2.twin = D.E[i].twin;
		D.E[i].twin = D.indexE;
		e1.twin = i;
		D.E[D.E[e2.twin].next].prev = D.indexE;
		e1.next = D.E[e2.twin].next;
		D.E[D.E[e1.twin].next].prev = D.indexE + 1;
		e2.next = D.E[e1.twin].next;
		D.E[e2.twin].next = D.indexE;
		e1.prev = e2.twin;
		D.E[e1.twin].next = D.indexE + 1;
		e2.prev = e1.twin;
		e1.incidentFace = D.E[e2.twin].incidentFace;
		e2.incidentFace = D.E[e1.twin].incidentFace;
		D.E[D.indexE++] = e1;
		D.E[D.indexE++] = e2;
		return D.indexV - 1;
	}
	void addEdge(DCEL& D, std::vector<int>& E, int u, int v) {
		E.push_back(D.indexE);
		HalfEdge e1, e2;
		e1.origin = u;
		e2.origin = v;
		e1.twin = D.indexE + 1;
		e2.twin = D.indexE;
		D.E[D.E[D.E[D.V[u].incidentEdge].twin].next].prev = D.indexE + 1;
		e2.next = D.E[D.E[D.V[u].incidentEdge].twin].next;
		D.E[D.E[D.V[u].incidentEdge].twin].next = D.indexE;
		e1.prev = D.E[D.V[u].incidentEdge].twin;
		D.E[D.E[D.E[D.V[v].incidentEdge].twin].next].prev = D.indexE;
		e1.next = D.E[D.E[D.V[v].incidentEdge].twin].next;
		D.E[D.E[D.V[v].incidentEdge].twin].next = D.indexE + 1;
		e2.prev = D.E[D.V[v].incidentEdge].twin;
		e2.incidentFace = D.E[D.E[D.V[u].incidentEdge].twin].incidentFace;
		D.E[D.indexE++] = e1;
		D.E[D.indexE++] = e2;
		int i = D.indexE - 2;
		do {
			D.E[i].incidentFace = D.indexF;
			i = D.E[i].next;
		} while (i != D.indexE - 2);
		Face f;
		D.F[e2.incidentFace].outerComponent = D.indexE - 1;
		f.outerComponent = D.indexE - 2;
		D.F[D.indexF++] = f;
	}
	void mergeFaces(DCEL& D, std::vector<int>& E, int i) {
		int size = E.size();
		for (int j = size - 1; 0 <= j; j--) {
			int k = D.E[E[j]].next;
			do {
				double x1 = D.V[D.E[k].origin].x;
				double y1 = D.V[D.E[k].origin].y;
				double x2 = D.V[D.E[D.E[k].next].origin].x;
				double y2 = D.V[D.E[D.E[k].next].origin].y;
				if ((abs(x1 - M) < EPS && abs(x2 - M) < EPS) || (abs(x1 + M) < EPS && abs(x2 + M) < EPS) || (abs(y1 - M) < EPS && abs(y2 - M) < EPS) || (abs(y1 + M) < EPS && abs(y2 + M) < EPS)) {
					E.push_back(k);
				}
				k = D.E[k].next;
			} while (k != E[j]);
		}
		for (int j = 0; j < size; j++) {
			int k = D.E[D.E[E[j]].next].next;
			do {
				double x = D.V[D.E[k].origin].x;
				double y = D.V[D.E[k].origin].y;
				if (abs(abs(x) - M) > EPS && abs(abs(y) - M) > EPS) {
					D.V.erase(D.E[k].origin);
				}
				k = D.E[k].next;
			} while (k != E[j]);
		}
		for (int j = 0; j < size; j++) {
			int k = D.E[E[j]].next;
			do {
				int l = D.E[k].next;
				if (D.V.find(D.E[k].origin) != D.V.end() && D.V.find(D.E[l].origin) != D.V.end()) {
					double x1 = D.V[D.E[k].origin].x;
					double y1 = D.V[D.E[k].origin].y;
					double x2 = D.V[D.E[l].origin].x;
					double y2 = D.V[D.E[l].origin].y;
					if ((abs(x1 - M) > EPS || abs(x2 - M) > EPS) && (abs(x1 + M) > EPS || abs(x2 + M) > EPS) && (abs(y1 - M) > EPS || abs(y2 - M) > EPS) && (abs(y1 + M) > EPS || abs(y2 + M) > EPS)) {
						D.E.erase(k);
					}
				}
				else {
					D.E.erase(k);
				}
				k = l;
			} while (k != E[j]);
		}
		for (int j = 1; j < size; j++) {
			D.F.erase(D.E[E[j]].incidentFace);
		}
		for (int j = 0; j < E.size(); j++) {
			D.V[D.E[E[j]].origin].incidentEdge = E[j];
		}
		for (int j = 0; j < E.size(); j++) {
			D.E[E[j]].next = E[(j + 1) % E.size()];
			D.E[E[(j + 1) % E.size()]].prev = E[j];
		}
		for (int j = 1; j < E.size(); j++) {
			D.E[E[j]].incidentFace = D.E[E[0]].incidentFace;
		}
		D.F[D.E[E[0]].incidentFace].site = i;
		int j = size;
		while (j < E.size()) {
			int k = j + 1;
			while (k < E.size()) {
				double x = D.V[D.E[E[k]].origin].x;
				double y = D.V[D.E[E[k]].origin].y;
				if (abs(abs(x) - M) < EPS && abs(abs(y) - M) < EPS) {
					break;
				}
				k++;
			}
			for (int l = j + 1; l < k; l++) {
				D.V.erase(D.E[E[l]].origin);
			}
			D.E[E[j]].next = E[k % E.size()];
			D.E[E[k % E.size()]].prev = E[j];
			D.E[D.E[E[k - 1]].twin].next = D.E[D.E[E[j]].twin].next;
			D.E[D.E[D.E[E[j]].twin].next].prev = D.E[E[k - 1]].twin;
			for (int l = j; l < k - 1; l++) {
				D.E.erase(D.E[E[l]].twin);
			}
			D.E[E[j]].twin = D.E[E[k - 1]].twin;
			D.E[D.E[E[k - 1]].twin].twin = E[j];
			for (int l = j + 1; l < k; l++) {
				D.E.erase(E[l]);
			}
			j = k;
		}
	}
	int partition(std::vector<Point>& P, std::vector<int>& Q, int low, int high) {
		int pivot = Q[low];
		while (low < high) {
			while (low < high && P[Q[high]].x >= P[pivot].x) {
				high--;
			}
			Q[low] = Q[high];
			while (low < high && P[Q[low]].x <= P[pivot].x) {
				low++;
			}
			Q[high] = Q[low];
		}
		Q[low] = pivot;
		return low;
	}
	void quickSort(std::vector<Point>& P, std::vector<int>& Q, int low, int high) {
		if (low < high) {
			int pivot = partition(P, Q, low, high);
			quickSort(P, Q, low, pivot - 1);
			quickSort(P, Q, pivot + 1, high);
		}
	}
	void divideAndConquer(DCEL& D, std::vector<int>& CH, int& rightmost, std::vector<Point>& P, std::vector<int>& Q, int left, int right) {
		if (left == right) {
			trivial(D, Q[left]);
			CH.push_back(Q[left]);
			rightmost = 0;
		}
		else {
			int mid = (left + right) / 2;
			DCEL DL, DR;
			std::vector<int> CHL, CHR;
			int rightmostL, rightmostR;
			divideAndConquer(DL, CHL, rightmostL, P, Q, left, mid);
			divideAndConquer(DR, CHR, rightmostR, P, Q, mid + 1, right);
			merge(D, CH, rightmost, DL, DR, CHL, CHR, rightmostL, rightmostR, P);
		}
	}
	void trivial(DCEL& D, int i) {
		D.indexV = 0;
		Vertex v;
		v.x = -M; v.y = -M; v.incidentEdge = 0; D.V[D.indexV++] = v;
		v.x = M; v.y = -M; v.incidentEdge = 2; D.V[D.indexV++] = v;
		v.x = M; v.y = M; v.incidentEdge = 4; D.V[D.indexV++] = v;
		v.x = -M; v.y = M; v.incidentEdge = 6; D.V[D.indexV++] = v;
		D.indexE = 0;
		HalfEdge e;
		e.origin = 0; e.twin = 1; e.prev = 6; e.next = 2; e.incidentFace = i; D.E[D.indexE++] = e;
		e.origin = 1; e.twin = 0; e.prev = 3; e.next = 7; e.incidentFace = -1; D.E[D.indexE++] = e;
		e.origin = 1; e.twin = 3; e.prev = 0; e.next = 4; e.incidentFace = i; D.E[D.indexE++] = e;
		e.origin = 2; e.twin = 2; e.prev = 5; e.next = 1; e.incidentFace = -1; D.E[D.indexE++] = e;
		e.origin = 2; e.twin = 5; e.prev = 2; e.next = 6; e.incidentFace = i; D.E[D.indexE++] = e;
		e.origin = 3; e.twin = 4; e.prev = 7; e.next = 3; e.incidentFace = -1; D.E[D.indexE++] = e;
		e.origin = 3; e.twin = 7; e.prev = 4; e.next = 0; e.incidentFace = i; D.E[D.indexE++] = e;
		e.origin = 0; e.twin = 6; e.prev = 1; e.next = 5; e.incidentFace = -1; D.E[D.indexE++] = e;
		Face f;
		f.outerComponent = -1; f.innerComponents.insert(1); f.site = -1; D.F[-1] = f;
		f.outerComponent = 0; f.innerComponents.clear(); f.site = i; D.F[i] = f;
	}
	void merge(DCEL& D, std::vector<int>& CH, int& rightmost, DCEL& DL, DCEL& DR, std::vector<int>& CHL, std::vector<int>& CHR, int rightmostL, int rightmostR, std::vector<Point>& P) {
		int sL, tL, sR, tR;
		stitch(CH, rightmost, CHL, CHR, rightmostL, rightmostR, sL, tL, sR, tR, P);
		std::vector<int> VL1, VR1, EL1, ER1;
		std::unordered_map<int, int> VL2, VR2, EL2, ER2;
		computeContourBetween(VL1, VR1, EL1, ER1, VL2, VR2, EL2, ER2, DL, DR, CHL, CHR, sL, tL, sR, tR, P);
		clip(DL, VL1, EL1, VL2, EL2, 'l');
		clip(DR, VR1, ER1, VR2, ER2, 'r');
		stitch(D, DL, DR, VL1, VR1, EL1, ER1, VL2, VR2, EL2, ER2);
	}
	void stitch(std::vector<int>& CH, int& rightmost, std::vector<int>& CHL, std::vector<int>& CHR, int rightmostL, int rightmostR, int& sL, int& tL, int& sR, int& tR, std::vector<Point>& P) {
		if (CHL.size() == 1 && CHR.size() == 1) {
			CH.push_back(CHL[0]);
			CH.push_back(CHR[0]);
			rightmost = 1;
			sL = 0;
			tL = 0;
			sR = 0;
			tR = 0;
		}
		else if (CHL.size() == 2 && CHR.size() == 1) {
			if (toLeft(P[CHL[0]], P[CHL[1]], P[CHR[0]])) {
				CH.push_back(CHL[0]);
				CH.push_back(CHL[1]);
				CH.push_back(CHR[0]);
				rightmost = 2;
				sL = 0;
				tL = 1;
				sR = 0;
				tR = 0;
			}
			else {
				CH.push_back(CHL[0]);
				CH.push_back(CHR[0]);
				CH.push_back(CHL[1]);
				rightmost = 1;
				sL = 1;
				tL = 0;
				sR = 0;
				tR = 0;
			}
		}
		else if (CHL.size() == 2 && CHR.size() == 2) {
			if (toLeft(P[CHL[1]], P[CHR[0]], P[CHR[1]])) {
				if (toLeft(P[CHR[1]], P[CHL[1]], P[CHL[0]])) {
					sL = 1;
					tR = 1;
				}
				else {
					sL = 0;
					tR = 1;
				}
			}
			else {
				if (toLeft(P[CHR[0]], P[CHL[1]], P[CHL[0]])) {
					sL = 1;
					tR = 0;
				}
				else {
					if (toLeft(P[CHL[0]], P[CHR[0]], P[CHR[1]])) {
						sL = 0;
						tR = 1;
					}
					else {
						sL = 0;
						tR = 0;
					}
				}
			}
			if (toLeft(P[CHL[1]], P[CHR[0]], P[CHR[1]])) {
				if (toLeft(P[CHR[0]], P[CHL[1]], P[CHL[0]])) {
					if (toLeft(P[CHL[0]], P[CHR[0]], P[CHR[1]])) {
						tL = 0;
						sR = 0;
					}
					else {
						tL = 0;
						sR = 1;
					}
				}
				else {
					tL = 1;
					sR = 0;
				}
			}
			else {
				if (toLeft(P[CHR[1]], P[CHL[1]], P[CHL[0]])) {
					tL = 0;
					sR = 1;
				}
				else {
					tL = 1;
					sR = 1;
				}
			}
			if (sL == 0 && tL == 0 && sR == 0 && tR == 1) {
				CH.push_back(CHL[0]);
				CH.push_back(CHR[0]);
				CH.push_back(CHR[1]);
				rightmost = 2;
			}
			else if (sL == 0 && tL == 0 && sR == 1 && tR == 0) {
				CH.push_back(CHL[0]);
				CH.push_back(CHR[1]);
				CH.push_back(CHR[0]);
				rightmost = 1;
			}
			else if (sL == 0 && tL == 1 && sR == 0 && tR == 1) {
				CH.push_back(CHL[0]);
				CH.push_back(CHL[1]);
				CH.push_back(CHR[0]);
				CH.push_back(CHR[1]);
				rightmost = 3;
			}
			else if (sL == 0 && tL == 1 && sR == 1 && tR == 0) {
				CH.push_back(CHL[0]);
				CH.push_back(CHL[1]);
				CH.push_back(CHR[1]);
				CH.push_back(CHR[0]);
				rightmost = 2;
			}
			else if (sL == 0 && tL == 1 && sR == 1 && tR == 1) {
				CH.push_back(CHL[0]);
				CH.push_back(CHL[1]);
				CH.push_back(CHR[1]);
				rightmost = 2;
			}
			else if (sL == 1 && tL == 0 && sR == 0 && tR == 1) {
				CH.push_back(CHL[0]);
				CH.push_back(CHR[0]);
				CH.push_back(CHR[1]);
				CH.push_back(CHL[1]);
				rightmost = 2;
			}
			else if (sL == 1 && tL == 0 && sR == 1 && tR == 0) {
				CH.push_back(CHL[0]);
				CH.push_back(CHR[1]);
				CH.push_back(CHR[0]);
				CH.push_back(CHL[1]);
				rightmost = 1;
			}
			else if (sL == 1 && tL == 0 && sR == 1 && tR == 1) {
				CH.push_back(CHL[0]);
				CH.push_back(CHR[1]);
				CH.push_back(CHL[1]);
				rightmost = 1;
			}
		}
		else if (CHL.size() > 2 && CHR.size() == 2) {
			sL = rightmostL;
			tR = 0;
			while (patternOfTurn(P, CHL, P[CHR[tR]], sL) != 's') {
				sL == CHL.size() - 1 ? sL = 0 : sL++;
			}
			if (toLeft(P[CHL[sL]], P[CHR[0]], P[CHR[1]])) {
				tR++;
				while (patternOfTurn(P, CHL, P[CHR[tR]], sL) != 's') {
					sL == CHL.size() - 1 ? sL = 0 : sL++;
				}
			}
			tL = rightmostL;
			sR = 0;
			while (patternOfTurn(P, CHL, P[CHR[sR]], tL) != 't') {
				tL == 0 ? tL = CHL.size() - 1 : tL--;
			}
			if (!toLeft(P[CHL[tL]], P[CHR[0]], P[CHR[1]])) {
				sR++;
				while (patternOfTurn(P, CHL, P[CHR[sR]], tL) != 't') {
					tL == 0 ? tL = CHL.size() - 1 : tL--;
				}
			}
			CH.insert(CH.end(), CHL.begin(), CHL.begin() + (tL + 1));
			if (sR == 0 && tR == 1) {
				CH.push_back(CHR[0]);
				CH.push_back(CHR[1]);
				rightmost = tL + 2;
			}
			else if (sR == 1 && tR == 0) {
				CH.push_back(CHR[1]);
				CH.push_back(CHR[0]);
				rightmost = tL + 1;
			}
			else if (sR == 1 && tR == 1) {
				CH.push_back(CHR[1]);
				rightmost = tL + 1;
			}
			if (sL > 0) {
				CH.insert(CH.end(), CHL.begin() + sL, CHL.end());
			}
		}
		else if (CHL.size() > 2 && CHR.size() > 2) {
			commonTangent(CHL, CHR, rightmostL, sL, tL, sR, tR, P);
			CH.insert(CH.end(), CHL.begin(), CHL.begin() + (tL + 1));
			CH.insert(CH.end(), CHR.begin() + sR, tR == 0 ? CHR.end() : CHR.begin() + (tR + 1));
			if (tR == 0) {
				CH.push_back(CHR[0]);
			}
			if (sL > 0) {
				CH.insert(CH.end(), CHL.begin() + sL, CHL.end());
			}
			rightmost = rightmostR - sR + tL + 1;
		}
	}
	double area2(Point& p, Point& q, Point& r) {
		return p.x * q.y - p.y * q.x + q.x * r.y - q.y * r.x + r.x * p.y - r.y * p.x;
	}
	bool toLeft(Point& p, Point& q, Point& r) {
		return area2(p, q, r) > 0;
	}
	char patternOfTurn(std::vector<Point>& P, std::vector<int>& CH, Point& x, int i) {
		if (toLeft(x, P[CH[i]], P[CH[i == 0 ? CH.size() - 1 : i - 1]]) && toLeft(x, P[CH[i]], P[CH[i == CH.size() - 1 ? 0 : i + 1]])) {
			return 's';
		}
		else if (!toLeft(x, P[CH[i]], P[CH[i == 0 ? CH.size() - 1 : i - 1]]) && !toLeft(x, P[CH[i]], P[CH[i == CH.size() - 1 ? 0 : i + 1]])) {
			return 't';
		}
		else {
			return 'v';
		}
	}
	void commonTangent(std::vector<int>& CHL, std::vector<int>& CHR, int rightmostL, int& sL, int& tL, int& sR, int& tR, std::vector<Point>& P) {
		sL = rightmostL;
		tR = 0;
		bool flag = true;
		while (true) {
			if (flag) {
				if (patternOfTurn(P, CHR, P[CHL[sL]], tR) == 't') {
					if (patternOfTurn(P, CHL, P[CHR[tR]], sL) == 's') {
						break;
					}
					sL == CHL.size() - 1 ? sL = 0 : sL++;
					flag = false;
				}
				else {
					tR == 0 ? tR = CHR.size() - 1 : tR--;
				}
			}
			else {
				if (patternOfTurn(P, CHL, P[CHR[tR]], sL) == 's') {
					if (patternOfTurn(P, CHR, P[CHL[sL]], tR) == 't') {
						break;
					}
					tR == 0 ? tR = CHR.size() - 1 : tR--;
					flag = true;
				}
				else {
					sL == CHL.size() - 1 ? sL = 0 : sL++;
				}
			}
		}
		tL = rightmostL;
		sR = 0;
		flag = true;
		while (true) {
			if (flag) {
				if (patternOfTurn(P, CHR, P[CHL[tL]], sR) == 's') {
					if (patternOfTurn(P, CHL, P[CHR[sR]], tL) == 't') {
						break;
					}
					tL == 0 ? tL = CHL.size() - 1 : tL--;
					flag = false;
				}
				else {
					sR == CHR.size() - 1 ? sR = 0 : sR++;
				}
			}
			else {
				if (patternOfTurn(P, CHL, P[CHR[sR]], tL) == 't') {
					if (patternOfTurn(P, CHR, P[CHL[tL]], sR) == 's') {
						break;
					}
					sR == CHR.size() - 1 ? sR = 0 : sR++;
					flag = true;
				}
				else {
					tL == 0 ? tL = CHL.size() - 1 : tL--;
				}
			}
		}
	}
	void computeContourBetween(std::vector<int>& VL1, std::vector<int>& VR1, std::vector<int>& EL1, std::vector<int>& ER1, std::unordered_map<int, int>& VL2, std::unordered_map<int, int>& VR2, std::unordered_map<int, int>& EL2, std::unordered_map<int, int>& ER2, DCEL& DL, DCEL& DR, std::vector<int>& CHL, std::vector<int>& CHR, int sL, int tL, int sR, int tR, std::vector<Point>& P) {
		int l = CHL[sL];
		int r = CHR[tR];
		double a, b, c;
		bisector(P[l], P[r], a, b, c);
		while (true) {
			bool flag1 = toLeft(DL.V[DL.E[DL.F[l].outerComponent].origin], a, b, c);
			bool flag2 = toLeft(DL.V[DL.E[DL.E[DL.F[l].outerComponent].twin].origin], a, b, c);
			if (!flag1 && flag2) {
				break;
			}
			DL.F[l].outerComponent = DL.E[DL.F[l].outerComponent].prev;
		}
		while (true) {
			bool flag1 = toLeft(DR.V[DR.E[DR.F[r].outerComponent].origin], a, b, c);
			bool flag2 = toLeft(DR.V[DR.E[DR.E[DR.F[r].outerComponent].twin].origin], a, b, c);
			if (!flag1 && flag2) {
				break;
			}
			DR.F[r].outerComponent = DR.E[DR.F[r].outerComponent].next;
		}
		int upperL = intersect(DL, DL.F[l].outerComponent, a, b, c);
		int upperR = intersect(DR, DR.F[r].outerComponent, a, b, c);
		VL1.push_back(upperL);
		VR1.push_back(upperR);
		VL2[upperL] = upperR;
		VR2[upperR] = upperL;
		int lowerL, lowerR;
		while (true) {
			while (true) {
				bool flag1 = toLeft(DL.V[DL.E[DL.F[l].outerComponent].origin], a, b, c);
				bool flag2 = toLeft(DL.V[DL.E[DL.E[DL.F[l].outerComponent].twin].origin], a, b, c);
				if (flag1 && !flag2) {
					break;
				}
				DL.F[l].outerComponent = DL.E[DL.F[l].outerComponent].prev;
			}
			while (true) {
				bool flag1 = toLeft(DR.V[DR.E[DR.F[r].outerComponent].origin], a, b, c);
				bool flag2 = toLeft(DR.V[DR.E[DR.E[DR.F[r].outerComponent].twin].origin], a, b, c);
				if (flag1 && !flag2) {
					break;
				}
				DR.F[r].outerComponent = DR.E[DR.F[r].outerComponent].next;
			}
			if (l == CHL[tL] && r == CHR[sR]) {
				lowerL = intersect(DL, DL.F[l].outerComponent, a, b, c);
				addEdge(DL, lowerL, upperL, l);
				lowerR = intersect(DR, DR.F[r].outerComponent, a, b, c);
				DR.F[r].outerComponent = DR.E[DR.F[r].outerComponent].next;
				addEdge(DR, lowerR, upperR, r);
				VL1.push_back(lowerL);
				VR1.push_back(lowerR);
				VL2[lowerL] = lowerR;
				VR2[lowerR] = lowerL;
				EL1.push_back(DL.indexE - 2);
				ER1.push_back(DR.indexE - 1);
				EL2[DL.indexE - 2] = DR.indexE - 1;
				ER2[DR.indexE - 1] = DL.indexE - 2;
				break;
			}
			Vertex pL = DL.V[DL.E[DL.F[l].outerComponent].origin];
			Vertex qL = DL.V[DL.E[DL.E[DL.F[l].outerComponent].twin].origin];
			double tL = -(a * qL.x + b * qL.y + c) / (a * (pL.x - qL.x) + b * (pL.y - qL.y));
			double yL = tL * pL.y + (1 - tL) * qL.y;
			Vertex pR = DR.V[DR.E[DR.F[r].outerComponent].origin];
			Vertex qR = DR.V[DR.E[DR.E[DR.F[r].outerComponent].twin].origin];
			double tR = -(a * qR.x + b * qR.y + c) / (a * (pR.x - qR.x) + b * (pR.y - qR.y));
			double yR = tR * pR.y + (1 - tR) * qR.y;
			if (yL < yR) {
				lowerR = intersect(DR, DR.F[r].outerComponent, a, b, c);
				DR.F[r].outerComponent = DR.E[DR.F[r].outerComponent].next;
				addEdge(DR, lowerR, upperR, r);
				lowerL = addEdge(DL, DR.V[lowerR], upperL, l);
				VL1.push_back(lowerL);
				VR1.push_back(lowerR);
				VL2[lowerL] = lowerR;
				VR2[lowerR] = lowerL;
				EL1.push_back(DL.indexE - 2);
				ER1.push_back(DR.indexE - 1);
				EL2[DL.indexE - 2] = DR.indexE - 1;
				ER2[DR.indexE - 1] = DL.indexE - 2;
				r = DR.E[DR.E[DR.F[r].outerComponent].twin].incidentFace;
			}
			else {
				lowerL = intersect(DL, DL.F[l].outerComponent, a, b, c);
				addEdge(DL, lowerL, upperL, l);
				lowerR = addEdge(DR, DL.V[lowerL], upperR, r);
				VL1.push_back(lowerL);
				VR1.push_back(lowerR);
				VL2[lowerL] = lowerR;
				VR2[lowerR] = lowerL;
				EL1.push_back(DL.indexE - 2);
				ER1.push_back(DR.indexE - 1);
				EL2[DL.indexE - 2] = DR.indexE - 1;
				ER2[DR.indexE - 1] = DL.indexE - 2;
				l = DL.E[DL.E[DL.F[l].outerComponent].twin].incidentFace;
			}
			bisector(P[l], P[r], a, b, c);
			upperL = lowerL;
			upperR = lowerR;
		}
	}
	void addEdge(DCEL& D, int lower, int upper, int f) {
		HalfEdge e1, e2;
		e1.origin = lower;
		e2.origin = upper;
		e1.twin = D.indexE + 1;
		e2.twin = D.indexE;
		D.E[D.E[D.E[D.V[lower].incidentEdge].twin].next].prev = D.indexE + 1;
		e2.next = D.E[D.E[D.V[lower].incidentEdge].twin].next;
		D.E[D.E[D.V[lower].incidentEdge].twin].next = D.indexE;
		e1.prev = D.E[D.V[lower].incidentEdge].twin;
		D.E[D.E[D.E[D.V[upper].incidentEdge].twin].next].prev = D.indexE;
		e1.next = D.E[D.E[D.V[upper].incidentEdge].twin].next;
		D.E[D.E[D.V[upper].incidentEdge].twin].next = D.indexE + 1;
		e2.prev = D.E[D.V[upper].incidentEdge].twin;
		e1.incidentFace = f;
		e2.incidentFace = f;
		D.E[D.indexE++] = e1;
		D.E[D.indexE++] = e2;
		D.V[lower].incidentEdge = D.E[D.E[D.V[lower].incidentEdge].prev].twin;
	}
	int addEdge(DCEL& D, Vertex& v, int upper, int f) {
		Vertex u;
		u.x = v.x;
		u.y = v.y;
		u.incidentEdge = D.indexE;
		D.V[D.indexV++] = u;
		HalfEdge e1, e2;
		e1.origin = D.indexV - 1;
		e2.origin = upper;
		e1.twin = D.indexE + 1;
		e2.twin = D.indexE;
		D.E[D.E[D.E[D.V[upper].incidentEdge].twin].next].prev = D.indexE;
		e1.next = D.E[D.E[D.V[upper].incidentEdge].twin].next;
		D.E[D.E[D.V[upper].incidentEdge].twin].next = D.indexE + 1;
		e2.prev = D.E[D.V[upper].incidentEdge].twin;
		e1.prev = D.indexE + 1;
		e2.next = D.indexE;
		e1.incidentFace = f;
		e2.incidentFace = f;
		D.E[D.indexE++] = e1;
		D.E[D.indexE++] = e2;
		return D.indexV - 1;
	}
	void clip(DCEL& D, std::vector<int>& V1, std::vector<int>& E1, std::unordered_map<int, int>& V2, std::unordered_map<int, int>& E2, char lr) {
		for (std::unordered_map<int, int>::iterator it = E2.begin(); it != E2.end(); it++) {
			int i = D.E[it->first].twin;
			do {
				if (D.E.find(i) == D.E.end()) {
					break;
				}
				if (V2.find(D.E[i].origin) == V2.end()) {
					D.V.erase(D.E[i].origin);
				}
				i = D.E[i].next;
			} while (i != D.E[it->first].twin);
			i = D.E[it->first].twin;
			do {
				if (D.E.find(i) == D.E.end()) {
					break;
				}
				int j = D.E[i].next;
				if (E2.find(D.E[i].twin) == E2.end()) {
					D.E.erase(i);
				}
				i = j;
			} while (i != D.E[it->first].twin);
			D.E[D.E[it->first].twin].incidentFace = -1;
		}
		int i = *D.F[-1].innerComponents.begin();
		do {
			int j = D.E[i].next;
			if (D.E.find(D.E[i].twin) == D.E.end()) {
				D.E.erase(i);
			}
			i = j;
		} while (i != *D.F[-1].innerComponents.begin());
		D.F[-1].innerComponents.clear();
		if (lr == 'l') {
			D.F[-1].innerComponents.insert(D.E[D.E[E1[0]].next].twin);
			D.V[V1[0]].incidentEdge = D.E[E1[0]].next;
			for (int i = 1; i < V1.size(); i++) {
				D.V[V1[i]].incidentEdge = E1[i - 1];
			}
			D.E[D.E[E1[0]].twin].prev = D.E[D.E[E1[0]].next].twin;
			D.E[D.E[D.E[E1[0]].next].twin].next = D.E[E1[0]].twin;
			for (int i = 0; i < E1.size() - 1; i++) {
				D.E[D.E[E1[i]].twin].next = D.E[E1[i + 1]].twin;
				D.E[D.E[E1[i + 1]].twin].prev = D.E[E1[i]].twin;
			}
			D.E[D.E[E1[E1.size() - 1]].twin].next = D.E[D.E[E1[E1.size() - 1]].prev].twin;
			D.E[D.E[D.E[E1[E1.size() - 1]].prev].twin].prev = D.E[E1[E1.size() - 1]].twin;
		}
		else if (lr == 'r') {
			D.F[-1].innerComponents.insert(D.E[D.E[E1[0]].prev].twin);
			for (int i = 0; i < V1.size() - 1; i++) {
				D.V[V1[i]].incidentEdge = E1[i];
			}
			D.V[V1[V1.size() - 1]].incidentEdge = D.E[E1[E1.size() - 1]].next;
			D.E[D.E[E1[0]].twin].next = D.E[D.E[E1[0]].prev].twin;
			D.E[D.E[D.E[E1[0]].prev].twin].prev = D.E[E1[0]].twin;
			for (int i = 0; i < E1.size() - 1; i++) {
				D.E[D.E[E1[i]].twin].prev = D.E[E1[i + 1]].twin;
				D.E[D.E[E1[i + 1]].twin].next = D.E[E1[i]].twin;
			}
			D.E[D.E[E1[E1.size() - 1]].twin].prev = D.E[D.E[E1[E1.size() - 1]].next].twin;
			D.E[D.E[D.E[E1[E1.size() - 1]].next].twin].next = D.E[E1[E1.size() - 1]].twin;
		}
	}
	void stitch(DCEL& D, DCEL& DL, DCEL& DR, std::vector<int>& VL1, std::vector<int>& VR1, std::vector<int>& EL1, std::vector<int>& ER1, std::unordered_map<int, int>& VL2, std::unordered_map<int, int>& VR2, std::unordered_map<int, int>& EL2, std::unordered_map<int, int>& ER2) {
		D.indexV = DL.indexV + DR.indexV;
		for (std::unordered_map<int, Vertex>::iterator it = DL.V.begin(); it != DL.V.end(); it++) {
			D.V[it->first] = it->second;
		}
		for (std::unordered_map<int, Vertex>::iterator it = DR.V.begin(); it != DR.V.end(); it++) {
			if (VR2.find(it->first) == VR2.end()) {
				D.V[DL.indexV + it->first] = it->second;
				D.V[DL.indexV + it->first].incidentEdge += DL.indexE;
			}
		}
		D.indexE = DL.indexE + DR.indexE;
		for (std::unordered_map<int, HalfEdge>::iterator it = DL.E.begin(); it != DL.E.end(); it++) {
			if (EL2.find(it->second.twin) == EL2.end()) {
				D.E[it->first] = it->second;
				if (EL2.find(it->first) != EL2.end()) {
					D.E[it->first].twin = DL.indexE + EL2[it->first];
				}
			}
		}
		for (std::unordered_map<int, HalfEdge>::iterator it = DR.E.begin(); it != DR.E.end(); it++) {
			if (ER2.find(it->second.twin) == ER2.end()) {
				D.E[DL.indexE + it->first] = it->second;
				if (VR2.find(it->second.origin) != VR2.end()) {
					D.E[DL.indexE + it->first].origin = VR2[it->second.origin];
				}
				else {
					D.E[DL.indexE + it->first].origin += DL.indexV;
				}
				if (ER2.find(it->first) != ER2.end()) {
					D.E[DL.indexE + it->first].twin = ER2[it->first];
				}
				else {
					D.E[DL.indexE + it->first].twin += DL.indexE;
				}
				if (it->first != DR.E[DR.E[ER1[0]].prev].twin) {
					D.E[DL.indexE + it->first].prev += DL.indexE;
				}
				if (it->first != DR.E[DR.E[ER1[ER1.size() - 1]].next].twin) {
					D.E[DL.indexE + it->first].next += DL.indexE;
				}
			}
		}
		D.E[D.E[D.E[EL1[0]].next].twin].next = D.E[D.E[DL.indexE + ER1[0]].prev].twin;
		D.E[D.E[D.E[DL.indexE + ER1[0]].prev].twin].prev = D.E[D.E[EL1[0]].next].twin;
		D.E[D.E[D.E[EL1[EL1.size() - 1]].prev].twin].prev = D.E[D.E[DL.indexE + ER1[ER1.size() - 1]].next].twin;
		D.E[D.E[D.E[DL.indexE + ER1[ER1.size() - 1]].next].twin].next = D.E[D.E[EL1[EL1.size() - 1]].prev].twin;
		for (std::unordered_map<int, Face>::iterator it = DL.F.begin(); it != DL.F.end(); it++) {
			D.F[it->first] = it->second;
		}
		for (std::unordered_map<int, Face>::iterator it = DR.F.begin(); it != DR.F.end(); it++) {
			if (0 <= it->first) {
				D.F[it->first] = it->second;
				D.F[it->first].outerComponent += DL.indexE;
			}
		}
	}
public:
	void incrementalConstruction(DCEL& D, std::vector<Point>& P) {
		initialize(D);
		for (int i = 1; i < P.size(); i++) {
			insert(D, P, i);
		}
	}
	void divideAndConquer(DCEL& D, std::vector<Point>& P) {
		std::vector<int> Q;
		for (int i = 0; i < P.size(); i++) {
			Q.push_back(i);
		}
		quickSort(P, Q, 0, Q.size() - 1);
		std::vector<int> CH;
		int rightmost;
		divideAndConquer(D, CH, rightmost, P, Q, 0, Q.size() - 1);
	}
};
