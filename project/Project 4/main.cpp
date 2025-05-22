#include "shortestP2P.hpp"

using namespace std;

int main () {
	ShortestP2P a;
	a.readGraph();// when you test locally, this function asks you to input the graph and will form a certain data structure.

	int A, B;

	cin >> A;
	do {
		cin >> B;
		a.distance(A, B);
		cin >> A;
	} while (A >= 0);// the loop will run endlessly unless you input a negative integer for "A".
	return 0;
}
