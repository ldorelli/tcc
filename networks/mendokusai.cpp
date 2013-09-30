#include <stdio.h>

int main (void) {
	int n, m;
	while (scanf ("%d%d", &n, &m) != EOF) {
		printf ("%d %d\n", n-1, m-1);
	}
	return 0;
}
