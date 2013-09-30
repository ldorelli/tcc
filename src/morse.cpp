#include <stdio.h>
#include <queue>
#include <iostream>
#include <wchar.h>

using namespace std;

int main (void) {
	wchar_t w;
	wchar_t s[1001];
	wscanf(L"%ls", s);	
	wprintf(L"%ls\n", s);
}
