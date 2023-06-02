#include <vector>
int fib(int n) {
    int t[n+1];
    t[0] = 0;
    t[1] = 1;
    for (int i = 2; i < n+1; i++) {
        t[i] = t[i-1] + t[i-2];
    }
    return t[n];
}

int fibb(int n) {
    int t0 = 0;
    int t1 = 1;
    
    for (int i = 2; i < n+1; i++) {
        int t = t0 + t1;
        t0 = t1;
        t1 = t;
    }

    return t1;
}
