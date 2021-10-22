#include <iostream>
#include <omp.h>
#include <stdio.h>


using namespace std;

int main()
{
    cout << "Hello World!" << endl;
#pragma omp parallel
{
    int ID = omp_get_thread_num();
    printf(" hello(%d) ", ID);
    printf(" world(%d) \n", ID);
}
    return 0;
}
