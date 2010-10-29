#include <cstdio>
#include <cstdlib>
#include <iostream>

const char sigma[] = {'a', 'c', 't', 'g'};

int main(int argc, char *argv[])
{
    int i;
    int bytes = (1 << 20) * atoi(argv[1]);

    for (i = 0 ; i < bytes ; ++i)
        std::cout << sigma[rand() % sizeof(sigma)];
}
