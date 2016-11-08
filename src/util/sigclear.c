#include <stdlib.h>
#include <signal.h>
void sigclear_(int *signum)
{
    signal(*signum, NULL);
}
