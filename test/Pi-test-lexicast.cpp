#include"Pi-Lexicast.h"
#include<iostream>
#include<iomanip>
#include"MML3-Timer.h"


int main()
{
    MML3::Timer timer;

    size_t NR=100000;

    timer.start();
    double acc=0.0;
    std::string a("1.3456");
    for(size_t i=0;i!=NR;++i)
    {
        acc+=Pitty::lexicast<double>(a);
    }
    double t=timer.now();
    std::cout << "time " << (t/NR) << "\t   ignore" << std::setprecision(5) << acc;


    return 0;
}
