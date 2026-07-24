#include "../include/isa.h"

int main()
{
    for(double z = 0.0; z < 7000.0; z += 500.0){
        double zft = z * M_2_FT;

        printf("%5.0lf, %5.0lf, %.3lf, %.4e, %.4e\n",
                z,
                ISA_H(&zft)*1000.0,
                ISA_T(&zft),
                ISA_P(&zft)/100,
                ISA_Rho(&zft));
    }
                 
    return EXIT_SUCCESS;
}

