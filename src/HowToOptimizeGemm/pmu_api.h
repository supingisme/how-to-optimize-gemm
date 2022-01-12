#ifndef PMU_API_H
#define PMU_API_H

#include <stdint.h>

typedef struct pmu_paras
{
    uint64_t evt[6];
    uint64_t clocks;
} pmu_paras_t;

double dclock();
int init_pmu();
int start_pmu();
int stop_pmu();
void get_pmu(pmu_paras_t *p);
#endif
