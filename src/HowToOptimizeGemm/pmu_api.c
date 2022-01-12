#include <stdio.h>
#include <sys/procmgr.h>
#include <sys/neutrino.h>
#include <stdint.h>
#include <stdlib.h>
#include <forksafe_mutex.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/times.h>
#include <time.h>
#include <pthread.h>

#include "armpmu_lib.h"
#include "pmu_api.h"


int enable_access_qnx_register();
uint64_t g_evt_over_flow[6] = {0};

static void * get_over_flow_thread(void* arg)
{
    uint i;
    while(1){
        for(i = 0; i < 6; i++){
            g_evt_over_flow[i] += getEvtOverflow(i);
        }
        usleep(20 * 1000);
    }
    pthread_exit(NULL);
    return NULL;
}

int init_pmu(){
        enable_access_qnx_register();
        enableUser();
        pthread_create(NULL, NULL, &get_over_flow_thread, NULL);
}

int start_pmu(){
	enableCounter();
	enablePMU();
	enableEvtCounter();
	setEvtCount();	
	isb();
}

void get_pmu(pmu_paras_t *p)
{
        double ipc_value, cpu_cycles_clock_rate;
        uint64_t clocks, evt_value[6];
        int n;

        clocks = getPMUCount();
	    for(n = 0; n<6; n++){
		    evt_value[n] = getEvtCount(n) + (g_evt_over_flow[n] << 32);
	    }

        if(p != NULL ){
            p->clocks = clocks;
            memcpy(p->evt, evt_value, sizeof(uint64_t)*6);
        }
            

        /* IPC = INST_RETIRED / CPU_CYCLES */
	    /* ipc_value = (double)evt_value[0]/(double)evt_value[1] ; */
        /* cpu_cycles_clock_rate = (double)evt_value[1]/(double)clocks ; */
	    /* printf("\n=======\ncpu clocks:%lld  cpu_cycles_clock_rate:%.3f\n",  clocks, cpu_cycles_clock_rate); */

		/* printf("inst_retired:       %lld\n", evt_value[0]); */
		/* printf("cpu cycles:         %lld\n", evt_value[1]); */
		/* printf("L1D cache:          %lld\n", evt_value[2]); */
		/* printf("L1D cache miss :    %lld\n", evt_value[3]); */
		/* printf("L2D cache:          %lld\n", evt_value[4]); */
		/* printf("L2D cache miss :    %lld\n", evt_value[5]); */

        /* /1* Speculative  accuracy =  INST_RETIRED / INST_SPEC *1/ */
	    /* printf("L1 Miss rate:           %.3f \n", (double)evt_value[3]/(double)(evt_value[2])); */
	    /* printf("L2 Miss rate:           %.3f \n", (double)evt_value[5]/(double)(evt_value[4])); */
	    /* printf("ipc=                    %.3f\n", ipc_value); */
}

int stop_pmu(){
	stopCounter();
	disableEvtCounter();
	disablePMU();
	isb();
    usleep(220*1000);
    for(int i = 0; i < 6; i++){
        g_evt_over_flow[i] = 0;
    }
}
