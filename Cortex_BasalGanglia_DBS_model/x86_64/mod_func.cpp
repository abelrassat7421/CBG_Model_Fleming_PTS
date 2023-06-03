#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _Cortical_Axon_I_Kd_reg(void);
extern void _Cortical_Axon_I_Kv_reg(void);
extern void _Cortical_Axon_I_Leak_reg(void);
extern void _Cortical_Axon_I_Na_reg(void);
extern void _Cortical_Soma_I_K_reg(void);
extern void _Cortical_Soma_I_Leak_reg(void);
extern void _Cortical_Soma_I_M_reg(void);
extern void _Cortical_Soma_I_Na_reg(void);
extern void _Destexhe_Static_AMPA_Synapse_reg(void);
extern void _Destexhe_Static_GABAA_Synapse_reg(void);
extern void _Interneuron_I_K_reg(void);
extern void _Interneuron_I_Leak_reg(void);
extern void _Interneuron_I_Na_reg(void);
extern void _myions_reg(void);
extern void _pGPeA_reg(void);
extern void _pSTN_reg(void);
extern void _SynNoise_reg(void);
extern void _Thalamic_I_leak_reg(void);
extern void _Thalamic_I_Na_K_reg(void);
extern void _Thalamic_I_T_reg(void);
extern void _xtra_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"Cortical_Axon_I_Kd.mod\"");
    fprintf(stderr, " \"Cortical_Axon_I_Kv.mod\"");
    fprintf(stderr, " \"Cortical_Axon_I_Leak.mod\"");
    fprintf(stderr, " \"Cortical_Axon_I_Na.mod\"");
    fprintf(stderr, " \"Cortical_Soma_I_K.mod\"");
    fprintf(stderr, " \"Cortical_Soma_I_Leak.mod\"");
    fprintf(stderr, " \"Cortical_Soma_I_M.mod\"");
    fprintf(stderr, " \"Cortical_Soma_I_Na.mod\"");
    fprintf(stderr, " \"Destexhe_Static_AMPA_Synapse.mod\"");
    fprintf(stderr, " \"Destexhe_Static_GABAA_Synapse.mod\"");
    fprintf(stderr, " \"Interneuron_I_K.mod\"");
    fprintf(stderr, " \"Interneuron_I_Leak.mod\"");
    fprintf(stderr, " \"Interneuron_I_Na.mod\"");
    fprintf(stderr, " \"myions.mod\"");
    fprintf(stderr, " \"pGPeA.mod\"");
    fprintf(stderr, " \"pSTN.mod\"");
    fprintf(stderr, " \"SynNoise.mod\"");
    fprintf(stderr, " \"Thalamic_I_leak.mod\"");
    fprintf(stderr, " \"Thalamic_I_Na_K.mod\"");
    fprintf(stderr, " \"Thalamic_I_T.mod\"");
    fprintf(stderr, " \"xtra.mod\"");
    fprintf(stderr, "\n");
  }
  _Cortical_Axon_I_Kd_reg();
  _Cortical_Axon_I_Kv_reg();
  _Cortical_Axon_I_Leak_reg();
  _Cortical_Axon_I_Na_reg();
  _Cortical_Soma_I_K_reg();
  _Cortical_Soma_I_Leak_reg();
  _Cortical_Soma_I_M_reg();
  _Cortical_Soma_I_Na_reg();
  _Destexhe_Static_AMPA_Synapse_reg();
  _Destexhe_Static_GABAA_Synapse_reg();
  _Interneuron_I_K_reg();
  _Interneuron_I_Leak_reg();
  _Interneuron_I_Na_reg();
  _myions_reg();
  _pGPeA_reg();
  _pSTN_reg();
  _SynNoise_reg();
  _Thalamic_I_leak_reg();
  _Thalamic_I_Na_K_reg();
  _Thalamic_I_T_reg();
  _xtra_reg();
}

#if defined(__cplusplus)
}
#endif
