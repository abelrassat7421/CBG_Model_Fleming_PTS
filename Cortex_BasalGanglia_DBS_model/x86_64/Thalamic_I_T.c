/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__thalamic_i_t
#define _nrn_initial _nrn_initial__thalamic_i_t
#define nrn_cur _nrn_cur__thalamic_i_t
#define _nrn_current _nrn_current__thalamic_i_t
#define nrn_jacob _nrn_jacob__thalamic_i_t
#define nrn_state _nrn_state__thalamic_i_t
#define _net_receive _net_receive__thalamic_i_t 
#define _f_settables _f_settables__thalamic_i_t 
#define settables settables__thalamic_i_t 
#define states states__thalamic_i_t 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define i_t _p[0]
#define i_t_columnindex 0
#define g_T _p[1]
#define g_T_columnindex 1
#define r _p[2]
#define r_columnindex 2
#define ica _p[3]
#define ica_columnindex 3
#define r_inf _p[4]
#define r_inf_columnindex 4
#define tau_r _p[5]
#define tau_r_columnindex 5
#define p_inf _p[6]
#define p_inf_columnindex 6
#define Dr _p[7]
#define Dr_columnindex 7
#define v _p[8]
#define v_columnindex 8
#define _g _p[9]
#define _g_columnindex 9
#define _ion_ica	*_ppvar[0]._pval
#define _ion_dicadv	*_ppvar[1]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_settables(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_thalamic_i_t", _hoc_setdata,
 "settables_thalamic_i_t", _hoc_settables,
 0, 0
};
 
static void _check_settables(double*, Datum*, Datum*, NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, int _type) {
   _check_settables(_p, _ppvar, _thread, _nt);
 }
 /* declare global and static user variables */
#define eca eca_thalamic_i_t
 double eca = 0;
#define usetable usetable_thalamic_i_t
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_thalamic_i_t", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "eca_thalamic_i_t", "mV",
 "i_t_thalamic_i_t", "mA/cm2",
 "g_T_thalamic_i_t", "S/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double r0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "eca_thalamic_i_t", &eca_thalamic_i_t,
 "usetable_thalamic_i_t", &usetable_thalamic_i_t,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"thalamic_i_t",
 "i_t_thalamic_i_t",
 "g_T_thalamic_i_t",
 0,
 0,
 "r_thalamic_i_t",
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 10, _prop);
 	/*initialize range parameters*/
 	i_t = 0;
 	g_T = 0.5;
 	_prop->param = _p;
 	_prop->param_size = 10;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 	_ppvar[0]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[1]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _Thalamic_I_T_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 10, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 thalamic_i_t /home/people/22213094/CBG_Model_Fleming_PTS/Cortex_BasalGanglia_DBS_model/Thalamic_I_T.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_r_inf;
 static double *_t_p_inf;
 static double *_t_tau_r;
static int _reset;
static char *modelname = "Low-threshold calcium current for Thalamic Neuron";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_settables(_threadargsprotocomma_ double);
static int settables(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_settables(_threadargsprotocomma_ double _lv);
 static int _slist1[1], _dlist1[1];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   settables ( _threadargscomma_ v ) ;
   Dr = ( r_inf - r ) / tau_r ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 settables ( _threadargscomma_ v ) ;
 Dr = Dr  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_r )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
   settables ( _threadargscomma_ v ) ;
    r = r + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_r)))*(- ( ( ( r_inf ) ) / tau_r ) / ( ( ( ( - 1.0 ) ) ) / tau_r ) - r) ;
   }
  return 0;
}
 static double _mfac_settables, _tmin_settables;
  static void _check_settables(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_settables =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_settables)/400.; _mfac_settables = 1./_dx;
   for (_i=0, _x=_tmin_settables; _i < 401; _x += _dx, _i++) {
    _f_settables(_p, _ppvar, _thread, _nt, _x);
    _t_r_inf[_i] = r_inf;
    _t_p_inf[_i] = p_inf;
    _t_tau_r[_i] = tau_r;
   }
  }
 }

 static int settables(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _lv) { 
#if 0
_check_settables(_p, _ppvar, _thread, _nt);
#endif
 _n_settables(_p, _ppvar, _thread, _nt, _lv);
 return 0;
 }

 static void _n_settables(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_settables(_p, _ppvar, _thread, _nt, _lv); return; 
}
 _xi = _mfac_settables * (_lv - _tmin_settables);
 if (isnan(_xi)) {
  r_inf = _xi;
  p_inf = _xi;
  tau_r = _xi;
  return;
 }
 if (_xi <= 0.) {
 r_inf = _t_r_inf[0];
 p_inf = _t_p_inf[0];
 tau_r = _t_tau_r[0];
 return; }
 if (_xi >= 400.) {
 r_inf = _t_r_inf[400];
 p_inf = _t_p_inf[400];
 tau_r = _t_tau_r[400];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 r_inf = _t_r_inf[_i] + _theta*(_t_r_inf[_i+1] - _t_r_inf[_i]);
 p_inf = _t_p_inf[_i] + _theta*(_t_p_inf[_i+1] - _t_p_inf[_i]);
 tau_r = _t_tau_r[_i] + _theta*(_t_tau_r[_i+1] - _t_tau_r[_i]);
 }

 
static int  _f_settables ( _threadargsprotocomma_ double _lv ) {
   r_inf = 1.0 / ( 1.0 + exp ( ( _lv + 84.0 ) / 4.0 ) ) ;
   tau_r = ( 28.0 + exp ( - ( _lv + 25.0 ) / 10.5 ) ) ;
   p_inf = 1.0 / ( 1.0 + exp ( - ( _lv + 60.0 ) / 6.2 ) ) ;
    return 0; }
 
static void _hoc_settables(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 
#if 1
 _check_settables(_p, _ppvar, _thread, _nt);
#endif
 _r = 1.;
 settables ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  r = r0;
 {
   settables ( _threadargscomma_ v ) ;
   r = r_inf ;
   }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];

#if 0
 _check_settables(_p, _ppvar, _thread, _nt);
#endif
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ica = g_T * p_inf * p_inf * r * ( v - eca ) ;
   i_t = ica ;
   }
 _current += ica;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = r_columnindex;  _dlist1[0] = Dr_columnindex;
   _t_r_inf = makevector(401*sizeof(double));
   _t_p_inf = makevector(401*sizeof(double));
   _t_tau_r = makevector(401*sizeof(double));
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/people/22213094/CBG_Model_Fleming_PTS/Cortex_BasalGanglia_DBS_model/Thalamic_I_T.mod";
static const char* nmodl_file_text = 
  "TITLE Low-threshold calcium current for Thalamic Neuron\n"
  "\n"
  "COMMENT\n"
  "\n"
  "  Model Reference:\n"
  "\n"
  "  Rubin, J.E. and Terman, D., 2004. \"High frequency stimulation\n"
  "  of the subthalamic nucleus eliminates pathological thalamic\n"
  "  rhythmicity in a computational model.\"\n"
  "  Journal of computational neuroscience, 16(3), pp.211-235.\n"
  "\n"
  "\n"
  "  Implemented by John Fleming - john.fleming@ucdconnect.ie - 06/12/18\n"
  "\n"
  "  Edits:\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "UNITS {\n"
  " (mV) = (millivolt)\n"
  " (mA) = (milliamp)\n"
  " (S) = (siemens)\n"
  "}\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX thalamic_i_t\n"
  "	USEION ca WRITE ica				: Using ca ion, treat the reversal potential as a parameter and write to ica so the total ca current can be tracked\n"
  "	RANGE g_T, i_t					: Calcium current, specific conductance and equilibrium potential\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	eca = 0 (mV)\n"
  "	i_t = 0.0 (mA/cm2)\n"
  "	g_T = 0.5 (S/cm2)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	ica (mA/cm2)\n"
  "	r_inf\n"
  "	tau_r (ms)\n"
  "	p_inf\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	r\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD cnexp\n"
  "	ica = g_T*p_inf*p_inf*r*(v - eca)\n"
  "	i_t = ica 							: Record i_t (just this calcium current) to check it is working\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "\n"
  "INITIAL {\n"
  "	settables(v)\n"
  "	r = r_inf\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "	settables(v)\n"
  "	r' = (r_inf - r)/tau_r\n"
  "}\n"
  "\n"
  "PROCEDURE settables(v) {\n"
  "	TABLE r_inf, p_inf, tau_r FROM -100 TO 100 WITH 400\n"
  "\n"
  "	r_inf = 1/(1+exp((v+84)/4))\n"
  "	tau_r = (28+exp(-(v+25)/10.5))\n"
  "	p_inf = 1/(1+exp(-(v+60)/6.2))\n"
  "}\n"
  "\n"
  "UNITSON\n"
  ;
#endif
