/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#define _pval pval
// clang-format off
#include "md1redef.h"
#include "section_fwd.hpp"
#include "nrniv_mf.h"
#include "md2redef.h"
#include "nrnconf.h"
// clang-format on
#include "neuron/cache/mechanism_range.hpp"
#include <vector>
using std::size_t;
static auto& std_cerr_stream = std::cerr;
static constexpr auto number_of_datum_variables = 0;
static constexpr auto number_of_floating_point_variables = 0;
namespace {
template <typename T>
using _nrn_mechanism_std_vector = std::vector<T>;
using _nrn_model_sorted_token = neuron::model_sorted_token;
using _nrn_mechanism_cache_range = neuron::cache::MechanismRange<number_of_floating_point_variables, number_of_datum_variables>;
using _nrn_mechanism_cache_instance = neuron::cache::MechanismInstance<number_of_floating_point_variables, number_of_datum_variables>;
using _nrn_non_owning_id_without_container = neuron::container::non_owning_identifier_without_container;
template <typename T>
using _nrn_mechanism_field = neuron::mechanism::field<T>;
template <typename... Args>
void _nrn_mechanism_register_data_fields(Args&&... args) {
  neuron::mechanism::register_data_fields(std::forward<Args>(args)...);
}
}
 
#if !NRNGPU
#undef exp
#define exp hoc_Exp
#if NRN_ENABLE_ARCH_INDEP_EXP_POW
#undef pow
#define pow hoc_pow
#endif
#endif
 
#define nrn_init _nrn_init__intfsw
#define _nrn_initial _nrn_initial__intfsw
#define nrn_cur _nrn_cur__intfsw
#define _nrn_current _nrn_current__intfsw
#define nrn_jacob _nrn_jacob__intfsw
#define nrn_state _nrn_state__intfsw
#define _net_receive _net_receive__intfsw 
#define install install__intfsw 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _internalthreadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
#define _internalthreadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *hoc_getarg(int);
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
 static _nrn_mechanism_cache_instance _ml_real{nullptr};
static _nrn_mechanism_cache_range *_ml{&_ml_real};
static size_t _iml{0};
static Datum *_ppvar;
 static int hoc_nrnpointerindex =  -1;
 static Prop* _extcall_prop;
 /* _prop_id kind of shadows _extcall_prop to allow validity checking. */
 static _nrn_non_owning_id_without_container _prop_id{};
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_CountNeighborsR(void);
 static void _hoc_Factorial(void);
 static void _hoc_GetPathEV(void);
 static void _hoc_GetLoopLength(void);
 static void _hoc_GetPathSubPop(void);
 static void _hoc_GetPairDist(void);
 static void _hoc_GetRecurCount(void);
 static void _hoc_GetCCSubPop(void);
 static void _hoc_GetPathR(void);
 static void _hoc_GetWPath(void);
 static void _hoc_GetCC(void);
 static void _hoc_GetCentrality(void);
 static void _hoc_GetCCR(void);
 static void _hoc_install(void);
 static void _hoc_perm(void);
 static void _hoc_predgefunc(void);
 static void _hoc_testmyq(void);
 static void _hoc_testmystack(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mechtype);
#endif
 static void _hoc_setdata();
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 {"setdata_intfsw", _hoc_setdata},
 {"CountNeighborsR_intfsw", _hoc_CountNeighborsR},
 {"Factorial_intfsw", _hoc_Factorial},
 {"GetPathEV_intfsw", _hoc_GetPathEV},
 {"GetLoopLength_intfsw", _hoc_GetLoopLength},
 {"GetPathSubPop_intfsw", _hoc_GetPathSubPop},
 {"GetPairDist_intfsw", _hoc_GetPairDist},
 {"GetRecurCount_intfsw", _hoc_GetRecurCount},
 {"GetCCSubPop_intfsw", _hoc_GetCCSubPop},
 {"GetPathR_intfsw", _hoc_GetPathR},
 {"GetWPath_intfsw", _hoc_GetWPath},
 {"GetCC_intfsw", _hoc_GetCC},
 {"GetCentrality_intfsw", _hoc_GetCentrality},
 {"GetCCR_intfsw", _hoc_GetCCR},
 {"install_intfsw", _hoc_install},
 {"perm_intfsw", _hoc_perm},
 {"predgefunc_intfsw", _hoc_predgefunc},
 {"testmyq_intfsw", _hoc_testmyq},
 {"testmystack_intfsw", _hoc_testmystack},
 {0, 0}
};
 
/* Direct Python call wrappers to density mechanism functions.*/
 static double _npy_CountNeighborsR(Prop*);
 static double _npy_Factorial(Prop*);
 static double _npy_GetPathEV(Prop*);
 static double _npy_GetLoopLength(Prop*);
 static double _npy_GetPathSubPop(Prop*);
 static double _npy_GetPairDist(Prop*);
 static double _npy_GetRecurCount(Prop*);
 static double _npy_GetCCSubPop(Prop*);
 static double _npy_GetPathR(Prop*);
 static double _npy_GetWPath(Prop*);
 static double _npy_GetCC(Prop*);
 static double _npy_GetCentrality(Prop*);
 static double _npy_GetCCR(Prop*);
 static double _npy_install(Prop*);
 static double _npy_perm(Prop*);
 static double _npy_predgefunc(Prop*);
 static double _npy_testmyq(Prop*);
 static double _npy_testmystack(Prop*);
 
static NPyDirectMechFunc npy_direct_func_proc[] = {
 {"CountNeighborsR", _npy_CountNeighborsR},
 {"Factorial", _npy_Factorial},
 {"GetPathEV", _npy_GetPathEV},
 {"GetLoopLength", _npy_GetLoopLength},
 {"GetPathSubPop", _npy_GetPathSubPop},
 {"GetPairDist", _npy_GetPairDist},
 {"GetRecurCount", _npy_GetRecurCount},
 {"GetCCSubPop", _npy_GetCCSubPop},
 {"GetPathR", _npy_GetPathR},
 {"GetWPath", _npy_GetWPath},
 {"GetCC", _npy_GetCC},
 {"GetCentrality", _npy_GetCentrality},
 {"GetCCR", _npy_GetCCR},
 {"install", _npy_install},
 {"perm", _npy_perm},
 {"predgefunc", _npy_predgefunc},
 {"testmyq", _npy_testmyq},
 {"testmystack", _npy_testmystack},
 {0, 0}
};
#define CountNeighborsR CountNeighborsR_intfsw
#define Factorial Factorial_intfsw
#define GetPathEV GetPathEV_intfsw
#define GetLoopLength GetLoopLength_intfsw
#define GetPathSubPop GetPathSubPop_intfsw
#define GetPairDist GetPairDist_intfsw
#define GetRecurCount GetRecurCount_intfsw
#define GetCCSubPop GetCCSubPop_intfsw
#define GetPathR GetPathR_intfsw
#define GetWPath GetWPath_intfsw
#define GetCC GetCC_intfsw
#define GetCentrality GetCentrality_intfsw
#define GetCCR GetCCR_intfsw
#define perm perm_intfsw
#define predgefunc predgefunc_intfsw
#define testmyq testmyq_intfsw
#define testmystack testmystack_intfsw
 extern double CountNeighborsR( );
 extern double Factorial( );
 extern double GetPathEV( );
 extern double GetLoopLength( );
 extern double GetPathSubPop( );
 extern double GetPairDist( );
 extern double GetRecurCount( );
 extern double GetCCSubPop( );
 extern double GetPathR( );
 extern double GetWPath( );
 extern double GetCC( );
 extern double GetCentrality( );
 extern double GetCCR( );
 extern double perm( );
 extern double predgefunc( );
 extern double testmyq( );
 extern double testmystack( );
 /* declare global and static user variables */
 #define gind 0
 #define _gth 0
#define INSTALLED INSTALLED_intfsw
 double INSTALLED = 0;
#define edgefuncid edgefuncid_intfsw
 double edgefuncid = 0;
#define verbose verbose_intfsw
 double verbose = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 {0, 0, 0}
};
 static HocParmUnits _hoc_parm_units[] = {
 {0, 0}
};
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 {"INSTALLED_intfsw", &INSTALLED_intfsw},
 {"verbose_intfsw", &verbose_intfsw},
 {"edgefuncid_intfsw", &edgefuncid_intfsw},
 {0, 0}
};
 static DoubVec hoc_vdoub[] = {
 {0, 0, 0}
};
 static double _sav_indep;
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 _prop_id = _nrn_get_prop_id(_prop);
 neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
_ppvar = _nrn_mechanism_access_dparam(_prop);
 Node * _node = _nrn_mechanism_access_node(_prop);
v = _nrn_mechanism_access_voltage(_node);
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 static void nrn_alloc(Prop*);
static void nrn_init(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
static void nrn_state(_nrn_model_sorted_token const&, NrnThread*, Memb_list*, int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"intfsw",
 0,
 0,
 0,
 0};
 
 /* Used by NrnProperty */
 static _nrn_mechanism_std_vector<double> _parm_default{
 }; 
 
 
extern Prop* need_memb(Symbol*);
static void nrn_alloc(Prop* _prop) {
  Prop *prop_ion{};
  Datum *_ppvar{};
     _nrn_mechanism_cache_instance _ml_real{_prop};
    auto* const _ml = &_ml_real;
    size_t const _iml{};
    assert(_nrn_mechanism_get_num_vars(_prop) == 0);
 	/*initialize range parameters*/
 	 assert(_nrn_mechanism_get_num_vars(_prop) == 0);
 
}
 static void _initlists();
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
void _nrn_thread_table_reg(int, nrn_thread_table_check_t);
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 extern "C" void _intfsw_reg() {
	int _vectorized = 0;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nullptr, nullptr, nullptr, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
 hoc_register_parm_default(_mechtype, &_parm_default);
         hoc_register_npy_direct(_mechtype, npy_direct_func_proc);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  register_nmodl_text_and_filename(_mechtype);
#endif
   _nrn_mechanism_register_data_fields(_mechtype);
  hoc_register_prop_size(_mechtype, 0, 0);
 
    hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 intfsw /Users/bs3667/LDDM/Froemke/NEURON/intfsw.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static const char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int install();
 
/*VERBATIM*/
#include "misc.h"

typedef struct {
  int isz;
  int imaxsz;
  double* p;  
} myvec;

myvec* allocmyvec (int maxsz){
  myvec* pv = (myvec*)malloc(sizeof(myvec));
  if(!pv) return 0x0;
  pv->isz=0;
  pv->imaxsz=maxsz;
  pv->p=(double*)malloc(sizeof(double)*maxsz);
  if(!pv->p) { free(pv); return 0x0; }
  return pv;
}

int freemyvec (myvec** pps) {
  if(!pps || !pps[0]) return 0;
  myvec* ps = pps[0];
  if(ps->p)free(ps->p);
  free(ps);
  pps[0]=0x0;
  return 1;
}

double popmyvec (myvec* pv) {
  if(pv->isz<1) {
    printf("popmyvec ERRA: can't pop empty stack!\n");
    return 0.0;
  }
  double d = pv->p[pv->isz-1]; pv->isz--;
  return d;
}

void popallmyvec (myvec* pv) {
  pv->isz=0;
}

double pushmyvec (myvec* ps,double d) {
  if(ps->isz==ps->imaxsz) {
    printf("pushmyvec realloc\n");
    ps->imaxsz*=2;
    ps->p=(double*)realloc(ps->p,sizeof(double)*ps->imaxsz);
    if(!ps->p){ printf("pushmyvec ERRA: myvec out of memory %d!!\n",ps->imaxsz); return 0.0; }
  }
  ps->p[ps->isz++]=d; 
  return 1.0;  
}

double appendmyvec (myvec* ps,double d) {
  return pushmyvec(ps,d);
}

typedef struct myqnode_ {
  struct myqnode_* pnext;  
  struct myqnode_* pprev;
  int dd;
} myqnode;

myqnode* allocmyqnode() {
  myqnode* p = (myqnode*)malloc(sizeof(myqnode));
  p->pnext=0x0;
  p->pprev=0x0;
  return p;
}

typedef struct {
  myqnode* pfront;
  myqnode* pback;
} myq;

myq* allocmyq() {
  myq* pq = (myq*)malloc(sizeof(myq));
  pq->pfront = pq->pback = 0x0;
  return pq;
}

int freemyq(myq** ppq) {
  myq* pq = *ppq;
  myqnode* ptmp=pq->pback;
  while(pq->pback){
    if(pq->pback->pprev==0x0){
      free(pq->pback);
      pq->pback=0x0;
      pq->pfront=0x0;
      break;
    } else {
      ptmp=pq->pback->pprev;
      free(pq->pback);    
    }
  }
  free(pq);
  ppq[0]=0;
  return 1;
}

int printfrontmyq (myq* pq) {
  if(pq && pq->pfront) {
    printf("front=%d  ",pq->pfront->dd);
    return 1;
  }
  printf("printfrontmyq ERRA: empty front!\n");
  return 0;
}

int printbackmyq (myq* pq) {
  if(pq && pq->pback) {
    printf("back=%d  ",pq->pback->dd);
    return 1;
  }
  printf("printbackmyq ERRA: empty back!\n");
  return 0;
}

int printmyq (myq* pq, int backwards) {
  if(pq){
    int i=0;
    if(backwards){
      myqnode* pnode = pq->pback;
      while(pnode){
        printf("val %d from back = %d\n",i++,pnode->dd);
        pnode = pnode->pprev;
      }
    } else {
      myqnode* pnode = pq->pfront;
      while(pnode){
        printf("val %d from front = %d\n",i++,pnode->dd);
        pnode = pnode->pnext;
      }
    }
    return 1;
  }
  printf("printmyq ERRA: null pointer!\n");
  return 0;
}

int enqmyq (myq* pq,int d) {
  if(pq->pfront==pq->pback) {
    if(!pq->pfront){
      pq->pfront = allocmyqnode();
      pq->pback = pq->pfront;
      pq->pfront->dd=d;
    } else {
      pq->pback = allocmyqnode();
      pq->pback->dd=d;
      pq->pback->pprev = pq->pfront;
      pq->pfront->pnext = pq->pback;
    }
  } else {
    myqnode* pnew = allocmyqnode();
    pnew->dd = d;
    pq->pback->pnext = pnew; 
    pnew->pprev = pq->pback;
    pq->pback = pnew;
  }
  return 1;
}

int emptymyq (myq* pq) {
  if(pq->pfront==0x0) return 1;
  return 0;
}

int deqmyq (myq* pq) {
  if(pq->pfront == pq->pback){
    if(!pq->pfront){
      printf("deqmyq ERRA: can't deq empty q!\n");
      return -1.0;
    } else {
      int d = pq->pfront->dd;
      free(pq->pfront);
      pq->pfront=pq->pback=0x0;
      return d;
    }
  } else {
    myqnode* tmp = pq->pfront;
    int d = tmp->dd;
    pq->pfront = pq->pfront->pnext;
    pq->pfront->pprev = 0x0;
    free(tmp);
    return d;
  }
}

 
double testmystack (  ) {
   double _ltestmystack;
 
/*VERBATIM*/
  myvec* pv = allocmyvec(10);
  printf("created stack with sz %d\n",pv->imaxsz);
  int i;
  for(i=0;i<pv->imaxsz;i++) {
    double d = 41.0 * (i%32) + rand()%100;
    printf("pushing %g onto stack of sz %d\n",d,pv->isz);
    pushmyvec(pv,d);
  }
  printf("test stack realloc by pushing 123.0\n");
  pushmyvec(pv,123.0);
  printf("stack now has %d elements, %d maxsz. contents:\n",pv->isz,pv->imaxsz);
  for(i=0;i<pv->isz;i++)printf("s[%d]=%g\n",i,pv->p[i]);
  printf("popping %d elements. contents:\n",pv->isz);
  while(pv->isz){
    double d = popmyvec(pv);
    printf("popped %g, new sz = %d\n",d,pv->isz);
  }
  printf("can't pop stack now, empty test: ");
  popmyvec(pv);
  freemyvec(&pv);
  printf("freed stack\n");
  return 1.0;
 
return _ltestmystack;
 }
 
static void _hoc_testmystack(void) {
  double _r;
    _r =  testmystack (  );
 hoc_retpushx(_r);
}
 
static double _npy_testmystack(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  testmystack (  );
 return(_r);
}
 
double testmyq (  ) {
   double _ltestmyq;
 
/*VERBATIM*/
  myq* pq = allocmyq();
  printf("created q, empty = %d\n",emptymyq(pq));
  printf("enqueing 10 values:\n");
  int i;
  for(i=0;i<10;i++){
    int d = 41 * (i%32) + rand()%252;
    printf("enqueuing %d...",d);
    enqmyq(pq,d);
    printfrontmyq(pq);
    printbackmyq(pq); printf("\n");
  }
  printf("printing q in forwards order:\n");
  printmyq(pq,0);
  printf("printing q in backwards order:\n");
  printmyq(pq,1);
  printf("testing deq:\n");
  while(!emptymyq(pq)){
    printf("b4 deq: ");
    printfrontmyq(pq); 
    printbackmyq(pq); printf("\n");
    int d = deqmyq(pq);
    printf("dequeued %d\n",d);
    printf("after deq: ");
    printfrontmyq(pq); 
    printbackmyq(pq); printf("\n");
  }
  freemyq(&pq);
  printf("freed myq\n");
  return 1.0;
 
return _ltestmyq;
 }
 
static void _hoc_testmyq(void) {
  double _r;
    _r =  testmyq (  );
 hoc_retpushx(_r);
}
 
static double _npy_testmyq(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  testmyq (  );
 return(_r);
}
 
/*VERBATIM*/
//copy values in valarray who's corresponding entry in binarray != 0 into this vector
//copynz(valvec,binvec)
static double copynz (void* vv) {
  double* pV;
  int n = vector_instance_px(vv,&pV) , iCount = 0 , idx=0;
  int iStartIDx = 0, iEndIDx = n - 1;
  if(ifarg(2)){
    iStartIDx = (int)*getarg(1);
    iEndIDx = (int) *getarg(2);
  }
  if(iEndIDx < iStartIDx || iStartIDx >= n || iEndIDx >= n
                         || iStartIDx<0    || iEndIDx < 0){
    printf("copynz ERRA: invalid indices start=%d end=%d size=%d\n",iStartIDx,iEndIDx,n);
    return -1.0;
  }

  double* pVal,*pBin;

  if(vector_arg_px(1,&pVal)!=n || vector_arg_px(2,&pBin)!=n){
    printf("copynz ERRB: vec args must have size %d!",n);
    return 0.0;
  }

  int iOutSz = 0;
  for(idx=iStartIDx;idx<=iEndIDx;idx++){
    if(pBin[idx]){
      pV[iOutSz++]=pVal[idx];
    }
  }

  vector_resize(pV,iOutSz);

  return (double)iOutSz;
}

//** nnmeandbl()
static double nnmeandbl (double* p,int iStartIDX,int iEndIDX) {
  int iCount=0,idx=0;
  double dSum = 0.0;
  for(idx=iStartIDX;idx<=iEndIDX;idx++){
    if(p[idx]>=0.0){
      dSum+=p[idx];
      iCount++;
    }
  }
  if(iCount>0) return dSum / iCount;
  return -1.0;
} 

//** gzmeandbl()
static double gzmeandbl (double* p,int iStartIDX,int iEndIDX) {
  int iCount=0,idx=0;
  double dSum = 0.0;
  for(idx=iStartIDX;idx<=iEndIDX;idx++){
    if(p[idx]>0.0){
      dSum+=p[idx];
      iCount++;
    }
  }
  if(iCount>0) return dSum / iCount;
  return -1.0;
}

//** gzmean() mean for elements in Vector > 0.0
static double gzmean (void* vv) {
  double* pV;
  int n = vector_instance_px(vv,&pV) , iCount = 0 , idx=0;
  int iStartIDx = 0, iEndIDx = n - 1;
  if(ifarg(2)){
    iStartIDx = (int)*getarg(1);
    iEndIDx = (int) *getarg(2);
  }
  if(iEndIDx < iStartIDx || iStartIDx >= n || iEndIDx >= n
                         || iStartIDx<0    || iEndIDx < 0){
    printf("gzmean ERRA: invalid indices start=%d end=%d size=%d\n",iStartIDx,iEndIDx,n);
    return -1.0;
  }
  return gzmeandbl(pV,iStartIDx,iEndIDx);
}


//** nnmean() mean for elements in Vector >= 0.0
static double nnmean (void* vv) {
  double* pV;
  int n = vector_instance_px(vv,&pV) , iCount = 0 , idx=0;
  int iStartIDx = 0, iEndIDx = n - 1;
  if(ifarg(2)){
    iStartIDx = (int)*getarg(1);
    iEndIDx = (int) *getarg(2);
  }
  if(iEndIDx < iStartIDx || iStartIDx >= n || iEndIDx >= n
                         || iStartIDx<0    || iEndIDx < 0){
    printf("nnmean ERRA: invalid indices start=%d end=%d size=%d\n",iStartIDx,iEndIDx,n);
    return -1.0;
  }
  return nnmeandbl(pV,iStartIDx,iEndIDx);
}
 
double GetCCR (  ) {
   double _lGetCCR;
 
/*VERBATIM*/
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetCC ERRA: problem initializing first arg!\n");
    return 0.0;
  }

  int iCells = pList->isz;
  if(iCells<2){
    printf("GetCC ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  unsigned int* pLen = pList->plen;

  //init vector of distances to each cell , 0 == no path found
  int* pNeighbors = (int*)calloc(iCells,sizeof(int));
  int i = 0, iNeighbors = 0;
  if(!pNeighbors){
    printf("GetCCR ERRE: out of memory!\n");
    FreeListVec(&pList);
    return 0.0;
  }  

  //init vector of avg distances to each cell , 0 == no path found
  double* pCC; 
  int iVecSz = vector_arg_px(2,&pCC);
  if(!pCC || iVecSz < iCells){
    printf("GetCCR ERRE: arg 2 must be a Vector with size %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }  
  memset(pCC,0,sizeof(double)*iVecSz);//init to 0

  //start/end id of cells to find path to
  int iStartID = ifarg(3) ? (int)*getarg(3) : 0,
      iEndID = ifarg(4) ? (int)*getarg(4) : iCells - 1;

  if(iStartID < 0 || iStartID >= iCells ||
     iEndID < 0 || iEndID >= iCells ||
     iStartID >= iEndID){
       printf("GetCCR ERRH: invalid ids start=%d end=%d numcells=%d\n",iStartID,iEndID,iCells);
       FreeListVec(&pList);
       free(pNeighbors);
       return 0.0;
  }

  double dSubsamp = ifarg(5)?*getarg(5):1.0;
  if(dSubsamp<0.0 || dSubsamp>1.0){
    printf("GetCCR ERRH: invalid subsamp = %g , must be btwn 0 and 1\n",dSubsamp);
    FreeListVec(&pList);
    free(pNeighbors);
    return 0.0;
  }

  unsigned int iSeed = ifarg(7)?(unsigned int)*getarg(7):INT_MAX-109754;

  double* pUse = 0; 
  
  if(dSubsamp<1.0){ //if using only a fraction of the cells
     pUse = (double*)malloc(iCells*sizeof(double));
     mcell_ran4(&iSeed, pUse, iCells, 1.0);
  }

  //get id of cell to find paths from
  int myID;

  int* pNeighborID = (int*)calloc(iCells,sizeof(int));

  if( verbose > 0 ) printf("searching from id: ");

  for(myID=0;myID<iCells;myID++) pCC[myID]=-1.0; //set invalid

  for(myID=iStartID;myID<=iEndID;myID++){

    if(verbose > 0 && myID%1000==0)printf("%d ",myID);

    //only use dSubSamp fraction of cells, skip rest
    if(pUse && pUse[myID]>=dSubsamp) continue;

    int idx = 0, youID = 0, youKidID=0 , iNeighbors = 0;

    //mark neighbors of distance == 1
    for(idx=0;idx<pLen[myID];idx++){
      youID = pLV[myID][idx];
      if(youID>=iStartID && youID<=iEndID){
        pNeighbors[youID]=1;      
        pNeighborID[iNeighbors++]=youID;
      }
    }

    if(iNeighbors < 2){
      for(i=0;i<iNeighbors;i++)pNeighbors[pNeighborID[i]]=0;
      continue;
    }

    int iConns = 0 ; 
  
    //this checks # of connections between neighbors of node
    for(i=0;i<iNeighbors;i++){
      if(!pNeighbors[pNeighborID[i]])continue;
      youID=pNeighborID[i];
      for(idx=0;idx<pLen[youID];idx++){
        youKidID=pLV[youID][idx];
        if(youKidID >= iStartID && youKidID <= iEndID && pNeighbors[youKidID]){
          iConns++;
        }
      }
    }
    pCC[myID]=(double)iConns/((double)iNeighbors*(iNeighbors-1));
    for(i=0;i<iNeighbors;i++)pNeighbors[pNeighborID[i]]=0;
  }
 
  free(pNeighborID);
  free(pNeighbors);
  FreeListVec(&pList);
  if(pUse)free(pUse);

  if( verbose > 0 ) printf("\n");

  return  1.0;
 
return _lGetCCR;
 }
 
static void _hoc_GetCCR(void) {
  double _r;
    _r =  GetCCR (  );
 hoc_retpushx(_r);
}
 
static double _npy_GetCCR(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  GetCCR (  );
 return(_r);
}
 
double GetCentrality (  ) {
   double _lGetCentrality;
 
/*VERBATIM*/

  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetCentrality ERRA: problem initializing first arg!\n");
    return 0.0;
  }

  int iCells = pList->isz;
  if(iCells<2){
    printf("GetCentrality ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  unsigned int* pLen = pList->plen;

  //init vector of avg distances to each cell , 0 == no path found
  double* pCE; 
  int iVecSz = vector_arg_px(2,&pCE);
  if(!pCE || iVecSz < iCells){
    printf("GetCCR ERRE: arg 2 must be a Vector with size %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }  
  memset(pCE,0,sizeof(double)*iVecSz);//init to 0

  double dSubsamp = ifarg(3)?*getarg(3):1.0;
  if(dSubsamp<0.0 || dSubsamp>1.0){
    printf("GetCCR ERRH: invalid subsamp = %g , must be btwn 0 and 1\n",dSubsamp);
    FreeListVec(&pList);
    return 0.0;
  }

  unsigned int iSeed = ifarg(4)?(unsigned int)*getarg(4):INT_MAX-109754;

  double* pUse = 0; 
  
  if(dSubsamp<1.0){ //if using only a fraction of the cells
     pUse = (double*)malloc(iCells*sizeof(double));
     mcell_ran4(&iSeed, pUse, iCells, 1.0);
  }

  int s,w,T,v,idx;

  myvec* S = allocmyvec(iCells*2);
  myvec** P = (myvec**)malloc(sizeof(myvec*)*iCells);
  myvec* d = allocmyvec(iCells);
  myvec* sigma = allocmyvec(iCells);
  myvec* di = allocmyvec(iCells);
  for(w=0;w<iCells;w++) P[w]=allocmyvec(iCells);
  for(s=0;s<iCells;s++){
    if(verbose && s%100==0) printf("s=%d\n",s);
    S->isz=0;//empty stack    
    for(w=0;w<iCells;w++) P[w]->isz=0;//empty list
    for(T=0;T<iCells;T++) sigma->p[T]=0; sigma->p[s]=1;
    for(T=0;T<iCells;T++) d->p[T]=-1; d->p[s]=0;
    myq* Q = allocmyq();
    enqmyq(Q,s);
    while(!emptymyq(Q)){
      v = deqmyq(Q);
      pushmyvec(S,v);
      for(idx=0;idx<pLen[v];idx++){
        w = (int) pLV[v][idx];
        if(d->p[w]<0){
          enqmyq(Q,w);
          d->p[w] = d->p[v] + 1;
        }
        if(d->p[w] == d->p[v] + 1){
          sigma->p[w] = sigma->p[w] + sigma->p[v];
          appendmyvec(P[w],v);
        }
      }
    }
    freemyq(&Q);
    for(v=0;v<iCells;v++) di->p[v]=0;
    while(S->isz){
      w = popmyvec(S);
      for(idx=0;idx<P[w]->isz;idx++){
        v=P[w]->p[idx];
        di->p[v] = di->p[v] + (sigma->p[v]/sigma->p[w])*(1.0+di->p[w]);
      }
      if(w!=s) pCE[w] = pCE[w] + di->p[w];
    }
  }

  int N = 0;
  for(s=0;s<iCells;s++) if(pLen[s]) N++;
  if(N>2){
    double scale = 1.0/( (N-1.0)*(N-2.0) );
    for(v=0;v<iCells;v++) if(pLen[v]) pCE[v] *= scale;
  }
  
CEFREE:
  freemyvec(&S);
  for(w=0;w<iCells;w++) freemyvec(&P[w]);
  free(P);
  freemyvec(&d);
  freemyvec(&sigma);
  freemyvec(&di);
  if(pUse)free(pUse);  
  return 1.0;

 
return _lGetCentrality;
 }
 
static void _hoc_GetCentrality(void) {
  double _r;
    _r =  GetCentrality (  );
 hoc_retpushx(_r);
}
 
static double _npy_GetCentrality(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  GetCentrality (  );
 return(_r);
}
 
double GetCC (  ) {
   double _lGetCC;
 
/*VERBATIM*/
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetCC ERRA: problem initializing first arg!\n");
    return -1.0;
  }

  int iCells = pList->isz;
  if(iCells<2){
    printf("GetCC ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return -1.0;
  }

  double** pLV = pList->pv;
  unsigned int* pLen = pList->plen;

  //init vector of distances to each cell , 0 == no path found
  int* pNeighbors = (int*)calloc(iCells,sizeof(int));
  int i = 0, iNeighbors = 0;
  if(!pNeighbors){
    printf("GetCC ERRE: out of memory!\n");
    FreeListVec(&pList);
    return -1.0;
  }  

  //get id of cell to find paths from
  int myID = (int) *getarg(2);
  if(myID < 0 || myID >= iCells){
    printf("GetCC ERRF: invalid id = %d\n",myID);
    FreeListVec(&pList);
    free(pNeighbors);
    return -1.0;
  }

  //start/end id of cells to find path to
  int iStartID = ifarg(3) ? (int)*getarg(3) : 0,
      iEndID = ifarg(4) ? (int)*getarg(4) : iCells - 1;

  if(iStartID < 0 || iStartID >= iCells ||
     iEndID < 0 || iEndID >= iCells ||
     iStartID >= iEndID){
       printf("GetCC ERRH: invalid ids start=%d end=%d numcells=%d\n",iStartID,iEndID,iCells);
       FreeListVec(&pList);
       free(pNeighbors);
       return -1.0;
     }

  int idx = 0, iDist = 1 , youID = 0, youKidID=0;

  int* pNeighborID = (int*)calloc(iCells,sizeof(int));

  //mark neighbors of distance == 1
  for(idx=0;idx<pLen[myID];idx++){
    youID = pLV[myID][idx];
    if(youID>=iStartID && youID<=iEndID){
      pNeighbors[youID]=1;      
      pNeighborID[iNeighbors++]=youID;
    }
  }

  if(iNeighbors < 2){
    FreeListVec(&pList);
    free(pNeighbors);
    return -1.0;
  }

  int iConns = 0; 

  //this checks # of connections between neighbors of node starting from
  for(i=0;i<iNeighbors;i++){
    if(!pNeighbors[pNeighborID[i]])continue;
    youID=pNeighborID[i];
    for(idx=0;idx<pLen[youID];idx++){
      youKidID=pLV[youID][idx];
      if(youKidID >= iStartID && youKidID <= iEndID && pNeighbors[youKidID]){
        iConns++;
      }
    }
  }
 
  free(pNeighborID);
  free(pNeighbors);
  FreeListVec(&pList);

  return  (double)iConns/((double)iNeighbors*(iNeighbors-1));
  
 
return _lGetCC;
 }
 
static void _hoc_GetCC(void) {
  double _r;
    _r =  GetCC (  );
 hoc_retpushx(_r);
}
 
static double _npy_GetCC(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  GetCC (  );
 return(_r);
}
 
double CountNeighborsR (  ) {
   double _lCountNeighborsR;
 
/*VERBATIM*/
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("CountNeighborsR ERRA: problem initializing first arg!\n");
    return 0.0;
  }
 
  int iCells = pList->isz; 
  if(iCells < 2){
    printf("CountNeighborsR ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  unsigned int* pLen = pList->plen;

  //init vector of avg distances to each cell , 0 == no path found
  double* pVD; 
  int iVecSz = vector_arg_px(2,&pVD) , i = 0;
  if(!pVD || iVecSz < iCells){
    printf("CountNeighborsR ERRE: arg 2 must be a Vector with size %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }  
  memset(pVD,0,sizeof(double)*iVecSz);//init to 0

  //get id of cell to find paths from
  int myID = (int) *getarg(3);
  if(myID < 0 || myID >= iCells){
    printf("CountNeighborsR ERRF: invalid id = %d\n",myID);
    FreeListVec(&pList);
    return 0.0;
  }

  //start/end id of cells to search for neighbors of degree iDist 
  int iStartID = (int)*getarg(3),
      iEndID =   (int)*getarg(4),
      iSearchDegree =    (int)*getarg(5);

  double dSubsamp = ifarg(6)?*getarg(6):1.0;

  unsigned int iSeed = ifarg(7)?(unsigned int)*getarg(7):INT_MAX-109754;

  if(iStartID < 0 || iStartID >= iCells ||
     iEndID < 0 || iEndID >= iCells ||
     iStartID >= iEndID){
       printf("CountNeighborsR ERRH: invalid ids start=%d end=%d numcells=%d\n",iStartID,iEndID,iCells);
       FreeListVec(&pList);
       return 0.0;
     }

  //check search degree
  if(iSearchDegree<=0){
    printf("CountNeighborsR ERRI: invalid searchdegree=%d\n",iSearchDegree);
    FreeListVec(&pList);
    return 0.0;
  }

  //init array of cells/neighbors to check
  int* pCheck = (int*)malloc(sizeof(int)*iCells);
  if(!pCheck){
    printf("CountNeighborsR ERRG: out of memory!\n");
    FreeListVec(&pList);
    return 0.0;
  }

  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0, iTmpSz = 0, jdx = 0, iMatches = 0;

  double* pVDTmp = 0, dgzt = 0.0; 
  int* pTmp = 0;
  double* pUse = 0; 
  
  if(dSubsamp<1.0){ //if using only a fraction of the cells
     pUse = (double*)malloc(iCells*sizeof(double));
     mcell_ran4(&iSeed, pUse, iCells, 1.0);
  }

  if( verbose > 0 ) printf("searching from id: ");

  pVDTmp = (double*)calloc(iCells,sizeof(double));
  pTmp = (int*)calloc(iCells,sizeof(int)); 

  for(myID=iStartID;myID<=iEndID;myID++){

    if(verbose > 0 && myID%1000==0)printf("%d ",myID); 

    //only use dSubSamp fraction of cells, skip rest
    if(pUse && pUse[myID]>=dSubsamp) continue;

    iMatches = 0;

    iCheckSz = 0; idx = 0; iDist = 1; youID = 0; youKidID = 0;

    //mark neighbors of distance == 1
    for(idx=0;idx<pLen[myID];idx++){
      youID = pLV[myID][idx];
      if(youID>=iStartID && youID<=iEndID && !pVDTmp[youID]){
        pVDTmp[youID]=(double)iDist;
        pCheck[iCheckSz++]=youID;
      }
    }

    if(iSearchDegree == iDist){
      pVD[myID] = iCheckSz;
      for(idx=0;idx<iCheckSz;idx++) pVDTmp[pCheck[idx]]=0; //reset for next cell
      continue;
    }

    pVDTmp[myID]=1;

    iTmpSz = 0;  jdx=0;

    iDist++;
  
    //this does a breadth-first search but avoids recursion
    while(iCheckSz>0 && iDist<=iSearchDegree){
      iTmpSz = 0;
      for(idx=0;idx<iCheckSz;idx++){
        youID=pCheck[idx];
        for(jdx=0;jdx<pLen[youID];jdx++){
          youKidID=pLV[youID][jdx];
          if(youKidID >= iStartID && youKidID <=iEndID && !pVDTmp[youKidID]){ 
            pTmp[iTmpSz++] = youKidID; //save id of cell to search it's kids on next iteration
            pVDTmp[youKidID]=(double)iDist; //this cell is at iDist away, even if it is also @ a shorter distance
          }
        }
      }
      iCheckSz = iTmpSz;
      
      if(iSearchDegree == iDist){
        pVD[myID] = iCheckSz;
        memset(pVDTmp,0,sizeof(double)*iCells); //reset to 0 for next cell
        break;
      } 

      if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);
      iDist++;
    }
  }

  if(pUse) free(pUse); 
  free(pCheck);
  FreeListVec(&pList);  
  free(pVDTmp); free(pTmp);

  if( verbose > 0 ) printf("\n");

  return 1.0;
 
return _lCountNeighborsR;
 }
 
static void _hoc_CountNeighborsR(void) {
  double _r;
    _r =  CountNeighborsR (  );
 hoc_retpushx(_r);
}
 
static double _npy_CountNeighborsR(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  CountNeighborsR (  );
 return(_r);
}
 
/*VERBATIM*/
double maxval(double* p,int sz)
{
  double dmax = p[0];
  int i = 1;
  for(;i<sz;i++) if(p[i]>dmax) dmax = p[i];
  return dmax;
}

double weightdelaydist(double w,double d)
{
  if(w < 0)
    return -w/d;
  if(w > 0)
    return d/w;
  return DBL_MAX; // no connection means infinite distance
}

double weightdist(double w,double d)
{
  if(w < 0)
    return -w;
  if(w > 0)
    return 1/w;
  return DBL_MAX; // no connection means infinite distance
}

double delaydist(double w,double d)
{
  return d;
}

void printedgefunc(int id)
{
  switch(id){
    case 0:
     printf("weightdelaydist\n");
     break;
    case 1:
     printf("weightdist\n");
     break;
    case 2:
     printf("delaydist\n");
     break;
    default:
     printf("unknown!\n");
     break;
  }
}

 
double predgefunc (  ) {
   double _lpredgefunc;
 
/*VERBATIM*/
  int i;
  if(ifarg(1)){ printf("%d=",(int)*getarg(1)); printedgefunc((int)*getarg(1)); printf("\n"); }    
  else for(i=0;i<3;i++){ printf("%d=",i); printedgefunc(i); printf("\n"); }
  return 0.0;
 
return _lpredgefunc;
 }
 
static void _hoc_predgefunc(void) {
  double _r;
    _r =  predgefunc (  );
 hoc_retpushx(_r);
}
 
static double _npy_predgefunc(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  predgefunc (  );
 return(_r);
}
 
double GetWPath (  ) {
   double _lGetWPath;
 
/*VERBATIM*/

  double* ppre = 0, *ppo = 0, *pwght = 0, *pdel = 0, *pout = 0;
  int iSz,iTmp,i,j,k,l;
  void* voi;

  iSz = vector_arg_px(1,&ppre);

  if(iSz < 1)
  { printf("GetWPath ERRO: invalid size for presynaptic ID Vector (arg 1) %d!\n",iSz);
    return -666.666;
  }

  if( (iTmp=vector_arg_px(2,&ppo)) != iSz)
  { printf("GetWPath ERRA: incorrectly sized postsynaptic ID Vector (arg 2) %d %d!",iSz,iTmp);
    return -666.666;
  }
  if( (iTmp=vector_arg_px(3,&pwght)) != iSz)
  { printf("GetWPath ERRB: incorrectly sized weight Vector (arg 3) %d %d!\n",iSz,iTmp);
    return -666.666;
  }
  if( (iTmp=vector_arg_px(4,&pdel)) != iSz)
  { printf("GetWPath ERRC: incorrectly sized delay Vector (arg 4) %d %d!\n",iSz,iTmp);
    return -666.666;
  }

  int maxid = maxval(ppre,iSz);

  iTmp = maxval(ppo,iSz);
  if(iTmp > maxid) maxid=iTmp;

  voi = vector_arg(5);

  if( (iTmp=vector_arg_px(5,&pout))!= maxid+1 && 0)
  { printf("GetWPath ERRD: incorrectly sized output Vector (arg 5) %d %d!\n",maxid+1,iTmp);
    return -666.666;
  }
  memset(pout,0,sizeof(double)*iTmp);//init to 0

  double (*EdgeFunc)(double,double) = &weightdelaydist;
  int iEdgeFuncID = (int)edgefuncid; 
  if(iEdgeFuncID < 0 || iEdgeFuncID > 2)
  {  printf("GetWPath ERRK: invalid edgedfunc id %d!\n",iEdgeFuncID);
     return -666.666;
  } else if(iEdgeFuncID == 1) EdgeFunc = &weightdist;
    else if(iEdgeFuncID == 2) EdgeFunc = &delaydist;
  if(verbose) printedgefunc(iEdgeFuncID);

 int** adj = (int**) calloc(maxid+1,sizeof(int*));
 if(!adj)
 { printf("GetWPath ERRE: out of memory!\n");
   return -666.666;
 }

 //stores weight of each edge
 //incident from edge is index into pdist
 //incident to edge id is stored in ppo
 double** pdist = (double**) calloc(maxid+1,sizeof(double*));

 int* pcounts = (int*) calloc(maxid+1,sizeof(int));

 //count divergence from each presynaptic cell
 for(i=0;i<iSz;i++)
 { //check for multiple synapses from same source to same target
   if(i+1<iSz && ppre[i]==ppre[i+1] && ppo[i]==ppo[i+1])
   { if(verbose>1) printf("first check double synapse i=%d\n",i);
     while(1)
     { if(i+1>=iSz) break;
       if(ppre[i]!=ppre[i+1] || ppo[i]!=ppo[i+1])
       { //new synapse?
         i--;//move back 1 so get this synapse on next for loop step
         break;
       }
       i++; //move to next synapse
     }      
   }
   pcounts[(int)ppre[i]]++;    //count this one and continue
 }

 //allocate memory for adjacency & distance lists
 for(i=0;i<maxid+1;i++){
   if(pcounts[i]){
     adj[i] = (int*)calloc(pcounts[i],sizeof(int));
     pdist[i] = (double*)calloc(pcounts[i],sizeof(double));
   }
 }

 //index for locations into adjacency lists
 int* pidx = (int*) calloc(maxid+1,sizeof(int));

 //set distance values based on weights and neighbors in adjacency lists based on postsynaptic ids
 for(i=0;i<iSz;i++)
 { int myID = (int)ppre[i];
   if(!pcounts[myID]) continue;//skip cells with 0 divergence
   double dist = EdgeFunc(pwght[i],pdel[i]);
   j=i; //store index of current synapse
   //check for multiple synapses from same source to same target
   if(i+1<iSz && ppre[i]==ppre[i+1] && ppo[i]==ppo[i+1])
   { if(verbose>1) printf("check double syn i=%d\n",i);
     while(1)
     { if(i+1>=iSz) break;
       if(ppre[i]!=ppre[i+1] || ppo[i]!=ppo[i+1])
       { //new synapse?
         i--;//move back 1 so get right synapse on next for loop step
         break;
       }
       if(j!=i) //if didn't count this synapse yet
         dist += EdgeFunc(pwght[i],pdel[i]);
       i++; //move to next synapse to see if it's the same pre,post pair
     }      
   }
   pdist[myID][pidx[myID]] = dist;
   adj[myID][pidx[myID]] = ppo[i];
   pidx[myID]++;
 }

 free(pidx);

 //perform bellman-ford single source shortest path algorithm once for each vertex
 //can improve efficiency by using johnson's algorithm, which uses dijkstra's alg  -- will do later
 double* d = (double*) malloc( (maxid+1)*sizeof(double) ); //distance vector for bellman ford algorithm
 for(i=0;i<=maxid;i++)
 { if(i%100==0) printf("%d ",i);
   if(!pcounts[i])continue;
   for(j=0;j<=maxid;j++) d[j] = DBL_MAX; //initialize distances to +infiniti
   d[i] = 0.0; //distance to self == 0.0
   int changed = 0;
   for(j=0;j<maxid;j++)//apply edge relaxation loop # of vertex-1 times
   { changed=0;
     for(k=0;k<=maxid;k++) //this is just to go thru all edges
     { for(l=0;l<pcounts[k];l++) //go thru all edges of vertex k
       {  if(d[adj[k][l]] > d[k] + pdist[k][l]){//perform edge relaxation
            d[adj[k][l]] = d[k] + pdist[k][l];
            changed=1;
          }
       }
     }
     if(!changed){ if(verbose>1) printf("early term @ j=%d\n",j); break; }
   }

//  int ok = 1;   //make sure no negative cycles
//  for(j=0;j<=maxid && ok;j++)
//  { for(k=0;k<=maxid && ok;k++)
//    { for(l=0;l<pcounts[k];l++)
//      { if( d[adj[k][l]] > d[k] + pdist[k][l] )
//        { ok = 0;
//          break;
//        }
//      }
//    }
//   }
   double avg = 0.0;   //get average distance from vertex i to all other vertices
   int N = 0;
   for(j=0;j<=maxid;j++)
   { if(j!=i && d[j] < DBL_MAX)
     { avg += d[j];
       N++;
     }
   }
   if(N) pout[i] = avg / (double) N;
 }

 free(d);

 //free memory
 free(pcounts);

 for(i=0;i<=maxid;i++){
   if(adj[i]) free(adj[i]);
   if(pdist[i]) free(pdist[i]);
 }

 free(adj);
 free(pdist);

 vector_resize(voi,maxid+1); // pass void* (Vect* ) instead of double*

 return gzmeandbl(pout,0,maxid);

 
return _lGetWPath;
 }
 
static void _hoc_GetWPath(void) {
  double _r;
    _r =  GetWPath (  );
 hoc_retpushx(_r);
}
 
static double _npy_GetWPath(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  GetWPath (  );
 return(_r);
}
 
double GetPathR (  ) {
   double _lGetPathR;
 
/*VERBATIM*/
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetPathEV ERRA: problem initializing first arg!\n");
    return 0.0;
  }
 
  int iCells = pList->isz; 
  if(iCells < 2){
    printf("GetPathEV ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  unsigned int* pLen = pList->plen;

  //init vector of avg distances to each cell , 0 == no path found
  double* pVD; 
  int iVecSz = vector_arg_px(2,&pVD) , i = 0;
  if(!pVD || iVecSz < iCells){
    printf("GetPathEV ERRE: arg 2 must be a Vector with size %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }  
  memset(pVD,0,sizeof(double)*iVecSz);//init to 0

  //start/end id of cells to find path to
  int iStartID = ifarg(3) ? (int)*getarg(3) : 0,
      iEndID = ifarg(4) ? (int)*getarg(4) : iCells - 1,
      iMaxDist = ifarg(5)? (int)*getarg(5): -1;

  double dSubsamp = ifarg(6)?*getarg(6):1.0;

  unsigned int iSeed = ifarg(7)?(unsigned int)*getarg(7):INT_MAX-109754;

  if(iStartID < 0 || iStartID >= iCells ||
     iEndID < 0 || iEndID >= iCells ||
     iStartID >= iEndID){
       printf("GetPathEV ERRH: invalid ids start=%d end=%d numcells=%d\n",iStartID,iEndID,iCells);
       FreeListVec(&pList);
       return 0.0;
     }

  //check max distance
  if(iMaxDist==0){
    printf("GetPathEV ERRI: invalid maxdist=%d\n",iMaxDist);
    FreeListVec(&pList);
    return 0.0;
  }

  //init array of cells/neighbors to check
  int* pCheck;
  pCheck = (int*)malloc(sizeof(int)*iCells);
  if(!pCheck){
    printf("GetPathEV ERRG: out of memory!\n");
    FreeListVec(&pList);
    return 0.0;
  }

  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0, iTmpSz = 0, jdx = 0;

  double* pVDTmp = 0, dgzt = 0.0; 
  int* pTmp = 0;
  double* pUse = 0; 
  
  if(dSubsamp<1.0){ //if using only a fraction of the cells
     pUse = (double*)malloc(iCells*sizeof(double));
     mcell_ran4(&iSeed, pUse, iCells, 1.0);
  }

  pTmp = (int*)calloc(iCells,sizeof(int)); 

  if( verbose > 0 ) printf("searching from id: ");

  pVDTmp = (double*)calloc(iCells,sizeof(double));

  int myID;

  for(myID=iStartID;myID<=iEndID;myID++){

    if(verbose > 0 && myID%1000==0)printf("%d ",myID); 

    //only use dSubSamp fraction of cells, skip rest
    if(pUse && pUse[myID]>=dSubsamp) continue;

    iCheckSz = 0; idx = 0; iDist = 1; youID = 0; youKidID = 0;

    pVDTmp[myID]=1;

    //mark neighbors of distance == 1
    for(idx=0;idx<pLen[myID];idx++){
      youID = pLV[myID][idx];
      if(youID>=iStartID && youID<=iEndID && !pVDTmp[youID]){
        pVDTmp[youID]=(double)iDist;
        pCheck[iCheckSz++]=youID;
      }
    }

    iTmpSz = 0;  jdx=0;

    iDist++;
  
    //this does a breadth-first search but avoids recursion
    while(iCheckSz>0 && (iMaxDist==-1 || iDist<=iMaxDist)){
      iTmpSz = 0;
      for(idx=0;idx<iCheckSz;idx++){
        youID=pCheck[idx];
        for(jdx=0;jdx<pLen[youID];jdx++){
          youKidID=pLV[youID][jdx];
          if(youKidID >= iStartID && youKidID <=iEndID && !pVDTmp[youKidID]){ //found a new connection
            pTmp[iTmpSz++] = youKidID; //save id of cell to search it's kids on next iteration
            pVDTmp[youKidID]=(double)iDist;
          }
        }
      }
      iCheckSz = iTmpSz;
      if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);
      iDist++;
    }

    pVDTmp[myID]=0.0; // distance to self == 0.0
    if((dgzt=gzmeandbl(pVDTmp,iStartID,iEndID))>0.0) pVD[myID]=dgzt;// save mean path length for given cell

    memset(pVDTmp,0,sizeof(double)*iCells);
  }
  
  free(pTmp);
  if(pUse) free(pUse); 
  free(pCheck);
  FreeListVec(&pList);  
  free(pVDTmp);

  if( verbose > 0 ) printf("\n");

  return 1.0;
 
return _lGetPathR;
 }
 
static void _hoc_GetPathR(void) {
  double _r;
    _r =  GetPathR (  );
 hoc_retpushx(_r);
}
 
static double _npy_GetPathR(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  GetPathR (  );
 return(_r);
}
 
double GetCCSubPop (  ) {
   double _lGetCCSubPop;
 
/*VERBATIM*/
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetCCSubPop ERRA: problem initializing first arg!\n");
    return 0.0;
  }
 
  int iCells = pList->isz; 
  if(iCells < 2){
    printf("GetCCSubPop ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  unsigned int* pLen = pList->plen;

  //init vector of distances to each cell , 0 == no path found
  int* pNeighbors = (int*)calloc(iCells,sizeof(int));
  int i = 0, iNeighbors = 0;
  if(!pNeighbors){
    printf("GetCCSubPop ERRE: out of memory!\n");
    FreeListVec(&pList);
    return 0.0;
  }  

  //init vector of avg distances to each cell , 0 == no path found
  double* pCC; 
  int iVecSz = vector_arg_px(2,&pCC);
  if(!pCC || iVecSz < iCells){
    printf("GetCCSubPop ERRE: arg 2 must be a Vector with size %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }  
  memset(pCC,0,sizeof(double)*iVecSz);

  double* pStart,  // bin vec of ids to search from 
          *pEnd;   // bin vec of ids to terminate search on

  if( vector_arg_px(3,&pStart) < iCells || vector_arg_px(4,&pEnd) < iCells){
    printf("GetCCSubPop ERRF: arg 3,4 must be Vectors with size >= %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }
  double dSubsamp = ifarg(5)?*getarg(5):1.0;

  unsigned int iSeed = ifarg(6)?(unsigned int)*getarg(6):INT_MAX-109754;

  double* pUse = 0; 
  
  if(dSubsamp<1.0){ //if using only a fraction of the cells
     pUse = (double*)malloc(iCells*sizeof(double));
     mcell_ran4(&iSeed, pUse, iCells, 1.0);
  }

  //get id of cell to find paths from
  int myID;

  int* pNeighborID = (int*)calloc(iCells,sizeof(int));

  if( verbose > 0 ) printf("searching from id: ");

  for(myID=0;myID<iCells;myID++) pCC[myID]=-1.0; //set invalid

  for(myID=0;myID<iCells;myID++){

    if(!pStart[myID]) continue;

    if(verbose > 0 && myID%1000==0)printf("%d ",myID);

    //only use dSubSamp fraction of cells, skip rest
    if(pUse && pUse[myID]>=dSubsamp) continue;

    int idx = 0, youID = 0, youKidID=0 , iNeighbors = 0;

    //mark neighbors of distance == 1
    for(idx=0;idx<pLen[myID];idx++){
      youID = pLV[myID][idx];
      if(pEnd[youID] && !pNeighbors[youID]){
        pNeighbors[youID]=1;      
        pNeighborID[iNeighbors++]=youID;
      }
    }

    if(iNeighbors < 2){
      for(i=0;i<iNeighbors;i++)pNeighbors[pNeighborID[i]]=0;
      continue;
    }

    int iConns = 0 ; 
  
    //this checks # of connections between neighbors of node
    for(i=0;i<iNeighbors;i++){
      if(!pNeighbors[pNeighborID[i]])continue;
      youID=pNeighborID[i];
      for(idx=0;idx<pLen[youID];idx++){
        youKidID=pLV[youID][idx];
        if(pEnd[youKidID] && pNeighbors[youKidID]){
          iConns++;
        }
      }
    }
    pCC[myID]=(double)iConns/((double)iNeighbors*(iNeighbors-1));
    for(i=0;i<iNeighbors;i++)pNeighbors[pNeighborID[i]]=0;
  }
 
  free(pNeighborID);
  free(pNeighbors);
  FreeListVec(&pList);
  if(pUse)free(pUse);

  if( verbose > 0 ) printf("\n");

  return  1.0;

 
return _lGetCCSubPop;
 }
 
static void _hoc_GetCCSubPop(void) {
  double _r;
    _r =  GetCCSubPop (  );
 hoc_retpushx(_r);
}
 
static double _npy_GetCCSubPop(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  GetCCSubPop (  );
 return(_r);
}
 
double GetRecurCount (  ) {
   double _lGetRecurCount;
 
/*VERBATIM*/
  ListVec* pList;
  int iCells,iFromSz,iThruSz,idx,myID,youID,jdx,iCheckSz,*pVisited,*pCheck;
  double **pLV,*pFrom,*pThru,*pR;

  pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetRecurCount ERRA: problem initializing first arg!\n");
    return 0.0;
  }
 
  iCells = pList->isz; 
  if(iCells < 2){
    printf("GetRecurCount ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  pLV = pList->pv;
  unsigned int* pLen = pList->plen;

  pFrom=pThru=0;
  iFromSz = vector_arg_px(3,&pFrom); iThruSz = vector_arg_px(4,&pThru);
  
  if( iFromSz <= 0 || iThruSz <= 0){
    printf("GetRecurCount ERRF: arg 3,4 bad (fromsz,thrusz)=(%d,%d)\n",iFromSz,iThruSz);
    FreeListVec(&pList);
    return 0.0;
  }

  pVisited = (int*)calloc(iCells,sizeof(int));//which vertices already marked to have children expanded

  pCheck = (int*)malloc(sizeof(int)*iCells);

  pR = vector_newsize(vector_arg(2),iCells);
  memset(pR,0,sizeof(double)*iCells); //zero out output first

  for(myID=0;myID<iCells;myID++) {
    if(!pFrom[myID]) continue;
    iCheckSz = 0; 
    for(idx=0;idx<pLen[myID];idx++){//mark neighbors of distance == 1
      youID = pLV[myID][idx];
      if(!pThru[youID] || pVisited[youID]) continue;
      pCheck[iCheckSz++]=youID;
      pVisited[youID]=1;
    }
    for(idx=0;idx<iCheckSz;idx++) {
      youID = pCheck[idx];
      for(jdx=0;jdx<pLen[youID];jdx++) {
        if(pLV[youID][jdx]==myID) pR[myID]++;
      }
    }
    memset(pVisited,0,sizeof(int)*iCells);
  }
  

  free(pCheck);
  FreeListVec(&pList);  
  free(pVisited);

  if( verbose > 0) printf("\n");

  return 1.0;

 
return _lGetRecurCount;
 }
 
static void _hoc_GetRecurCount(void) {
  double _r;
    _r =  GetRecurCount (  );
 hoc_retpushx(_r);
}
 
static double _npy_GetRecurCount(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  GetRecurCount (  );
 return(_r);
}
 
double GetPairDist (  ) {
   double _lGetPairDist;
 
/*VERBATIM*/
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetPairDist ERRA: problem initializing first arg!\n");
    return 0.0;
  }
 
  int iCells = pList->isz; 
  if(iCells < 2){
    printf("GetPairDist ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  unsigned int* pLen = pList->plen;

  double* pFrom = 0, *pTo = 0;
  int iFromSz = vector_arg_px(3,&pFrom) , iToSz = vector_arg_px(4,&pTo);
  
  if( iFromSz <= 0 || iToSz <= 0){
    printf("GetPairDist ERRF: arg 3,4 bad (fromsz,tosz)=(%d,%d)\n",iFromSz,iToSz);
    FreeListVec(&pList);
    return 0.0;
  }

  int iMinSz = iFromSz * iToSz;

  //init vector of avg distances to each cell , 0 == no path found
  double* pVD; 
  pVD = vector_newsize(vector_arg(2),iMinSz);
  memset(pVD,0,sizeof(double)*iMinSz); //zero out output first

  //init array of cells/neighbors to check
  int* pCheck;
  pCheck = (int*)malloc(sizeof(int)*iCells);
  if(!pCheck){
    printf("GetPairDist ERRG: out of memory!\n");
    FreeListVec(&pList);
    return 0.0;
  }

  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0, iTmpSz = 0, jdx = 0;

  int* pTmp = (int*)calloc(iCells,sizeof(int)); 

  if( verbose > 0 ) printf("searching from id: ");

  int myID , iOff = 0 , kdx = 0;

  int* pVisited = (int*)calloc(iCells,sizeof(int)); //which vertices already marked to have children expanded
  int* pUse = (int*)calloc(iCells,sizeof(int)); //which 'TO' vertices
  int* pMap = (int*)calloc(iCells,sizeof(int)); //index of 'TO' vertices to output index
  for(idx=0;idx<iToSz;idx++){
    pUse[(int)pTo[idx]]=1;
    pMap[(int)pTo[idx]]=idx;
  }

  for(kdx=0;kdx<iFromSz;kdx++,iOff+=iToSz){
    myID=pFrom[kdx];
    if(verbose > 0 && myID%100==0)printf("%d\n",myID);

    iCheckSz = 0; idx = 0; iDist = 1; youID = 0; youKidID = 0;
      
    //mark neighbors of distance == 1
    for(idx=0;idx<pLen[myID];idx++){
      youID = pLV[myID][idx];
      if(pUse[youID]) pVD[ iOff + pMap[youID]  ] = 1; //mark 1st degree neighbor distance as 1
      if(!pVisited[youID]){ 
        pCheck[iCheckSz++]=youID;
        pVisited[youID]=1;
      }
    }

    iTmpSz = 0;  jdx=0;
      
    iDist++;
  
    //this does a breadth-first search but avoids recursion
    while(iCheckSz>0){
      iTmpSz = 0;
      for(idx=0;idx<iCheckSz;idx++){
        youID=pCheck[idx];
        for(jdx=0;jdx<pLen[youID];jdx++){
          youKidID=pLV[youID][jdx];
          if(pUse[youKidID] && !pVD[iOff + pMap[youKidID]])
            pVD[iOff + pMap[youKidID]] = iDist; 
          if(!pVisited[youKidID]){ //found a new connection
            pTmp[iTmpSz++] = youKidID; //save id of cell to search it's kids on next iteration
            pVisited[youKidID]=1;
          }
        }
      }
      iCheckSz = iTmpSz;
      if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);
      iDist++;
    }
    memset(pVisited,0,sizeof(int)*iCells);
  }
  
  free(pTmp);
  free(pCheck);
  FreeListVec(&pList);  
  free(pUse);
  free(pMap);
  free(pVisited);

  if( verbose > 0) printf("\n");

  return 1.0;
 
return _lGetPairDist;
 }
 
static void _hoc_GetPairDist(void) {
  double _r;
    _r =  GetPairDist (  );
 hoc_retpushx(_r);
}
 
static double _npy_GetPairDist(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  GetPairDist (  );
 return(_r);
}
 
double GetPathSubPop (  ) {
   double _lGetPathSubPop;
 
/*VERBATIM*/
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetPathEV ERRA: problem initializing first arg!\n");
    return 0.0;
  }
 
  int iCells = pList->isz; 
  if(iCells < 2){
    printf("GetPathEV ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  unsigned int* pLen = pList->plen;

  //init vector of avg distances to each cell , 0 == no path found
  double* pVD; 
  int iVecSz = vector_arg_px(2,&pVD) , i = 0;
  if(!pVD || iVecSz < iCells){
    printf("GetPathEV ERRE: arg 2 must be a Vector with size %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }  
  memset(pVD,0,sizeof(double)*iVecSz);

  double* pStart,  // bin vec of ids to search from 
          *pEnd;   // bin vec of ids to terminate search on

  if( vector_arg_px(3,&pStart) < iCells || vector_arg_px(4,&pEnd) < iCells){
    printf("GetPathSubPop ERRF: arg 3,4 must be Vectors with size >= %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }
  double dSubsamp = ifarg(5)?*getarg(5):1.0;

  int bSelfLoop = ifarg(6)?(int)*getarg(6):0;

  unsigned int iSeed = ifarg(7)?(unsigned int)*getarg(7):INT_MAX-109754;

  //init array of cells/neighbors to check
  int* pCheck = (int*)malloc(sizeof(int)*iCells);
  if(!pCheck){
    printf("GetPathEV ERRG: out of memory!\n");
    FreeListVec(&pList);
    return 0.0;
  }

  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0, iTmpSz = 0, jdx = 0;

  double  dgzt = 0.0; 
  int* pTmp = 0;
  double* pUse = 0; 
  
  if(dSubsamp<1.0){ //if using only a fraction of the cells
     pUse = (double*)malloc(iCells*sizeof(double));
     mcell_ran4(&iSeed, pUse, iCells, 1.0);
  }

  pTmp = (int*)calloc(iCells,sizeof(int)); 

  if( verbose > 0 ) printf("searching from id: ");

  int* pVDTmp = (int*)calloc(iCells,sizeof(int)) , myID;

  for(myID=0;myID<iCells;myID++){

    if(!pStart[myID]) continue;

    if(verbose > 0 && myID%1000==0)printf("%d ",myID); 

    //only use dSubSamp fraction of cells, skip rest
    if(pUse && pUse[myID]>=dSubsamp) continue;

    unsigned long int iSelfLoopDist = LONG_MAX;
    int bFindThisSelfLoop = bSelfLoop && pEnd[myID]; // search for self loop for this vertex?

    iCheckSz = 0; idx = 0; iDist = 1; youID = 0; youKidID = 0;

    pVDTmp[myID]=1;

    //mark neighbors of distance == 1
    for(idx=0;idx<pLen[myID];idx++){
      youID = pLV[myID][idx];
      if(bFindThisSelfLoop && youID==myID && iDist<iSelfLoopDist) iSelfLoopDist = iDist; //found a self-loop? 
      if(!pVDTmp[youID]){
        pVDTmp[youID]=iDist;
        pCheck[iCheckSz++]=youID;
      }
    }

    iTmpSz = 0;  jdx=0;

    iDist++;
  
    //this does a breadth-first search but avoids recursion
    while(iCheckSz>0){
      iTmpSz = 0;
      for(idx=0;idx<iCheckSz;idx++){
        youID=pCheck[idx];
        for(jdx=0;jdx<pLen[youID];jdx++){
          youKidID=pLV[youID][jdx];
          if(bFindThisSelfLoop && youKidID==myID && iDist<iSelfLoopDist) iSelfLoopDist = iDist; //found a self-loop? 
          if(!pVDTmp[youKidID]){ //found a new connection
            pTmp[iTmpSz++] = youKidID; //save id of cell to search it's kids on next iteration
            pVDTmp[youKidID]=iDist;
          }
        }
      }
      iCheckSz = iTmpSz;
      if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);
      iDist++;
    }

    if(bFindThisSelfLoop && iSelfLoopDist<LONG_MAX){//if checking for this vertex's self-loop dist. and found a self-loop
      pVDTmp[myID] = iSelfLoopDist;
    } else {
      pVDTmp[myID]=0; // distance to self == 0.0
    }
    pVD[myID] = 0.0;
    int N = 0; //take average path length (+ self-loop length if needed) from myID to pEnd cells
    for(idx=0;idx<iCells;idx++){
      if(pEnd[idx] && pVDTmp[idx]){
        pVD[myID] += pVDTmp[idx];
        N++;
      }
    }

    if(N) pVD[myID] /= (double) N; // save mean path (and maybe self-loop) length for given cell

    memset(pVDTmp,0,sizeof(int)*iCells);
  }
  
  free(pTmp);
  if(pUse) free(pUse); 
  free(pCheck);
  FreeListVec(&pList);  
  free(pVDTmp);

  if( verbose > 0 ) printf("\n");

  return 1.0;
 
return _lGetPathSubPop;
 }
 
static void _hoc_GetPathSubPop(void) {
  double _r;
    _r =  GetPathSubPop (  );
 hoc_retpushx(_r);
}
 
static double _npy_GetPathSubPop(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  GetPathSubPop (  );
 return(_r);
}
 
double GetLoopLength (  ) {
   double _lGetLoopLength;
 
/*VERBATIM*/
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetLoopLength ERRA: problem initializing first arg!\n");
    return 0.0;
  }
 
  int iCells = pList->isz; 
  if(iCells < 2){
    printf("GetLoopLength ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  unsigned int* pLen = pList->plen;

  //init vector of avg distances to each cell , 0 == no path found
  double* pVD; 
  int iVecSz = vector_arg_px(2,&pVD) , i = 0;
  if(!pVD || iVecSz < iCells){
    printf("GetLoopLength ERRE: arg 2 must be a Vector with size %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }  
  memset(pVD,0,sizeof(double)*iVecSz);//init to 0

  double* pLoop,  // bin vec of ids to search from 
          *pThru;   // bin vec of ids to terminate search on

  if( vector_arg_px(3,&pLoop) < iCells || vector_arg_px(4,&pThru) < iCells){
    printf("GetLoopLength ERRF: arg 3,4 must be Vectors with size >= %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }
  double dSubsamp = ifarg(5)?*getarg(5):1.0;

  unsigned int iSeed = ifarg(6)?(unsigned int)*getarg(6):INT_MAX-109754;

  //init array of cells/neighbors to check
  int* pCheck = (int*)malloc(sizeof(int)*iCells);
  if(!pCheck){
    printf("GetLoopLength ERRG: out of memory!\n");
    FreeListVec(&pList);
    return 0.0;
  }

  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0, iTmpSz = 0, jdx = 0;

  double  dgzt = 0.0; 
  int* pTmp = 0 , found = 0;
  double* pUse = 0; 
  
  if(dSubsamp<1.0){ //if using only a fraction of the cells
     pUse = (double*)malloc(iCells*sizeof(double));
     mcell_ran4(&iSeed, pUse, iCells, 1.0);
  }

  pTmp = (int*)calloc(iCells,sizeof(int)); 

  if( verbose > 0 ) printf("searching loops from id: ");

  int* pVDTmp = (int*)calloc(iCells,sizeof(int)) , myID;

  for(myID=0;myID<iCells;myID++){

    if(!pLoop[myID]) continue;

    if(verbose > 0 && myID%1000==0)printf("%d ",myID); 

    //only use dSubSamp fraction of cells, skip rest
    if(pUse && pUse[myID]>=dSubsamp) continue;

    iCheckSz = 0; idx = 0; iDist = 1; youID = 0; youKidID = 0; found = 0;

    pVDTmp[myID]=1;

    //mark neighbors of distance == 1
    for(idx=0;idx<pLen[myID];idx++){
      youID = pLV[myID][idx];
      if(youID==myID) {
        found = 1;
        pVD[myID]=iDist;
        iCheckSz=0;
        break;
      }
      if(pThru[youID] && !pVDTmp[youID]){
        pVDTmp[youID]=iDist;
        pCheck[iCheckSz++]=youID;
      }
    }

    iTmpSz = 0;  jdx=0;

    iDist++;
  
    //this does a breadth-first search but avoids recursion
    while(iCheckSz>0){
      iTmpSz = 0;
      for(idx=0;idx<iCheckSz;idx++){
        youID=pCheck[idx];
        for(jdx=0;jdx<pLen[youID];jdx++){
          youKidID=pLV[youID][jdx];
          if(youKidID==myID){
            pVD[myID]=iDist;
            found = 1;
            break;
          }
          if(pThru[youKidID] && !pVDTmp[youKidID]){ //found a new connection
            pTmp[iTmpSz++] = youKidID; //save id of cell to search it's kids on next iteration
            pVDTmp[youKidID]=iDist;
          }
        }
      }
      if(found) break;
      iCheckSz = iTmpSz;
      if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);
      iDist++;
    }
    memset(pVDTmp,0,sizeof(int)*iCells);
  }
  
  free(pTmp);
  if(pUse) free(pUse); 
  free(pCheck);
  FreeListVec(&pList);  
  free(pVDTmp);

  if( verbose > 0 ) printf("\n");

  return 1.0;
 
return _lGetLoopLength;
 }
 
static void _hoc_GetLoopLength(void) {
  double _r;
    _r =  GetLoopLength (  );
 hoc_retpushx(_r);
}
 
static double _npy_GetLoopLength(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  GetLoopLength (  );
 return(_r);
}
 
double GetPathEV (  ) {
   double _lGetPathEV;
 
/*VERBATIM*/
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetPathEV ERRA: problem initializing first arg!\n");
    return 0.0;
  }
 
  int iCells = pList->isz; 
  if(iCells < 2){
    printf("GetPathEV ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  unsigned int* pLen = pList->plen;

  //init vector of distances to each cell , 0 == no path found
  double* pVD; 
  int iVecSz = vector_arg_px(2,&pVD) , i = 0;
  if(!pVD || iVecSz < iCells){
    printf("GetPathEV ERRE: arg 2 must be a Vector with size %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }  
  memset(pVD,0,sizeof(double)*iVecSz);//init to 0

  //get id of cell to find paths from
  int myID = (int) *getarg(3);
  if(myID < 0 || myID >= iCells){
    printf("GetPathEV ERRF: invalid id = %d\n",myID);
    FreeListVec(&pList);
    return 0.0;
  }

  //start/end id of cells to find path to
  int iStartID = ifarg(4) ? (int)*getarg(4) : 0,
      iEndID = ifarg(5) ? (int)*getarg(5) : iCells - 1,
      iMaxDist = ifarg(6)? (int)*getarg(6): -1;

  if(iStartID < 0 || iStartID >= iCells ||
     iEndID < 0 || iEndID >= iCells ||
     iStartID >= iEndID){
       printf("GetPathEV ERRH: invalid ids start=%d end=%d numcells=%d\n",iStartID,iEndID,iCells);
       FreeListVec(&pList);
       return 0.0;
     }

  //check max distance
  if(iMaxDist==0){
    printf("GetPathEV ERRI: invalid maxdist=%d\n",iMaxDist);
    FreeListVec(&pList);
    return 0.0;
  }

  //init array of cells/neighbors to check
  int* pCheck = (int*)malloc(sizeof(int)*iCells);
  if(!pCheck){
    printf("GetPathEV ERRG: out of memory!\n");
    FreeListVec(&pList);
    return 0.0;
  }
  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0;

  pVD[myID]=1;

  //mark neighbors of distance == 1
  for(idx=0;idx<pLen[myID];idx++){
    youID = pLV[myID][idx];
    if(youID>=iStartID && youID<=iEndID && !pVD[youID]){
      pVD[youID]=(double)iDist;
      pCheck[iCheckSz++]=youID;
    }
  }

  int* pTmp = (int*)malloc(sizeof(int)*iCells);
  int iTmpSz = 0 , jdx=0;

  iDist++;
  
  //this does a breadth-first search but avoids deep nesting of recursive version
  while(iCheckSz>0 && (iMaxDist==-1 || iDist<=iMaxDist)){
    iTmpSz = 0;
    for(idx=0;idx<iCheckSz;idx++){
      youID=pCheck[idx];
      for(jdx=0;jdx<pLen[youID];jdx++){
        youKidID=pLV[youID][jdx];
        if(youKidID >= iStartID && youKidID <=iEndID && !pVD[youKidID]){ //found a new connection
          pTmp[iTmpSz++] = youKidID; //save id of cell to search it's kids on next iteration
          pVD[youKidID]=(double)iDist;
        }
      }
    }
    iCheckSz = iTmpSz;
    if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);
    iDist++;
  }

  pVD[myID]=0.0;
 
  free(pCheck);
  free(pTmp);
  FreeListVec(&pList);

  return 1.0;
 
return _lGetPathEV;
 }
 
static void _hoc_GetPathEV(void) {
  double _r;
    _r =  GetPathEV (  );
 hoc_retpushx(_r);
}
 
static double _npy_GetPathEV(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  GetPathEV (  );
 return(_r);
}
 
double Factorial (  ) {
   double _lFactorial;
 
/*VERBATIM*/
  double N = (int)*getarg(1) , i = 0.0;
  double val = 1.0;
  if(N<=1) return 1.0;
  if(N>=171){
    double PI=3.1415926535897932384626433832795;
    double E=2.71828183;
    val=sqrt(2*PI*N)*(pow(N,N)/pow(E,N));
  } else {
    for(i=2.0;i<=N;i++) val*=i;
  }
  return (double) val;  
 
return _lFactorial;
 }
 
static void _hoc_Factorial(void) {
  double _r;
    _r =  Factorial (  );
 hoc_retpushx(_r);
}
 
static double _npy_Factorial(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  Factorial (  );
 return(_r);
}
 
double perm (  ) {
   double _lperm;
 
/*VERBATIM*/
  if(ifarg(3)){
    double N = (int)*getarg(1);
    double R = (int)*getarg(2);
    double b = *getarg(3);
    double val = N/b;
    int i = 0;
    for(i=1;i<R;i++){
      N--;
      val*=(N/b);
    }
    return val;
  } else {
    int N = (int)*getarg(1);
    int R = (int)*getarg(2);
    int val = N;
    int i = 0;
    for(i=1;i<R;i++){
      N--;
      val*=N;
    }
    return (double)val;
  }
 
return _lperm;
 }
 
static void _hoc_perm(void) {
  double _r;
    _r =  perm (  );
 hoc_retpushx(_r);
}
 
static double _npy_perm(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r =  perm (  );
 return(_r);
}
 
static int  install (  ) {
   if ( INSTALLED  == 1.0 ) {
     printf ( "Already installed $Id: intfsw.mod,v 1.50 2009/02/26 18:24:34 samn Exp $ \n" ) ;
     }
   else {
     INSTALLED = 1.0 ;
     
/*VERBATIM*/
 install_vector_method("gzmean" ,gzmean);
 install_vector_method("nnmean" ,nnmean);
 install_vector_method("copynz" ,copynz);
 printf ( "Installed $Id: intfsw.mod,v 1.50 2009/02/26 18:24:34 samn Exp $ \n" ) ;
     }
    return 0; }
 
static void _hoc_install(void) {
  double _r;
    _r = 1.;
 install (  );
 hoc_retpushx(_r);
}
 
static double _npy_install(Prop* _prop) {
    double _r{0.0};
    neuron::legacy::set_globals_from_prop(_prop, _ml_real, _ml, _iml);
  _ppvar = _nrn_mechanism_access_dparam(_prop);
 _r = 1.;
 install (  );
 return(_r);
}

static void initmodel() {
  int _i; double _save;_ninits++;
{

}
}

static void nrn_init(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type){
Node *_nd; double _v; int* _ni; int _cntml;
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto* const _vec_v = _nt->node_voltage_storage();
_ml = &_lmr;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
   _v = _vec_v[_ni[_iml]];
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{
} return _current;
}

static void nrn_state(_nrn_model_sorted_token const& _sorted_token, NrnThread* _nt, Memb_list* _ml_arg, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _cntml;
_nrn_mechanism_cache_range _lmr{_sorted_token, *_nt, *_ml_arg, _type};
auto* const _vec_v = _nt->node_voltage_storage();
_ml = &_lmr;
_ni = _ml_arg->_nodeindices;
_cntml = _ml_arg->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _ppvar = _ml_arg->_pdata[_iml];
 _nd = _ml_arg->_nodelist[_iml];
   _v = _vec_v[_ni[_iml]];
 v=_v;
{
}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if NMODL_TEXT
static void register_nmodl_text_and_filename(int mech_type) {
    const char* nmodl_filename = "/Users/bs3667/LDDM/Froemke/NEURON/intfsw.mod";
    const char* nmodl_file_text = 
  ": $Id: intfsw.mod,v 1.50 2009/02/26 18:24:34 samn Exp $ \n"
  "\n"
  ":* COMMENT\n"
  "COMMENT\n"
  "this file contains functions/utilities for computing the network/graph-theoretic\n"
  "properties of INTF and other networks represented as adjacency lists\n"
  ":** clustering coefficient functions\n"
  "\n"
  "FUNCTION GetCCR(adj,outvec,[startid,endid,subsamp]) gets the clustering coefficient on a range\n"
  "         of cells\n"
  "FUNCTION GetCC -- gets clustering coefficient\n"
  "FUNCTION GetCCSubPop -- get the clustering coefficient between 'sub-populations' of vertices\n"
  "\n"
  ":** path length related functions\n"
  "\n"
  "FUNCTION GetPathR -- gets path length on a range of cells at a time\n"
  "FUNCTION GetWPath -- gets weighted path length , which may be weighted by synaptic weights &\n"
  "         delays\n"
  "FUNCTION GetPairDist -- computes distances between all pairs of vertices, self->self distance==\n"
  "           distance of shortest loop\n"
  "FUNCTION GetPathSubPop -- computes path lengths between sub-populations\n"
  "FUNCTION GetLoopLength -- computes distance to loop back to each node\n"
  "FUNCTION GetPathEV -- gets path length\n"
  "FUNCTION CountNeighborsR -- counts the # of neighbors/outputs of a specified degree on a range\n"
  "         of cells\n"
  "\n"
  ":** miscellaneous functions\n"
  "FUNCTION GetRecurCount -- counts # of recurrent connections\n"
  "FUNCTION Factorial -- computes factorial, if input is too large uses approximation\n"
  "FUNCTION perm - count # of permutations from set of N elements with R selections\n"
  "ENDCOMMENT\n"
  "\n"
  ":* NEURON blocks\n"
  "NEURON {\n"
  "  SUFFIX intfsw\n"
  "  GLOBAL INSTALLED\n"
  "  GLOBAL verbose\n"
  "  GLOBAL edgefuncid : edge-weight-function for GetWPath,0=weightdelaydist,1=weightdist,2=delaydist\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "  INSTALLED=0\n"
  "  verbose=0\n"
  "  edgefuncid=0\n"
  "}\n"
  "\n"
  "VERBATIM\n"
  "#include \"misc.h\"\n"
  "\n"
  "typedef struct {\n"
  "  int isz;\n"
  "  int imaxsz;\n"
  "  double* p;  \n"
  "} myvec;\n"
  "\n"
  "myvec* allocmyvec (int maxsz){\n"
  "  myvec* pv = (myvec*)malloc(sizeof(myvec));\n"
  "  if(!pv) return 0x0;\n"
  "  pv->isz=0;\n"
  "  pv->imaxsz=maxsz;\n"
  "  pv->p=(double*)malloc(sizeof(double)*maxsz);\n"
  "  if(!pv->p) { free(pv); return 0x0; }\n"
  "  return pv;\n"
  "}\n"
  "\n"
  "int freemyvec (myvec** pps) {\n"
  "  if(!pps || !pps[0]) return 0;\n"
  "  myvec* ps = pps[0];\n"
  "  if(ps->p)free(ps->p);\n"
  "  free(ps);\n"
  "  pps[0]=0x0;\n"
  "  return 1;\n"
  "}\n"
  "\n"
  "double popmyvec (myvec* pv) {\n"
  "  if(pv->isz<1) {\n"
  "    printf(\"popmyvec ERRA: can't pop empty stack!\\n\");\n"
  "    return 0.0;\n"
  "  }\n"
  "  double d = pv->p[pv->isz-1]; pv->isz--;\n"
  "  return d;\n"
  "}\n"
  "\n"
  "void popallmyvec (myvec* pv) {\n"
  "  pv->isz=0;\n"
  "}\n"
  "\n"
  "double pushmyvec (myvec* ps,double d) {\n"
  "  if(ps->isz==ps->imaxsz) {\n"
  "    printf(\"pushmyvec realloc\\n\");\n"
  "    ps->imaxsz*=2;\n"
  "    ps->p=(double*)realloc(ps->p,sizeof(double)*ps->imaxsz);\n"
  "    if(!ps->p){ printf(\"pushmyvec ERRA: myvec out of memory %d!!\\n\",ps->imaxsz); return 0.0; }\n"
  "  }\n"
  "  ps->p[ps->isz++]=d; \n"
  "  return 1.0;  \n"
  "}\n"
  "\n"
  "double appendmyvec (myvec* ps,double d) {\n"
  "  return pushmyvec(ps,d);\n"
  "}\n"
  "\n"
  "typedef struct myqnode_ {\n"
  "  struct myqnode_* pnext;  \n"
  "  struct myqnode_* pprev;\n"
  "  int dd;\n"
  "} myqnode;\n"
  "\n"
  "myqnode* allocmyqnode() {\n"
  "  myqnode* p = (myqnode*)malloc(sizeof(myqnode));\n"
  "  p->pnext=0x0;\n"
  "  p->pprev=0x0;\n"
  "  return p;\n"
  "}\n"
  "\n"
  "typedef struct {\n"
  "  myqnode* pfront;\n"
  "  myqnode* pback;\n"
  "} myq;\n"
  "\n"
  "myq* allocmyq() {\n"
  "  myq* pq = (myq*)malloc(sizeof(myq));\n"
  "  pq->pfront = pq->pback = 0x0;\n"
  "  return pq;\n"
  "}\n"
  "\n"
  "int freemyq(myq** ppq) {\n"
  "  myq* pq = *ppq;\n"
  "  myqnode* ptmp=pq->pback;\n"
  "  while(pq->pback){\n"
  "    if(pq->pback->pprev==0x0){\n"
  "      free(pq->pback);\n"
  "      pq->pback=0x0;\n"
  "      pq->pfront=0x0;\n"
  "      break;\n"
  "    } else {\n"
  "      ptmp=pq->pback->pprev;\n"
  "      free(pq->pback);    \n"
  "    }\n"
  "  }\n"
  "  free(pq);\n"
  "  ppq[0]=0;\n"
  "  return 1;\n"
  "}\n"
  "\n"
  "int printfrontmyq (myq* pq) {\n"
  "  if(pq && pq->pfront) {\n"
  "    printf(\"front=%d  \",pq->pfront->dd);\n"
  "    return 1;\n"
  "  }\n"
  "  printf(\"printfrontmyq ERRA: empty front!\\n\");\n"
  "  return 0;\n"
  "}\n"
  "\n"
  "int printbackmyq (myq* pq) {\n"
  "  if(pq && pq->pback) {\n"
  "    printf(\"back=%d  \",pq->pback->dd);\n"
  "    return 1;\n"
  "  }\n"
  "  printf(\"printbackmyq ERRA: empty back!\\n\");\n"
  "  return 0;\n"
  "}\n"
  "\n"
  "int printmyq (myq* pq, int backwards) {\n"
  "  if(pq){\n"
  "    int i=0;\n"
  "    if(backwards){\n"
  "      myqnode* pnode = pq->pback;\n"
  "      while(pnode){\n"
  "        printf(\"val %d from back = %d\\n\",i++,pnode->dd);\n"
  "        pnode = pnode->pprev;\n"
  "      }\n"
  "    } else {\n"
  "      myqnode* pnode = pq->pfront;\n"
  "      while(pnode){\n"
  "        printf(\"val %d from front = %d\\n\",i++,pnode->dd);\n"
  "        pnode = pnode->pnext;\n"
  "      }\n"
  "    }\n"
  "    return 1;\n"
  "  }\n"
  "  printf(\"printmyq ERRA: null pointer!\\n\");\n"
  "  return 0;\n"
  "}\n"
  "\n"
  "int enqmyq (myq* pq,int d) {\n"
  "  if(pq->pfront==pq->pback) {\n"
  "    if(!pq->pfront){\n"
  "      pq->pfront = allocmyqnode();\n"
  "      pq->pback = pq->pfront;\n"
  "      pq->pfront->dd=d;\n"
  "    } else {\n"
  "      pq->pback = allocmyqnode();\n"
  "      pq->pback->dd=d;\n"
  "      pq->pback->pprev = pq->pfront;\n"
  "      pq->pfront->pnext = pq->pback;\n"
  "    }\n"
  "  } else {\n"
  "    myqnode* pnew = allocmyqnode();\n"
  "    pnew->dd = d;\n"
  "    pq->pback->pnext = pnew; \n"
  "    pnew->pprev = pq->pback;\n"
  "    pq->pback = pnew;\n"
  "  }\n"
  "  return 1;\n"
  "}\n"
  "\n"
  "int emptymyq (myq* pq) {\n"
  "  if(pq->pfront==0x0) return 1;\n"
  "  return 0;\n"
  "}\n"
  "\n"
  "int deqmyq (myq* pq) {\n"
  "  if(pq->pfront == pq->pback){\n"
  "    if(!pq->pfront){\n"
  "      printf(\"deqmyq ERRA: can't deq empty q!\\n\");\n"
  "      return -1.0;\n"
  "    } else {\n"
  "      int d = pq->pfront->dd;\n"
  "      free(pq->pfront);\n"
  "      pq->pfront=pq->pback=0x0;\n"
  "      return d;\n"
  "    }\n"
  "  } else {\n"
  "    myqnode* tmp = pq->pfront;\n"
  "    int d = tmp->dd;\n"
  "    pq->pfront = pq->pfront->pnext;\n"
  "    pq->pfront->pprev = 0x0;\n"
  "    free(tmp);\n"
  "    return d;\n"
  "  }\n"
  "}\n"
  "\n"
  "ENDVERBATIM\n"
  "\n"
  "FUNCTION testmystack () {\n"
  "VERBATIM\n"
  "  myvec* pv = allocmyvec(10);\n"
  "  printf(\"created stack with sz %d\\n\",pv->imaxsz);\n"
  "  int i;\n"
  "  for(i=0;i<pv->imaxsz;i++) {\n"
  "    double d = 41.0 * (i%32) + rand()%100;\n"
  "    printf(\"pushing %g onto stack of sz %d\\n\",d,pv->isz);\n"
  "    pushmyvec(pv,d);\n"
  "  }\n"
  "  printf(\"test stack realloc by pushing 123.0\\n\");\n"
  "  pushmyvec(pv,123.0);\n"
  "  printf(\"stack now has %d elements, %d maxsz. contents:\\n\",pv->isz,pv->imaxsz);\n"
  "  for(i=0;i<pv->isz;i++)printf(\"s[%d]=%g\\n\",i,pv->p[i]);\n"
  "  printf(\"popping %d elements. contents:\\n\",pv->isz);\n"
  "  while(pv->isz){\n"
  "    double d = popmyvec(pv);\n"
  "    printf(\"popped %g, new sz = %d\\n\",d,pv->isz);\n"
  "  }\n"
  "  printf(\"can't pop stack now, empty test: \");\n"
  "  popmyvec(pv);\n"
  "  freemyvec(&pv);\n"
  "  printf(\"freed stack\\n\");\n"
  "  return 1.0;\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  "FUNCTION testmyq () {\n"
  "VERBATIM\n"
  "  myq* pq = allocmyq();\n"
  "  printf(\"created q, empty = %d\\n\",emptymyq(pq));\n"
  "  printf(\"enqueing 10 values:\\n\");\n"
  "  int i;\n"
  "  for(i=0;i<10;i++){\n"
  "    int d = 41 * (i%32) + rand()%252;\n"
  "    printf(\"enqueuing %d...\",d);\n"
  "    enqmyq(pq,d);\n"
  "    printfrontmyq(pq);\n"
  "    printbackmyq(pq); printf(\"\\n\");\n"
  "  }\n"
  "  printf(\"printing q in forwards order:\\n\");\n"
  "  printmyq(pq,0);\n"
  "  printf(\"printing q in backwards order:\\n\");\n"
  "  printmyq(pq,1);\n"
  "  printf(\"testing deq:\\n\");\n"
  "  while(!emptymyq(pq)){\n"
  "    printf(\"b4 deq: \");\n"
  "    printfrontmyq(pq); \n"
  "    printbackmyq(pq); printf(\"\\n\");\n"
  "    int d = deqmyq(pq);\n"
  "    printf(\"dequeued %d\\n\",d);\n"
  "    printf(\"after deq: \");\n"
  "    printfrontmyq(pq); \n"
  "    printbackmyq(pq); printf(\"\\n\");\n"
  "  }\n"
  "  freemyq(&pq);\n"
  "  printf(\"freed myq\\n\");\n"
  "  return 1.0;\n"
  "ENDVERBATIM\n"
  "}\n"
  "\n"
  ":* utility functions: copynz(), nnmeandbl(), gzmeandbl(), gzmean(), nnmean() \n"
  "VERBATIM\n"
  "//copy values in valarray who's corresponding entry in binarray != 0 into this vector\n"
  "//copynz(valvec,binvec)\n"
  "static double copynz (void* vv) {\n"
  "  double* pV;\n"
  "  int n = vector_instance_px(vv,&pV) , iCount = 0 , idx=0;\n"
  "  int iStartIDx = 0, iEndIDx = n - 1;\n"
  "  if(ifarg(2)){\n"
  "    iStartIDx = (int)*getarg(1);\n"
  "    iEndIDx = (int) *getarg(2);\n"
  "  }\n"
  "  if(iEndIDx < iStartIDx || iStartIDx >= n || iEndIDx >= n\n"
  "                         || iStartIDx<0    || iEndIDx < 0){\n"
  "    printf(\"copynz ERRA: invalid indices start=%d end=%d size=%d\\n\",iStartIDx,iEndIDx,n);\n"
  "    return -1.0;\n"
  "  }\n"
  "\n"
  "  double* pVal,*pBin;\n"
  "\n"
  "  if(vector_arg_px(1,&pVal)!=n || vector_arg_px(2,&pBin)!=n){\n"
  "    printf(\"copynz ERRB: vec args must have size %d!\",n);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  int iOutSz = 0;\n"
  "  for(idx=iStartIDx;idx<=iEndIDx;idx++){\n"
  "    if(pBin[idx]){\n"
  "      pV[iOutSz++]=pVal[idx];\n"
  "    }\n"
  "  }\n"
  "\n"
  "  vector_resize(pV,iOutSz);\n"
  "\n"
  "  return (double)iOutSz;\n"
  "}\n"
  "\n"
  "//** nnmeandbl()\n"
  "static double nnmeandbl (double* p,int iStartIDX,int iEndIDX) {\n"
  "  int iCount=0,idx=0;\n"
  "  double dSum = 0.0;\n"
  "  for(idx=iStartIDX;idx<=iEndIDX;idx++){\n"
  "    if(p[idx]>=0.0){\n"
  "      dSum+=p[idx];\n"
  "      iCount++;\n"
  "    }\n"
  "  }\n"
  "  if(iCount>0) return dSum / iCount;\n"
  "  return -1.0;\n"
  "} \n"
  "\n"
  "//** gzmeandbl()\n"
  "static double gzmeandbl (double* p,int iStartIDX,int iEndIDX) {\n"
  "  int iCount=0,idx=0;\n"
  "  double dSum = 0.0;\n"
  "  for(idx=iStartIDX;idx<=iEndIDX;idx++){\n"
  "    if(p[idx]>0.0){\n"
  "      dSum+=p[idx];\n"
  "      iCount++;\n"
  "    }\n"
  "  }\n"
  "  if(iCount>0) return dSum / iCount;\n"
  "  return -1.0;\n"
  "}\n"
  "\n"
  "//** gzmean() mean for elements in Vector > 0.0\n"
  "static double gzmean (void* vv) {\n"
  "  double* pV;\n"
  "  int n = vector_instance_px(vv,&pV) , iCount = 0 , idx=0;\n"
  "  int iStartIDx = 0, iEndIDx = n - 1;\n"
  "  if(ifarg(2)){\n"
  "    iStartIDx = (int)*getarg(1);\n"
  "    iEndIDx = (int) *getarg(2);\n"
  "  }\n"
  "  if(iEndIDx < iStartIDx || iStartIDx >= n || iEndIDx >= n\n"
  "                         || iStartIDx<0    || iEndIDx < 0){\n"
  "    printf(\"gzmean ERRA: invalid indices start=%d end=%d size=%d\\n\",iStartIDx,iEndIDx,n);\n"
  "    return -1.0;\n"
  "  }\n"
  "  return gzmeandbl(pV,iStartIDx,iEndIDx);\n"
  "}\n"
  "\n"
  "\n"
  "//** nnmean() mean for elements in Vector >= 0.0\n"
  "static double nnmean (void* vv) {\n"
  "  double* pV;\n"
  "  int n = vector_instance_px(vv,&pV) , iCount = 0 , idx=0;\n"
  "  int iStartIDx = 0, iEndIDx = n - 1;\n"
  "  if(ifarg(2)){\n"
  "    iStartIDx = (int)*getarg(1);\n"
  "    iEndIDx = (int) *getarg(2);\n"
  "  }\n"
  "  if(iEndIDx < iStartIDx || iStartIDx >= n || iEndIDx >= n\n"
  "                         || iStartIDx<0    || iEndIDx < 0){\n"
  "    printf(\"nnmean ERRA: invalid indices start=%d end=%d size=%d\\n\",iStartIDx,iEndIDx,n);\n"
  "    return -1.0;\n"
  "  }\n"
  "  return nnmeandbl(pV,iStartIDx,iEndIDx);\n"
  "}\n"
  "ENDVERBATIM\n"
  "\n"
  ":* GetCCR(adj,outvec,[startid,endid,subsamp]) \n"
  "FUNCTION GetCCR () {\n"
  "  VERBATIM\n"
  "  ListVec* pList = AllocListVec(*hoc_objgetarg(1));\n"
  "  if(!pList){\n"
  "    printf(\"GetCC ERRA: problem initializing first arg!\\n\");\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  int iCells = pList->isz;\n"
  "  if(iCells<2){\n"
  "    printf(\"GetCC ERRB: size of List < 2 !\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  double** pLV = pList->pv;\n"
  "  unsigned int* pLen = pList->plen;\n"
  "\n"
  "  //init vector of distances to each cell , 0 == no path found\n"
  "  int* pNeighbors = (int*)calloc(iCells,sizeof(int));\n"
  "  int i = 0, iNeighbors = 0;\n"
  "  if(!pNeighbors){\n"
  "    printf(\"GetCCR ERRE: out of memory!\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }  \n"
  "\n"
  "  //init vector of avg distances to each cell , 0 == no path found\n"
  "  double* pCC; \n"
  "  int iVecSz = vector_arg_px(2,&pCC);\n"
  "  if(!pCC || iVecSz < iCells){\n"
  "    printf(\"GetCCR ERRE: arg 2 must be a Vector with size %d\\n\",iCells);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }  \n"
  "  memset(pCC,0,sizeof(double)*iVecSz);//init to 0\n"
  "\n"
  "  //start/end id of cells to find path to\n"
  "  int iStartID = ifarg(3) ? (int)*getarg(3) : 0,\n"
  "      iEndID = ifarg(4) ? (int)*getarg(4) : iCells - 1;\n"
  "\n"
  "  if(iStartID < 0 || iStartID >= iCells ||\n"
  "     iEndID < 0 || iEndID >= iCells ||\n"
  "     iStartID >= iEndID){\n"
  "       printf(\"GetCCR ERRH: invalid ids start=%d end=%d numcells=%d\\n\",iStartID,iEndID,iCells);\n"
  "       FreeListVec(&pList);\n"
  "       free(pNeighbors);\n"
  "       return 0.0;\n"
  "  }\n"
  "\n"
  "  double dSubsamp = ifarg(5)?*getarg(5):1.0;\n"
  "  if(dSubsamp<0.0 || dSubsamp>1.0){\n"
  "    printf(\"GetCCR ERRH: invalid subsamp = %g , must be btwn 0 and 1\\n\",dSubsamp);\n"
  "    FreeListVec(&pList);\n"
  "    free(pNeighbors);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  unsigned int iSeed = ifarg(7)?(unsigned int)*getarg(7):INT_MAX-109754;\n"
  "\n"
  "  double* pUse = 0; \n"
  "  \n"
  "  if(dSubsamp<1.0){ //if using only a fraction of the cells\n"
  "     pUse = (double*)malloc(iCells*sizeof(double));\n"
  "     mcell_ran4(&iSeed, pUse, iCells, 1.0);\n"
  "  }\n"
  "\n"
  "  //get id of cell to find paths from\n"
  "  int myID;\n"
  "\n"
  "  int* pNeighborID = (int*)calloc(iCells,sizeof(int));\n"
  "\n"
  "  if( verbose > 0 ) printf(\"searching from id: \");\n"
  "\n"
  "  for(myID=0;myID<iCells;myID++) pCC[myID]=-1.0; //set invalid\n"
  "\n"
  "  for(myID=iStartID;myID<=iEndID;myID++){\n"
  "\n"
  "    if(verbose > 0 && myID%1000==0)printf(\"%d \",myID);\n"
  "\n"
  "    //only use dSubSamp fraction of cells, skip rest\n"
  "    if(pUse && pUse[myID]>=dSubsamp) continue;\n"
  "\n"
  "    int idx = 0, youID = 0, youKidID=0 , iNeighbors = 0;\n"
  "\n"
  "    //mark neighbors of distance == 1\n"
  "    for(idx=0;idx<pLen[myID];idx++){\n"
  "      youID = pLV[myID][idx];\n"
  "      if(youID>=iStartID && youID<=iEndID){\n"
  "        pNeighbors[youID]=1;      \n"
  "        pNeighborID[iNeighbors++]=youID;\n"
  "      }\n"
  "    }\n"
  "\n"
  "    if(iNeighbors < 2){\n"
  "      for(i=0;i<iNeighbors;i++)pNeighbors[pNeighborID[i]]=0;\n"
  "      continue;\n"
  "    }\n"
  "\n"
  "    int iConns = 0 ; \n"
  "  \n"
  "    //this checks # of connections between neighbors of node\n"
  "    for(i=0;i<iNeighbors;i++){\n"
  "      if(!pNeighbors[pNeighborID[i]])continue;\n"
  "      youID=pNeighborID[i];\n"
  "      for(idx=0;idx<pLen[youID];idx++){\n"
  "        youKidID=pLV[youID][idx];\n"
  "        if(youKidID >= iStartID && youKidID <= iEndID && pNeighbors[youKidID]){\n"
  "          iConns++;\n"
  "        }\n"
  "      }\n"
  "    }\n"
  "    pCC[myID]=(double)iConns/((double)iNeighbors*(iNeighbors-1));\n"
  "    for(i=0;i<iNeighbors;i++)pNeighbors[pNeighborID[i]]=0;\n"
  "  }\n"
  " \n"
  "  free(pNeighborID);\n"
  "  free(pNeighbors);\n"
  "  FreeListVec(&pList);\n"
  "  if(pUse)free(pUse);\n"
  "\n"
  "  if( verbose > 0 ) printf(\"\\n\");\n"
  "\n"
  "  return  1.0;\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ":* usage GetCentrality(adjlist,outvec)\n"
  ": based on code from http://www.inf.uni-konstanz.de/algo/publications/b-fabc-01.pdf\n"
  ": and python networkx centrality.py implementation (brandes betweenness centrality)\n"
  "FUNCTION GetCentrality () {\n"
  "  VERBATIM\n"
  "\n"
  "  ListVec* pList = AllocListVec(*hoc_objgetarg(1));\n"
  "  if(!pList){\n"
  "    printf(\"GetCentrality ERRA: problem initializing first arg!\\n\");\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  int iCells = pList->isz;\n"
  "  if(iCells<2){\n"
  "    printf(\"GetCentrality ERRB: size of List < 2 !\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  double** pLV = pList->pv;\n"
  "  unsigned int* pLen = pList->plen;\n"
  "\n"
  "  //init vector of avg distances to each cell , 0 == no path found\n"
  "  double* pCE; \n"
  "  int iVecSz = vector_arg_px(2,&pCE);\n"
  "  if(!pCE || iVecSz < iCells){\n"
  "    printf(\"GetCCR ERRE: arg 2 must be a Vector with size %d\\n\",iCells);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }  \n"
  "  memset(pCE,0,sizeof(double)*iVecSz);//init to 0\n"
  "\n"
  "  double dSubsamp = ifarg(3)?*getarg(3):1.0;\n"
  "  if(dSubsamp<0.0 || dSubsamp>1.0){\n"
  "    printf(\"GetCCR ERRH: invalid subsamp = %g , must be btwn 0 and 1\\n\",dSubsamp);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  unsigned int iSeed = ifarg(4)?(unsigned int)*getarg(4):INT_MAX-109754;\n"
  "\n"
  "  double* pUse = 0; \n"
  "  \n"
  "  if(dSubsamp<1.0){ //if using only a fraction of the cells\n"
  "     pUse = (double*)malloc(iCells*sizeof(double));\n"
  "     mcell_ran4(&iSeed, pUse, iCells, 1.0);\n"
  "  }\n"
  "\n"
  "  int s,w,T,v,idx;\n"
  "\n"
  "  myvec* S = allocmyvec(iCells*2);\n"
  "  myvec** P = (myvec**)malloc(sizeof(myvec*)*iCells);\n"
  "  myvec* d = allocmyvec(iCells);\n"
  "  myvec* sigma = allocmyvec(iCells);\n"
  "  myvec* di = allocmyvec(iCells);\n"
  "  for(w=0;w<iCells;w++) P[w]=allocmyvec(iCells);\n"
  "  for(s=0;s<iCells;s++){\n"
  "    if(verbose && s%100==0) printf(\"s=%d\\n\",s);\n"
  "    S->isz=0;//empty stack    \n"
  "    for(w=0;w<iCells;w++) P[w]->isz=0;//empty list\n"
  "    for(T=0;T<iCells;T++) sigma->p[T]=0; sigma->p[s]=1;\n"
  "    for(T=0;T<iCells;T++) d->p[T]=-1; d->p[s]=0;\n"
  "    myq* Q = allocmyq();\n"
  "    enqmyq(Q,s);\n"
  "    while(!emptymyq(Q)){\n"
  "      v = deqmyq(Q);\n"
  "      pushmyvec(S,v);\n"
  "      for(idx=0;idx<pLen[v];idx++){\n"
  "        w = (int) pLV[v][idx];\n"
  "        if(d->p[w]<0){\n"
  "          enqmyq(Q,w);\n"
  "          d->p[w] = d->p[v] + 1;\n"
  "        }\n"
  "        if(d->p[w] == d->p[v] + 1){\n"
  "          sigma->p[w] = sigma->p[w] + sigma->p[v];\n"
  "          appendmyvec(P[w],v);\n"
  "        }\n"
  "      }\n"
  "    }\n"
  "    freemyq(&Q);\n"
  "    for(v=0;v<iCells;v++) di->p[v]=0;\n"
  "    while(S->isz){\n"
  "      w = popmyvec(S);\n"
  "      for(idx=0;idx<P[w]->isz;idx++){\n"
  "        v=P[w]->p[idx];\n"
  "        di->p[v] = di->p[v] + (sigma->p[v]/sigma->p[w])*(1.0+di->p[w]);\n"
  "      }\n"
  "      if(w!=s) pCE[w] = pCE[w] + di->p[w];\n"
  "    }\n"
  "  }\n"
  "\n"
  "  int N = 0;\n"
  "  for(s=0;s<iCells;s++) if(pLen[s]) N++;\n"
  "  if(N>2){\n"
  "    double scale = 1.0/( (N-1.0)*(N-2.0) );\n"
  "    for(v=0;v<iCells;v++) if(pLen[v]) pCE[v] *= scale;\n"
  "  }\n"
  "  \n"
  "CEFREE:\n"
  "  freemyvec(&S);\n"
  "  for(w=0;w<iCells;w++) freemyvec(&P[w]);\n"
  "  free(P);\n"
  "  freemyvec(&d);\n"
  "  freemyvec(&sigma);\n"
  "  freemyvec(&di);\n"
  "  if(pUse)free(pUse);  \n"
  "  return 1.0;\n"
  "\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ":* usage GetCC(adjlist,myid,[startid,endid])\n"
  ": adjlist == list of vectors specifying connectivity - adjacency list : from row -> to entry in column\n"
  ": myid == id of cell to get clustering coefficient for\n"
  ": startid == min id of cells search can terminate on or go through\n"
  ": endid   == max  '    '   '  '   '  '  '  '  ' '  '  '  '  '  ' \n"
  "FUNCTION GetCC () {\n"
  "  VERBATIM\n"
  "  ListVec* pList = AllocListVec(*hoc_objgetarg(1));\n"
  "  if(!pList){\n"
  "    printf(\"GetCC ERRA: problem initializing first arg!\\n\");\n"
  "    return -1.0;\n"
  "  }\n"
  "\n"
  "  int iCells = pList->isz;\n"
  "  if(iCells<2){\n"
  "    printf(\"GetCC ERRB: size of List < 2 !\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return -1.0;\n"
  "  }\n"
  "\n"
  "  double** pLV = pList->pv;\n"
  "  unsigned int* pLen = pList->plen;\n"
  "\n"
  "  //init vector of distances to each cell , 0 == no path found\n"
  "  int* pNeighbors = (int*)calloc(iCells,sizeof(int));\n"
  "  int i = 0, iNeighbors = 0;\n"
  "  if(!pNeighbors){\n"
  "    printf(\"GetCC ERRE: out of memory!\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return -1.0;\n"
  "  }  \n"
  "\n"
  "  //get id of cell to find paths from\n"
  "  int myID = (int) *getarg(2);\n"
  "  if(myID < 0 || myID >= iCells){\n"
  "    printf(\"GetCC ERRF: invalid id = %d\\n\",myID);\n"
  "    FreeListVec(&pList);\n"
  "    free(pNeighbors);\n"
  "    return -1.0;\n"
  "  }\n"
  "\n"
  "  //start/end id of cells to find path to\n"
  "  int iStartID = ifarg(3) ? (int)*getarg(3) : 0,\n"
  "      iEndID = ifarg(4) ? (int)*getarg(4) : iCells - 1;\n"
  "\n"
  "  if(iStartID < 0 || iStartID >= iCells ||\n"
  "     iEndID < 0 || iEndID >= iCells ||\n"
  "     iStartID >= iEndID){\n"
  "       printf(\"GetCC ERRH: invalid ids start=%d end=%d numcells=%d\\n\",iStartID,iEndID,iCells);\n"
  "       FreeListVec(&pList);\n"
  "       free(pNeighbors);\n"
  "       return -1.0;\n"
  "     }\n"
  "\n"
  "  int idx = 0, iDist = 1 , youID = 0, youKidID=0;\n"
  "\n"
  "  int* pNeighborID = (int*)calloc(iCells,sizeof(int));\n"
  "\n"
  "  //mark neighbors of distance == 1\n"
  "  for(idx=0;idx<pLen[myID];idx++){\n"
  "    youID = pLV[myID][idx];\n"
  "    if(youID>=iStartID && youID<=iEndID){\n"
  "      pNeighbors[youID]=1;      \n"
  "      pNeighborID[iNeighbors++]=youID;\n"
  "    }\n"
  "  }\n"
  "\n"
  "  if(iNeighbors < 2){\n"
  "    FreeListVec(&pList);\n"
  "    free(pNeighbors);\n"
  "    return -1.0;\n"
  "  }\n"
  "\n"
  "  int iConns = 0; \n"
  "\n"
  "  //this checks # of connections between neighbors of node starting from\n"
  "  for(i=0;i<iNeighbors;i++){\n"
  "    if(!pNeighbors[pNeighborID[i]])continue;\n"
  "    youID=pNeighborID[i];\n"
  "    for(idx=0;idx<pLen[youID];idx++){\n"
  "      youKidID=pLV[youID][idx];\n"
  "      if(youKidID >= iStartID && youKidID <= iEndID && pNeighbors[youKidID]){\n"
  "        iConns++;\n"
  "      }\n"
  "    }\n"
  "  }\n"
  " \n"
  "  free(pNeighborID);\n"
  "  free(pNeighbors);\n"
  "  FreeListVec(&pList);\n"
  "\n"
  "  return  (double)iConns/((double)iNeighbors*(iNeighbors-1));\n"
  "  \n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ":* usage CountNeighborsR(adjlist,outvec,startid,endid,degree,subsamp])\n"
  ": adjlist == list of vectors specifying connectivity - adjacency list : from row -> to entry in column\n"
  ": outvec == vector of distances\n"
  ": startid == min id of cells search can terminate on or go through\n"
  ": endid   == max  '    '   '  '   '  '  '  '  ' '  '  '  '  '  ' \n"
  ": degree == distance of neighbors -- counts # of neighbors of EXACT distance specified ONLY\n"
  ": subsamp == specifies fraction btwn 0 and 1 of starting nodes to search\n"
  "FUNCTION CountNeighborsR () {\n"
  "  VERBATIM\n"
  "  ListVec* pList = AllocListVec(*hoc_objgetarg(1));\n"
  "  if(!pList){\n"
  "    printf(\"CountNeighborsR ERRA: problem initializing first arg!\\n\");\n"
  "    return 0.0;\n"
  "  }\n"
  " \n"
  "  int iCells = pList->isz; \n"
  "  if(iCells < 2){\n"
  "    printf(\"CountNeighborsR ERRB: size of List < 2 !\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  double** pLV = pList->pv;\n"
  "  unsigned int* pLen = pList->plen;\n"
  "\n"
  "  //init vector of avg distances to each cell , 0 == no path found\n"
  "  double* pVD; \n"
  "  int iVecSz = vector_arg_px(2,&pVD) , i = 0;\n"
  "  if(!pVD || iVecSz < iCells){\n"
  "    printf(\"CountNeighborsR ERRE: arg 2 must be a Vector with size %d\\n\",iCells);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }  \n"
  "  memset(pVD,0,sizeof(double)*iVecSz);//init to 0\n"
  "\n"
  "  //get id of cell to find paths from\n"
  "  int myID = (int) *getarg(3);\n"
  "  if(myID < 0 || myID >= iCells){\n"
  "    printf(\"CountNeighborsR ERRF: invalid id = %d\\n\",myID);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  //start/end id of cells to search for neighbors of degree iDist \n"
  "  int iStartID = (int)*getarg(3),\n"
  "      iEndID =   (int)*getarg(4),\n"
  "      iSearchDegree =    (int)*getarg(5);\n"
  "\n"
  "  double dSubsamp = ifarg(6)?*getarg(6):1.0;\n"
  "\n"
  "  unsigned int iSeed = ifarg(7)?(unsigned int)*getarg(7):INT_MAX-109754;\n"
  "\n"
  "  if(iStartID < 0 || iStartID >= iCells ||\n"
  "     iEndID < 0 || iEndID >= iCells ||\n"
  "     iStartID >= iEndID){\n"
  "       printf(\"CountNeighborsR ERRH: invalid ids start=%d end=%d numcells=%d\\n\",iStartID,iEndID,iCells);\n"
  "       FreeListVec(&pList);\n"
  "       return 0.0;\n"
  "     }\n"
  "\n"
  "  //check search degree\n"
  "  if(iSearchDegree<=0){\n"
  "    printf(\"CountNeighborsR ERRI: invalid searchdegree=%d\\n\",iSearchDegree);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  //init array of cells/neighbors to check\n"
  "  int* pCheck = (int*)malloc(sizeof(int)*iCells);\n"
  "  if(!pCheck){\n"
  "    printf(\"CountNeighborsR ERRG: out of memory!\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0, iTmpSz = 0, jdx = 0, iMatches = 0;\n"
  "\n"
  "  double* pVDTmp = 0, dgzt = 0.0; \n"
  "  int* pTmp = 0;\n"
  "  double* pUse = 0; \n"
  "  \n"
  "  if(dSubsamp<1.0){ //if using only a fraction of the cells\n"
  "     pUse = (double*)malloc(iCells*sizeof(double));\n"
  "     mcell_ran4(&iSeed, pUse, iCells, 1.0);\n"
  "  }\n"
  "\n"
  "  if( verbose > 0 ) printf(\"searching from id: \");\n"
  "\n"
  "  pVDTmp = (double*)calloc(iCells,sizeof(double));\n"
  "  pTmp = (int*)calloc(iCells,sizeof(int)); \n"
  "\n"
  "  for(myID=iStartID;myID<=iEndID;myID++){\n"
  "\n"
  "    if(verbose > 0 && myID%1000==0)printf(\"%d \",myID); \n"
  "\n"
  "    //only use dSubSamp fraction of cells, skip rest\n"
  "    if(pUse && pUse[myID]>=dSubsamp) continue;\n"
  "\n"
  "    iMatches = 0;\n"
  "\n"
  "    iCheckSz = 0; idx = 0; iDist = 1; youID = 0; youKidID = 0;\n"
  "\n"
  "    //mark neighbors of distance == 1\n"
  "    for(idx=0;idx<pLen[myID];idx++){\n"
  "      youID = pLV[myID][idx];\n"
  "      if(youID>=iStartID && youID<=iEndID && !pVDTmp[youID]){\n"
  "        pVDTmp[youID]=(double)iDist;\n"
  "        pCheck[iCheckSz++]=youID;\n"
  "      }\n"
  "    }\n"
  "\n"
  "    if(iSearchDegree == iDist){\n"
  "      pVD[myID] = iCheckSz;\n"
  "      for(idx=0;idx<iCheckSz;idx++) pVDTmp[pCheck[idx]]=0; //reset for next cell\n"
  "      continue;\n"
  "    }\n"
  "\n"
  "    pVDTmp[myID]=1;\n"
  "\n"
  "    iTmpSz = 0;  jdx=0;\n"
  "\n"
  "    iDist++;\n"
  "  \n"
  "    //this does a breadth-first search but avoids recursion\n"
  "    while(iCheckSz>0 && iDist<=iSearchDegree){\n"
  "      iTmpSz = 0;\n"
  "      for(idx=0;idx<iCheckSz;idx++){\n"
  "        youID=pCheck[idx];\n"
  "        for(jdx=0;jdx<pLen[youID];jdx++){\n"
  "          youKidID=pLV[youID][jdx];\n"
  "          if(youKidID >= iStartID && youKidID <=iEndID && !pVDTmp[youKidID]){ \n"
  "            pTmp[iTmpSz++] = youKidID; //save id of cell to search it's kids on next iteration\n"
  "            pVDTmp[youKidID]=(double)iDist; //this cell is at iDist away, even if it is also @ a shorter distance\n"
  "          }\n"
  "        }\n"
  "      }\n"
  "      iCheckSz = iTmpSz;\n"
  "      \n"
  "      if(iSearchDegree == iDist){\n"
  "        pVD[myID] = iCheckSz;\n"
  "        memset(pVDTmp,0,sizeof(double)*iCells); //reset to 0 for next cell\n"
  "        break;\n"
  "      } \n"
  "\n"
  "      if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);\n"
  "      iDist++;\n"
  "    }\n"
  "  }\n"
  "\n"
  "  if(pUse) free(pUse); \n"
  "  free(pCheck);\n"
  "  FreeListVec(&pList);  \n"
  "  free(pVDTmp); free(pTmp);\n"
  "\n"
  "  if( verbose > 0 ) printf(\"\\n\");\n"
  "\n"
  "  return 1.0;\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ":* utility functions: maxval(), weightdelaydist(), weightdist(), delaydist(), printedgefunc()\n"
  "VERBATIM\n"
  "double maxval(double* p,int sz)\n"
  "{\n"
  "  double dmax = p[0];\n"
  "  int i = 1;\n"
  "  for(;i<sz;i++) if(p[i]>dmax) dmax = p[i];\n"
  "  return dmax;\n"
  "}\n"
  "\n"
  "double weightdelaydist(double w,double d)\n"
  "{\n"
  "  if(w < 0)\n"
  "    return -w/d;\n"
  "  if(w > 0)\n"
  "    return d/w;\n"
  "  return DBL_MAX; // no connection means infinite distance\n"
  "}\n"
  "\n"
  "double weightdist(double w,double d)\n"
  "{\n"
  "  if(w < 0)\n"
  "    return -w;\n"
  "  if(w > 0)\n"
  "    return 1/w;\n"
  "  return DBL_MAX; // no connection means infinite distance\n"
  "}\n"
  "\n"
  "double delaydist(double w,double d)\n"
  "{\n"
  "  return d;\n"
  "}\n"
  "\n"
  "void printedgefunc(int id)\n"
  "{\n"
  "  switch(id){\n"
  "    case 0:\n"
  "     printf(\"weightdelaydist\\n\");\n"
  "     break;\n"
  "    case 1:\n"
  "     printf(\"weightdist\\n\");\n"
  "     break;\n"
  "    case 2:\n"
  "     printf(\"delaydist\\n\");\n"
  "     break;\n"
  "    default:\n"
  "     printf(\"unknown!\\n\");\n"
  "     break;\n"
  "  }\n"
  "}\n"
  "\n"
  "ENDVERBATIM\n"
  "\n"
  ":* FUNCTION predgefunc()\n"
  "FUNCTION predgefunc () {\n"
  "  VERBATIM\n"
  "  int i;\n"
  "  if(ifarg(1)){ printf(\"%d=\",(int)*getarg(1)); printedgefunc((int)*getarg(1)); printf(\"\\n\"); }    \n"
  "  else for(i=0;i<3;i++){ printf(\"%d=\",i); printedgefunc(i); printf(\"\\n\"); }\n"
  "  return 0.0;\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ":* usage GetWPath(preid,poid,weights,delays,outvec,[subsamp])\n"
  ": preid == list of presynaptic IDs\n"
  ": poid == list of postsynaptic IDs\n"
  ": weights == list of weights, excit > 0 , inhib < 0\n"
  ": delays == list of delays \n"
  ": outvec == vector of distances\n"
  ": subsamp == only use specified fraction of synapses , optional\n"
  "FUNCTION GetWPath () {\n"
  "  VERBATIM\n"
  "\n"
  "  double* ppre = 0, *ppo = 0, *pwght = 0, *pdel = 0, *pout = 0;\n"
  "  int iSz,iTmp,i,j,k,l;\n"
  "  void* voi;\n"
  "\n"
  "  iSz = vector_arg_px(1,&ppre);\n"
  "\n"
  "  if(iSz < 1)\n"
  "  { printf(\"GetWPath ERRO: invalid size for presynaptic ID Vector (arg 1) %d!\\n\",iSz);\n"
  "    return -666.666;\n"
  "  }\n"
  "\n"
  "  if( (iTmp=vector_arg_px(2,&ppo)) != iSz)\n"
  "  { printf(\"GetWPath ERRA: incorrectly sized postsynaptic ID Vector (arg 2) %d %d!\",iSz,iTmp);\n"
  "    return -666.666;\n"
  "  }\n"
  "  if( (iTmp=vector_arg_px(3,&pwght)) != iSz)\n"
  "  { printf(\"GetWPath ERRB: incorrectly sized weight Vector (arg 3) %d %d!\\n\",iSz,iTmp);\n"
  "    return -666.666;\n"
  "  }\n"
  "  if( (iTmp=vector_arg_px(4,&pdel)) != iSz)\n"
  "  { printf(\"GetWPath ERRC: incorrectly sized delay Vector (arg 4) %d %d!\\n\",iSz,iTmp);\n"
  "    return -666.666;\n"
  "  }\n"
  "\n"
  "  int maxid = maxval(ppre,iSz);\n"
  "\n"
  "  iTmp = maxval(ppo,iSz);\n"
  "  if(iTmp > maxid) maxid=iTmp;\n"
  "\n"
  "  voi = vector_arg(5);\n"
  "\n"
  "  if( (iTmp=vector_arg_px(5,&pout))!= maxid+1 && 0)\n"
  "  { printf(\"GetWPath ERRD: incorrectly sized output Vector (arg 5) %d %d!\\n\",maxid+1,iTmp);\n"
  "    return -666.666;\n"
  "  }\n"
  "  memset(pout,0,sizeof(double)*iTmp);//init to 0\n"
  "\n"
  "  double (*EdgeFunc)(double,double) = &weightdelaydist;\n"
  "  int iEdgeFuncID = (int)edgefuncid; \n"
  "  if(iEdgeFuncID < 0 || iEdgeFuncID > 2)\n"
  "  {  printf(\"GetWPath ERRK: invalid edgedfunc id %d!\\n\",iEdgeFuncID);\n"
  "     return -666.666;\n"
  "  } else if(iEdgeFuncID == 1) EdgeFunc = &weightdist;\n"
  "    else if(iEdgeFuncID == 2) EdgeFunc = &delaydist;\n"
  "  if(verbose) printedgefunc(iEdgeFuncID);\n"
  "\n"
  " int** adj = (int**) calloc(maxid+1,sizeof(int*));\n"
  " if(!adj)\n"
  " { printf(\"GetWPath ERRE: out of memory!\\n\");\n"
  "   return -666.666;\n"
  " }\n"
  "\n"
  " //stores weight of each edge\n"
  " //incident from edge is index into pdist\n"
  " //incident to edge id is stored in ppo\n"
  " double** pdist = (double**) calloc(maxid+1,sizeof(double*));\n"
  "\n"
  " int* pcounts = (int*) calloc(maxid+1,sizeof(int));\n"
  "\n"
  " //count divergence from each presynaptic cell\n"
  " for(i=0;i<iSz;i++)\n"
  " { //check for multiple synapses from same source to same target\n"
  "   if(i+1<iSz && ppre[i]==ppre[i+1] && ppo[i]==ppo[i+1])\n"
  "   { if(verbose>1) printf(\"first check double synapse i=%d\\n\",i);\n"
  "     while(1)\n"
  "     { if(i+1>=iSz) break;\n"
  "       if(ppre[i]!=ppre[i+1] || ppo[i]!=ppo[i+1])\n"
  "       { //new synapse?\n"
  "         i--;//move back 1 so get this synapse on next for loop step\n"
  "         break;\n"
  "       }\n"
  "       i++; //move to next synapse\n"
  "     }      \n"
  "   }\n"
  "   pcounts[(int)ppre[i]]++;    //count this one and continue\n"
  " }\n"
  "\n"
  " //allocate memory for adjacency & distance lists\n"
  " for(i=0;i<maxid+1;i++){\n"
  "   if(pcounts[i]){\n"
  "     adj[i] = (int*)calloc(pcounts[i],sizeof(int));\n"
  "     pdist[i] = (double*)calloc(pcounts[i],sizeof(double));\n"
  "   }\n"
  " }\n"
  "\n"
  " //index for locations into adjacency lists\n"
  " int* pidx = (int*) calloc(maxid+1,sizeof(int));\n"
  "\n"
  " //set distance values based on weights and neighbors in adjacency lists based on postsynaptic ids\n"
  " for(i=0;i<iSz;i++)\n"
  " { int myID = (int)ppre[i];\n"
  "   if(!pcounts[myID]) continue;//skip cells with 0 divergence\n"
  "   double dist = EdgeFunc(pwght[i],pdel[i]);\n"
  "   j=i; //store index of current synapse\n"
  "   //check for multiple synapses from same source to same target\n"
  "   if(i+1<iSz && ppre[i]==ppre[i+1] && ppo[i]==ppo[i+1])\n"
  "   { if(verbose>1) printf(\"check double syn i=%d\\n\",i);\n"
  "     while(1)\n"
  "     { if(i+1>=iSz) break;\n"
  "       if(ppre[i]!=ppre[i+1] || ppo[i]!=ppo[i+1])\n"
  "       { //new synapse?\n"
  "         i--;//move back 1 so get right synapse on next for loop step\n"
  "         break;\n"
  "       }\n"
  "       if(j!=i) //if didn't count this synapse yet\n"
  "         dist += EdgeFunc(pwght[i],pdel[i]);\n"
  "       i++; //move to next synapse to see if it's the same pre,post pair\n"
  "     }      \n"
  "   }\n"
  "   pdist[myID][pidx[myID]] = dist;\n"
  "   adj[myID][pidx[myID]] = ppo[i];\n"
  "   pidx[myID]++;\n"
  " }\n"
  "\n"
  " free(pidx);\n"
  "\n"
  " //perform bellman-ford single source shortest path algorithm once for each vertex\n"
  " //can improve efficiency by using johnson's algorithm, which uses dijkstra's alg  -- will do later\n"
  " double* d = (double*) malloc( (maxid+1)*sizeof(double) ); //distance vector for bellman ford algorithm\n"
  " for(i=0;i<=maxid;i++)\n"
  " { if(i%100==0) printf(\"%d \",i);\n"
  "   if(!pcounts[i])continue;\n"
  "   for(j=0;j<=maxid;j++) d[j] = DBL_MAX; //initialize distances to +infiniti\n"
  "   d[i] = 0.0; //distance to self == 0.0\n"
  "   int changed = 0;\n"
  "   for(j=0;j<maxid;j++)//apply edge relaxation loop # of vertex-1 times\n"
  "   { changed=0;\n"
  "     for(k=0;k<=maxid;k++) //this is just to go thru all edges\n"
  "     { for(l=0;l<pcounts[k];l++) //go thru all edges of vertex k\n"
  "       {  if(d[adj[k][l]] > d[k] + pdist[k][l]){//perform edge relaxation\n"
  "            d[adj[k][l]] = d[k] + pdist[k][l];\n"
  "            changed=1;\n"
  "          }\n"
  "       }\n"
  "     }\n"
  "     if(!changed){ if(verbose>1) printf(\"early term @ j=%d\\n\",j); break; }\n"
  "   }\n"
  "\n"
  "//  int ok = 1;   //make sure no negative cycles\n"
  "//  for(j=0;j<=maxid && ok;j++)\n"
  "//  { for(k=0;k<=maxid && ok;k++)\n"
  "//    { for(l=0;l<pcounts[k];l++)\n"
  "//      { if( d[adj[k][l]] > d[k] + pdist[k][l] )\n"
  "//        { ok = 0;\n"
  "//          break;\n"
  "//        }\n"
  "//      }\n"
  "//    }\n"
  "//   }\n"
  "   double avg = 0.0;   //get average distance from vertex i to all other vertices\n"
  "   int N = 0;\n"
  "   for(j=0;j<=maxid;j++)\n"
  "   { if(j!=i && d[j] < DBL_MAX)\n"
  "     { avg += d[j];\n"
  "       N++;\n"
  "     }\n"
  "   }\n"
  "   if(N) pout[i] = avg / (double) N;\n"
  " }\n"
  "\n"
  " free(d);\n"
  "\n"
  " //free memory\n"
  " free(pcounts);\n"
  "\n"
  " for(i=0;i<=maxid;i++){\n"
  "   if(adj[i]) free(adj[i]);\n"
  "   if(pdist[i]) free(pdist[i]);\n"
  " }\n"
  "\n"
  " free(adj);\n"
  " free(pdist);\n"
  "\n"
  " vector_resize(voi,maxid+1); // pass void* (Vect* ) instead of double*\n"
  "\n"
  " return gzmeandbl(pout,0,maxid);\n"
  "\n"
  " ENDVERBATIM\n"
  "}\n"
  "\n"
  ":* usage GetPathR(adjlist,outvec,[startid,endid,maxdist,subsamp])\n"
  ": adjlist == list of vectors specifying connectivity - adjacency list : from row -> to entry in column\n"
  ": outvec == vector of distances\n"
  ": startid == min id of cells search can terminate on or go through\n"
  ": endid   == max  '    '   '  '   '  '  '  '  ' '  '  '  '  '  ' \n"
  ": maxdist == max # of connections to allow hops over\n"
  ": subsamp == perform calculation on % of cells, default == 1\n"
  "FUNCTION GetPathR () {\n"
  "  VERBATIM\n"
  "  ListVec* pList = AllocListVec(*hoc_objgetarg(1));\n"
  "  if(!pList){\n"
  "    printf(\"GetPathEV ERRA: problem initializing first arg!\\n\");\n"
  "    return 0.0;\n"
  "  }\n"
  " \n"
  "  int iCells = pList->isz; \n"
  "  if(iCells < 2){\n"
  "    printf(\"GetPathEV ERRB: size of List < 2 !\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  double** pLV = pList->pv;\n"
  "  unsigned int* pLen = pList->plen;\n"
  "\n"
  "  //init vector of avg distances to each cell , 0 == no path found\n"
  "  double* pVD; \n"
  "  int iVecSz = vector_arg_px(2,&pVD) , i = 0;\n"
  "  if(!pVD || iVecSz < iCells){\n"
  "    printf(\"GetPathEV ERRE: arg 2 must be a Vector with size %d\\n\",iCells);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }  \n"
  "  memset(pVD,0,sizeof(double)*iVecSz);//init to 0\n"
  "\n"
  "  //start/end id of cells to find path to\n"
  "  int iStartID = ifarg(3) ? (int)*getarg(3) : 0,\n"
  "      iEndID = ifarg(4) ? (int)*getarg(4) : iCells - 1,\n"
  "      iMaxDist = ifarg(5)? (int)*getarg(5): -1;\n"
  "\n"
  "  double dSubsamp = ifarg(6)?*getarg(6):1.0;\n"
  "\n"
  "  unsigned int iSeed = ifarg(7)?(unsigned int)*getarg(7):INT_MAX-109754;\n"
  "\n"
  "  if(iStartID < 0 || iStartID >= iCells ||\n"
  "     iEndID < 0 || iEndID >= iCells ||\n"
  "     iStartID >= iEndID){\n"
  "       printf(\"GetPathEV ERRH: invalid ids start=%d end=%d numcells=%d\\n\",iStartID,iEndID,iCells);\n"
  "       FreeListVec(&pList);\n"
  "       return 0.0;\n"
  "     }\n"
  "\n"
  "  //check max distance\n"
  "  if(iMaxDist==0){\n"
  "    printf(\"GetPathEV ERRI: invalid maxdist=%d\\n\",iMaxDist);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  //init array of cells/neighbors to check\n"
  "  int* pCheck;\n"
  "  pCheck = (int*)malloc(sizeof(int)*iCells);\n"
  "  if(!pCheck){\n"
  "    printf(\"GetPathEV ERRG: out of memory!\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0, iTmpSz = 0, jdx = 0;\n"
  "\n"
  "  double* pVDTmp = 0, dgzt = 0.0; \n"
  "  int* pTmp = 0;\n"
  "  double* pUse = 0; \n"
  "  \n"
  "  if(dSubsamp<1.0){ //if using only a fraction of the cells\n"
  "     pUse = (double*)malloc(iCells*sizeof(double));\n"
  "     mcell_ran4(&iSeed, pUse, iCells, 1.0);\n"
  "  }\n"
  "\n"
  "  pTmp = (int*)calloc(iCells,sizeof(int)); \n"
  "\n"
  "  if( verbose > 0 ) printf(\"searching from id: \");\n"
  "\n"
  "  pVDTmp = (double*)calloc(iCells,sizeof(double));\n"
  "\n"
  "  int myID;\n"
  "\n"
  "  for(myID=iStartID;myID<=iEndID;myID++){\n"
  "\n"
  "    if(verbose > 0 && myID%1000==0)printf(\"%d \",myID); \n"
  "\n"
  "    //only use dSubSamp fraction of cells, skip rest\n"
  "    if(pUse && pUse[myID]>=dSubsamp) continue;\n"
  "\n"
  "    iCheckSz = 0; idx = 0; iDist = 1; youID = 0; youKidID = 0;\n"
  "\n"
  "    pVDTmp[myID]=1;\n"
  "\n"
  "    //mark neighbors of distance == 1\n"
  "    for(idx=0;idx<pLen[myID];idx++){\n"
  "      youID = pLV[myID][idx];\n"
  "      if(youID>=iStartID && youID<=iEndID && !pVDTmp[youID]){\n"
  "        pVDTmp[youID]=(double)iDist;\n"
  "        pCheck[iCheckSz++]=youID;\n"
  "      }\n"
  "    }\n"
  "\n"
  "    iTmpSz = 0;  jdx=0;\n"
  "\n"
  "    iDist++;\n"
  "  \n"
  "    //this does a breadth-first search but avoids recursion\n"
  "    while(iCheckSz>0 && (iMaxDist==-1 || iDist<=iMaxDist)){\n"
  "      iTmpSz = 0;\n"
  "      for(idx=0;idx<iCheckSz;idx++){\n"
  "        youID=pCheck[idx];\n"
  "        for(jdx=0;jdx<pLen[youID];jdx++){\n"
  "          youKidID=pLV[youID][jdx];\n"
  "          if(youKidID >= iStartID && youKidID <=iEndID && !pVDTmp[youKidID]){ //found a new connection\n"
  "            pTmp[iTmpSz++] = youKidID; //save id of cell to search it's kids on next iteration\n"
  "            pVDTmp[youKidID]=(double)iDist;\n"
  "          }\n"
  "        }\n"
  "      }\n"
  "      iCheckSz = iTmpSz;\n"
  "      if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);\n"
  "      iDist++;\n"
  "    }\n"
  "\n"
  "    pVDTmp[myID]=0.0; // distance to self == 0.0\n"
  "    if((dgzt=gzmeandbl(pVDTmp,iStartID,iEndID))>0.0) pVD[myID]=dgzt;// save mean path length for given cell\n"
  "\n"
  "    memset(pVDTmp,0,sizeof(double)*iCells);\n"
  "  }\n"
  "  \n"
  "  free(pTmp);\n"
  "  if(pUse) free(pUse); \n"
  "  free(pCheck);\n"
  "  FreeListVec(&pList);  \n"
  "  free(pVDTmp);\n"
  "\n"
  "  if( verbose > 0 ) printf(\"\\n\");\n"
  "\n"
  "  return 1.0;\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ":* usage GetCCSubPop(adjlist,outvec,startids,endids[,subsamp])\n"
  ": computes clustering cofficient between sub-populations\n"
  ": adjlist == list of vectors specifying connectivity - adjacency list : from row -> to entry in column\n"
  ": outvec == vector of distances\n"
  ": startid == binary vector of ids of cells to start search from (from population)\n"
  ": endid   == binary vector of ids of cells to terminate search on (to population)\n"
  ": subsamp == perform calculation on ratio of cells btwn 0-1, default == 1\n"
  "FUNCTION GetCCSubPop () {\n"
  "  VERBATIM\n"
  "  ListVec* pList = AllocListVec(*hoc_objgetarg(1));\n"
  "  if(!pList){\n"
  "    printf(\"GetCCSubPop ERRA: problem initializing first arg!\\n\");\n"
  "    return 0.0;\n"
  "  }\n"
  " \n"
  "  int iCells = pList->isz; \n"
  "  if(iCells < 2){\n"
  "    printf(\"GetCCSubPop ERRB: size of List < 2 !\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  double** pLV = pList->pv;\n"
  "  unsigned int* pLen = pList->plen;\n"
  "\n"
  "  //init vector of distances to each cell , 0 == no path found\n"
  "  int* pNeighbors = (int*)calloc(iCells,sizeof(int));\n"
  "  int i = 0, iNeighbors = 0;\n"
  "  if(!pNeighbors){\n"
  "    printf(\"GetCCSubPop ERRE: out of memory!\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }  \n"
  "\n"
  "  //init vector of avg distances to each cell , 0 == no path found\n"
  "  double* pCC; \n"
  "  int iVecSz = vector_arg_px(2,&pCC);\n"
  "  if(!pCC || iVecSz < iCells){\n"
  "    printf(\"GetCCSubPop ERRE: arg 2 must be a Vector with size %d\\n\",iCells);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }  \n"
  "  memset(pCC,0,sizeof(double)*iVecSz);\n"
  "\n"
  "  double* pStart,  // bin vec of ids to search from \n"
  "          *pEnd;   // bin vec of ids to terminate search on\n"
  "\n"
  "  if( vector_arg_px(3,&pStart) < iCells || vector_arg_px(4,&pEnd) < iCells){\n"
  "    printf(\"GetCCSubPop ERRF: arg 3,4 must be Vectors with size >= %d\\n\",iCells);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "  double dSubsamp = ifarg(5)?*getarg(5):1.0;\n"
  "\n"
  "  unsigned int iSeed = ifarg(6)?(unsigned int)*getarg(6):INT_MAX-109754;\n"
  "\n"
  "  double* pUse = 0; \n"
  "  \n"
  "  if(dSubsamp<1.0){ //if using only a fraction of the cells\n"
  "     pUse = (double*)malloc(iCells*sizeof(double));\n"
  "     mcell_ran4(&iSeed, pUse, iCells, 1.0);\n"
  "  }\n"
  "\n"
  "  //get id of cell to find paths from\n"
  "  int myID;\n"
  "\n"
  "  int* pNeighborID = (int*)calloc(iCells,sizeof(int));\n"
  "\n"
  "  if( verbose > 0 ) printf(\"searching from id: \");\n"
  "\n"
  "  for(myID=0;myID<iCells;myID++) pCC[myID]=-1.0; //set invalid\n"
  "\n"
  "  for(myID=0;myID<iCells;myID++){\n"
  "\n"
  "    if(!pStart[myID]) continue;\n"
  "\n"
  "    if(verbose > 0 && myID%1000==0)printf(\"%d \",myID);\n"
  "\n"
  "    //only use dSubSamp fraction of cells, skip rest\n"
  "    if(pUse && pUse[myID]>=dSubsamp) continue;\n"
  "\n"
  "    int idx = 0, youID = 0, youKidID=0 , iNeighbors = 0;\n"
  "\n"
  "    //mark neighbors of distance == 1\n"
  "    for(idx=0;idx<pLen[myID];idx++){\n"
  "      youID = pLV[myID][idx];\n"
  "      if(pEnd[youID] && !pNeighbors[youID]){\n"
  "        pNeighbors[youID]=1;      \n"
  "        pNeighborID[iNeighbors++]=youID;\n"
  "      }\n"
  "    }\n"
  "\n"
  "    if(iNeighbors < 2){\n"
  "      for(i=0;i<iNeighbors;i++)pNeighbors[pNeighborID[i]]=0;\n"
  "      continue;\n"
  "    }\n"
  "\n"
  "    int iConns = 0 ; \n"
  "  \n"
  "    //this checks # of connections between neighbors of node\n"
  "    for(i=0;i<iNeighbors;i++){\n"
  "      if(!pNeighbors[pNeighborID[i]])continue;\n"
  "      youID=pNeighborID[i];\n"
  "      for(idx=0;idx<pLen[youID];idx++){\n"
  "        youKidID=pLV[youID][idx];\n"
  "        if(pEnd[youKidID] && pNeighbors[youKidID]){\n"
  "          iConns++;\n"
  "        }\n"
  "      }\n"
  "    }\n"
  "    pCC[myID]=(double)iConns/((double)iNeighbors*(iNeighbors-1));\n"
  "    for(i=0;i<iNeighbors;i++)pNeighbors[pNeighborID[i]]=0;\n"
  "  }\n"
  " \n"
  "  free(pNeighborID);\n"
  "  free(pNeighbors);\n"
  "  FreeListVec(&pList);\n"
  "  if(pUse)free(pUse);\n"
  "\n"
  "  if( verbose > 0 ) printf(\"\\n\");\n"
  "\n"
  "  return  1.0;\n"
  "\n"
  "  ENDVERBATIM\n"
  "}\n"
  ":* usage GetRecurCount(adjlist,outvec,fromids,thruids)\n"
  ": counts # of A -> B -> A patterns in adj adjacency list , using from ids as A\n"
  ": and thruids as B. fromids/thruids should have size of adjacency list and have a \n"
  ": 1 in index iff using that cell, same with thruids\n"
  "FUNCTION GetRecurCount () {\n"
  "  VERBATIM\n"
  "  ListVec* pList;\n"
  "  int iCells,iFromSz,iThruSz,idx,myID,youID,jdx,iCheckSz,*pVisited,*pCheck;\n"
  "  double **pLV,*pFrom,*pThru,*pR;\n"
  "\n"
  "  pList = AllocListVec(*hoc_objgetarg(1));\n"
  "  if(!pList){\n"
  "    printf(\"GetRecurCount ERRA: problem initializing first arg!\\n\");\n"
  "    return 0.0;\n"
  "  }\n"
  " \n"
  "  iCells = pList->isz; \n"
  "  if(iCells < 2){\n"
  "    printf(\"GetRecurCount ERRB: size of List < 2 !\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  pLV = pList->pv;\n"
  "  unsigned int* pLen = pList->plen;\n"
  "\n"
  "  pFrom=pThru=0;\n"
  "  iFromSz = vector_arg_px(3,&pFrom); iThruSz = vector_arg_px(4,&pThru);\n"
  "  \n"
  "  if( iFromSz <= 0 || iThruSz <= 0){\n"
  "    printf(\"GetRecurCount ERRF: arg 3,4 bad (fromsz,thrusz)=(%d,%d)\\n\",iFromSz,iThruSz);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  pVisited = (int*)calloc(iCells,sizeof(int));//which vertices already marked to have children expanded\n"
  "\n"
  "  pCheck = (int*)malloc(sizeof(int)*iCells);\n"
  "\n"
  "  pR = vector_newsize(vector_arg(2),iCells);\n"
  "  memset(pR,0,sizeof(double)*iCells); //zero out output first\n"
  "\n"
  "  for(myID=0;myID<iCells;myID++) {\n"
  "    if(!pFrom[myID]) continue;\n"
  "    iCheckSz = 0; \n"
  "    for(idx=0;idx<pLen[myID];idx++){//mark neighbors of distance == 1\n"
  "      youID = pLV[myID][idx];\n"
  "      if(!pThru[youID] || pVisited[youID]) continue;\n"
  "      pCheck[iCheckSz++]=youID;\n"
  "      pVisited[youID]=1;\n"
  "    }\n"
  "    for(idx=0;idx<iCheckSz;idx++) {\n"
  "      youID = pCheck[idx];\n"
  "      for(jdx=0;jdx<pLen[youID];jdx++) {\n"
  "        if(pLV[youID][jdx]==myID) pR[myID]++;\n"
  "      }\n"
  "    }\n"
  "    memset(pVisited,0,sizeof(int)*iCells);\n"
  "  }\n"
  "  \n"
  "\n"
  "  free(pCheck);\n"
  "  FreeListVec(&pList);  \n"
  "  free(pVisited);\n"
  "\n"
  "  if( verbose > 0) printf(\"\\n\");\n"
  "\n"
  "  return 1.0;\n"
  "\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ":* usage GetPairDist(adjlist,outvec,startid,endid[subsamp,seed])\n"
  ": computes distances between all pairs of vertices, self->self distance == distance of shortest loop\n"
  ": adjlist == list of vectors specifying connectivity - adjacency list : from row -> to entry in column\n"
  ": outvec == vector of distances from vertex i in outvec.x(i)\n"
  ": startid == first id to check\n"
  ": endid   == last id to check\n"
  "FUNCTION GetPairDist () {\n"
  "  VERBATIM\n"
  "  ListVec* pList = AllocListVec(*hoc_objgetarg(1));\n"
  "  if(!pList){\n"
  "    printf(\"GetPairDist ERRA: problem initializing first arg!\\n\");\n"
  "    return 0.0;\n"
  "  }\n"
  " \n"
  "  int iCells = pList->isz; \n"
  "  if(iCells < 2){\n"
  "    printf(\"GetPairDist ERRB: size of List < 2 !\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  double** pLV = pList->pv;\n"
  "  unsigned int* pLen = pList->plen;\n"
  "\n"
  "  double* pFrom = 0, *pTo = 0;\n"
  "  int iFromSz = vector_arg_px(3,&pFrom) , iToSz = vector_arg_px(4,&pTo);\n"
  "  \n"
  "  if( iFromSz <= 0 || iToSz <= 0){\n"
  "    printf(\"GetPairDist ERRF: arg 3,4 bad (fromsz,tosz)=(%d,%d)\\n\",iFromSz,iToSz);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  int iMinSz = iFromSz * iToSz;\n"
  "\n"
  "  //init vector of avg distances to each cell , 0 == no path found\n"
  "  double* pVD; \n"
  "  pVD = vector_newsize(vector_arg(2),iMinSz);\n"
  "  memset(pVD,0,sizeof(double)*iMinSz); //zero out output first\n"
  "\n"
  "  //init array of cells/neighbors to check\n"
  "  int* pCheck;\n"
  "  pCheck = (int*)malloc(sizeof(int)*iCells);\n"
  "  if(!pCheck){\n"
  "    printf(\"GetPairDist ERRG: out of memory!\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0, iTmpSz = 0, jdx = 0;\n"
  "\n"
  "  int* pTmp = (int*)calloc(iCells,sizeof(int)); \n"
  "\n"
  "  if( verbose > 0 ) printf(\"searching from id: \");\n"
  "\n"
  "  int myID , iOff = 0 , kdx = 0;\n"
  "\n"
  "  int* pVisited = (int*)calloc(iCells,sizeof(int)); //which vertices already marked to have children expanded\n"
  "  int* pUse = (int*)calloc(iCells,sizeof(int)); //which 'TO' vertices\n"
  "  int* pMap = (int*)calloc(iCells,sizeof(int)); //index of 'TO' vertices to output index\n"
  "  for(idx=0;idx<iToSz;idx++){\n"
  "    pUse[(int)pTo[idx]]=1;\n"
  "    pMap[(int)pTo[idx]]=idx;\n"
  "  }\n"
  "\n"
  "  for(kdx=0;kdx<iFromSz;kdx++,iOff+=iToSz){\n"
  "    myID=pFrom[kdx];\n"
  "    if(verbose > 0 && myID%100==0)printf(\"%d\\n\",myID);\n"
  "\n"
  "    iCheckSz = 0; idx = 0; iDist = 1; youID = 0; youKidID = 0;\n"
  "      \n"
  "    //mark neighbors of distance == 1\n"
  "    for(idx=0;idx<pLen[myID];idx++){\n"
  "      youID = pLV[myID][idx];\n"
  "      if(pUse[youID]) pVD[ iOff + pMap[youID]  ] = 1; //mark 1st degree neighbor distance as 1\n"
  "      if(!pVisited[youID]){ \n"
  "        pCheck[iCheckSz++]=youID;\n"
  "        pVisited[youID]=1;\n"
  "      }\n"
  "    }\n"
  "\n"
  "    iTmpSz = 0;  jdx=0;\n"
  "      \n"
  "    iDist++;\n"
  "  \n"
  "    //this does a breadth-first search but avoids recursion\n"
  "    while(iCheckSz>0){\n"
  "      iTmpSz = 0;\n"
  "      for(idx=0;idx<iCheckSz;idx++){\n"
  "        youID=pCheck[idx];\n"
  "        for(jdx=0;jdx<pLen[youID];jdx++){\n"
  "          youKidID=pLV[youID][jdx];\n"
  "          if(pUse[youKidID] && !pVD[iOff + pMap[youKidID]])\n"
  "            pVD[iOff + pMap[youKidID]] = iDist; \n"
  "          if(!pVisited[youKidID]){ //found a new connection\n"
  "            pTmp[iTmpSz++] = youKidID; //save id of cell to search it's kids on next iteration\n"
  "            pVisited[youKidID]=1;\n"
  "          }\n"
  "        }\n"
  "      }\n"
  "      iCheckSz = iTmpSz;\n"
  "      if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);\n"
  "      iDist++;\n"
  "    }\n"
  "    memset(pVisited,0,sizeof(int)*iCells);\n"
  "  }\n"
  "  \n"
  "  free(pTmp);\n"
  "  free(pCheck);\n"
  "  FreeListVec(&pList);  \n"
  "  free(pUse);\n"
  "  free(pMap);\n"
  "  free(pVisited);\n"
  "\n"
  "  if( verbose > 0) printf(\"\\n\");\n"
  "\n"
  "  return 1.0;\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ":* usage GetPathSubPop(adjlist,outvec,startids,endids[subsamp,loop,seed])\n"
  ": computes path lengths between sub-populations\n"
  ": adjlist == list of vectors specifying connectivity - adjacency list : from row -> to entry in column\n"
  ": outvec == vector of distances from vertex i in outvec.x(i)\n"
  ": startid == binary vector of ids of cells to start search from (from population)\n"
  ": endid   == binary vector of ids of cells to terminate search on (to population)\n"
  ": subsamp == perform calculation on ratio of cells btwn 0-1, default == 1\n"
  ": loop == check self-loops , default == 0\n"
  ": seed == random # seed when using subsampling\n"
  "FUNCTION GetPathSubPop () {\n"
  "  VERBATIM\n"
  "  ListVec* pList = AllocListVec(*hoc_objgetarg(1));\n"
  "  if(!pList){\n"
  "    printf(\"GetPathEV ERRA: problem initializing first arg!\\n\");\n"
  "    return 0.0;\n"
  "  }\n"
  " \n"
  "  int iCells = pList->isz; \n"
  "  if(iCells < 2){\n"
  "    printf(\"GetPathEV ERRB: size of List < 2 !\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  double** pLV = pList->pv;\n"
  "  unsigned int* pLen = pList->plen;\n"
  "\n"
  "  //init vector of avg distances to each cell , 0 == no path found\n"
  "  double* pVD; \n"
  "  int iVecSz = vector_arg_px(2,&pVD) , i = 0;\n"
  "  if(!pVD || iVecSz < iCells){\n"
  "    printf(\"GetPathEV ERRE: arg 2 must be a Vector with size %d\\n\",iCells);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }  \n"
  "  memset(pVD,0,sizeof(double)*iVecSz);\n"
  "\n"
  "  double* pStart,  // bin vec of ids to search from \n"
  "          *pEnd;   // bin vec of ids to terminate search on\n"
  "\n"
  "  if( vector_arg_px(3,&pStart) < iCells || vector_arg_px(4,&pEnd) < iCells){\n"
  "    printf(\"GetPathSubPop ERRF: arg 3,4 must be Vectors with size >= %d\\n\",iCells);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "  double dSubsamp = ifarg(5)?*getarg(5):1.0;\n"
  "\n"
  "  int bSelfLoop = ifarg(6)?(int)*getarg(6):0;\n"
  "\n"
  "  unsigned int iSeed = ifarg(7)?(unsigned int)*getarg(7):INT_MAX-109754;\n"
  "\n"
  "  //init array of cells/neighbors to check\n"
  "  int* pCheck = (int*)malloc(sizeof(int)*iCells);\n"
  "  if(!pCheck){\n"
  "    printf(\"GetPathEV ERRG: out of memory!\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0, iTmpSz = 0, jdx = 0;\n"
  "\n"
  "  double  dgzt = 0.0; \n"
  "  int* pTmp = 0;\n"
  "  double* pUse = 0; \n"
  "  \n"
  "  if(dSubsamp<1.0){ //if using only a fraction of the cells\n"
  "     pUse = (double*)malloc(iCells*sizeof(double));\n"
  "     mcell_ran4(&iSeed, pUse, iCells, 1.0);\n"
  "  }\n"
  "\n"
  "  pTmp = (int*)calloc(iCells,sizeof(int)); \n"
  "\n"
  "  if( verbose > 0 ) printf(\"searching from id: \");\n"
  "\n"
  "  int* pVDTmp = (int*)calloc(iCells,sizeof(int)) , myID;\n"
  "\n"
  "  for(myID=0;myID<iCells;myID++){\n"
  "\n"
  "    if(!pStart[myID]) continue;\n"
  "\n"
  "    if(verbose > 0 && myID%1000==0)printf(\"%d \",myID); \n"
  "\n"
  "    //only use dSubSamp fraction of cells, skip rest\n"
  "    if(pUse && pUse[myID]>=dSubsamp) continue;\n"
  "\n"
  "    unsigned long int iSelfLoopDist = LONG_MAX;\n"
  "    int bFindThisSelfLoop = bSelfLoop && pEnd[myID]; // search for self loop for this vertex?\n"
  "\n"
  "    iCheckSz = 0; idx = 0; iDist = 1; youID = 0; youKidID = 0;\n"
  "\n"
  "    pVDTmp[myID]=1;\n"
  "\n"
  "    //mark neighbors of distance == 1\n"
  "    for(idx=0;idx<pLen[myID];idx++){\n"
  "      youID = pLV[myID][idx];\n"
  "      if(bFindThisSelfLoop && youID==myID && iDist<iSelfLoopDist) iSelfLoopDist = iDist; //found a self-loop? \n"
  "      if(!pVDTmp[youID]){\n"
  "        pVDTmp[youID]=iDist;\n"
  "        pCheck[iCheckSz++]=youID;\n"
  "      }\n"
  "    }\n"
  "\n"
  "    iTmpSz = 0;  jdx=0;\n"
  "\n"
  "    iDist++;\n"
  "  \n"
  "    //this does a breadth-first search but avoids recursion\n"
  "    while(iCheckSz>0){\n"
  "      iTmpSz = 0;\n"
  "      for(idx=0;idx<iCheckSz;idx++){\n"
  "        youID=pCheck[idx];\n"
  "        for(jdx=0;jdx<pLen[youID];jdx++){\n"
  "          youKidID=pLV[youID][jdx];\n"
  "          if(bFindThisSelfLoop && youKidID==myID && iDist<iSelfLoopDist) iSelfLoopDist = iDist; //found a self-loop? \n"
  "          if(!pVDTmp[youKidID]){ //found a new connection\n"
  "            pTmp[iTmpSz++] = youKidID; //save id of cell to search it's kids on next iteration\n"
  "            pVDTmp[youKidID]=iDist;\n"
  "          }\n"
  "        }\n"
  "      }\n"
  "      iCheckSz = iTmpSz;\n"
  "      if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);\n"
  "      iDist++;\n"
  "    }\n"
  "\n"
  "    if(bFindThisSelfLoop && iSelfLoopDist<LONG_MAX){//if checking for this vertex's self-loop dist. and found a self-loop\n"
  "      pVDTmp[myID] = iSelfLoopDist;\n"
  "    } else {\n"
  "      pVDTmp[myID]=0; // distance to self == 0.0\n"
  "    }\n"
  "    pVD[myID] = 0.0;\n"
  "    int N = 0; //take average path length (+ self-loop length if needed) from myID to pEnd cells\n"
  "    for(idx=0;idx<iCells;idx++){\n"
  "      if(pEnd[idx] && pVDTmp[idx]){\n"
  "        pVD[myID] += pVDTmp[idx];\n"
  "        N++;\n"
  "      }\n"
  "    }\n"
  "\n"
  "    if(N) pVD[myID] /= (double) N; // save mean path (and maybe self-loop) length for given cell\n"
  "\n"
  "    memset(pVDTmp,0,sizeof(int)*iCells);\n"
  "  }\n"
  "  \n"
  "  free(pTmp);\n"
  "  if(pUse) free(pUse); \n"
  "  free(pCheck);\n"
  "  FreeListVec(&pList);  \n"
  "  free(pVDTmp);\n"
  "\n"
  "  if( verbose > 0 ) printf(\"\\n\");\n"
  "\n"
  "  return 1.0;\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ":* usage GetLoopLength(adjlist,outvec,loopids,thruids[,subsamp,seed])\n"
  ": computes distance to loop back to each node\n"
  ": adjlist == list of vectors specifying connectivity - adjacency list : from row -> to entry in column\n"
  ": outvec == vector of distances\n"
  ": loopids == binary vector of ids of cells to start/end search from/to\n"
  ": thruids == binary vector of ids of cells thru which loop can pass\n"
  ": subsamp == perform calculation on ratio of cells btwn 0-1, default == 1\n"
  ": seed == random # seed when using subsampling\n"
  "FUNCTION GetLoopLength () {\n"
  "  VERBATIM\n"
  "  ListVec* pList = AllocListVec(*hoc_objgetarg(1));\n"
  "  if(!pList){\n"
  "    printf(\"GetLoopLength ERRA: problem initializing first arg!\\n\");\n"
  "    return 0.0;\n"
  "  }\n"
  " \n"
  "  int iCells = pList->isz; \n"
  "  if(iCells < 2){\n"
  "    printf(\"GetLoopLength ERRB: size of List < 2 !\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  double** pLV = pList->pv;\n"
  "  unsigned int* pLen = pList->plen;\n"
  "\n"
  "  //init vector of avg distances to each cell , 0 == no path found\n"
  "  double* pVD; \n"
  "  int iVecSz = vector_arg_px(2,&pVD) , i = 0;\n"
  "  if(!pVD || iVecSz < iCells){\n"
  "    printf(\"GetLoopLength ERRE: arg 2 must be a Vector with size %d\\n\",iCells);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }  \n"
  "  memset(pVD,0,sizeof(double)*iVecSz);//init to 0\n"
  "\n"
  "  double* pLoop,  // bin vec of ids to search from \n"
  "          *pThru;   // bin vec of ids to terminate search on\n"
  "\n"
  "  if( vector_arg_px(3,&pLoop) < iCells || vector_arg_px(4,&pThru) < iCells){\n"
  "    printf(\"GetLoopLength ERRF: arg 3,4 must be Vectors with size >= %d\\n\",iCells);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "  double dSubsamp = ifarg(5)?*getarg(5):1.0;\n"
  "\n"
  "  unsigned int iSeed = ifarg(6)?(unsigned int)*getarg(6):INT_MAX-109754;\n"
  "\n"
  "  //init array of cells/neighbors to check\n"
  "  int* pCheck = (int*)malloc(sizeof(int)*iCells);\n"
  "  if(!pCheck){\n"
  "    printf(\"GetLoopLength ERRG: out of memory!\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0, iTmpSz = 0, jdx = 0;\n"
  "\n"
  "  double  dgzt = 0.0; \n"
  "  int* pTmp = 0 , found = 0;\n"
  "  double* pUse = 0; \n"
  "  \n"
  "  if(dSubsamp<1.0){ //if using only a fraction of the cells\n"
  "     pUse = (double*)malloc(iCells*sizeof(double));\n"
  "     mcell_ran4(&iSeed, pUse, iCells, 1.0);\n"
  "  }\n"
  "\n"
  "  pTmp = (int*)calloc(iCells,sizeof(int)); \n"
  "\n"
  "  if( verbose > 0 ) printf(\"searching loops from id: \");\n"
  "\n"
  "  int* pVDTmp = (int*)calloc(iCells,sizeof(int)) , myID;\n"
  "\n"
  "  for(myID=0;myID<iCells;myID++){\n"
  "\n"
  "    if(!pLoop[myID]) continue;\n"
  "\n"
  "    if(verbose > 0 && myID%1000==0)printf(\"%d \",myID); \n"
  "\n"
  "    //only use dSubSamp fraction of cells, skip rest\n"
  "    if(pUse && pUse[myID]>=dSubsamp) continue;\n"
  "\n"
  "    iCheckSz = 0; idx = 0; iDist = 1; youID = 0; youKidID = 0; found = 0;\n"
  "\n"
  "    pVDTmp[myID]=1;\n"
  "\n"
  "    //mark neighbors of distance == 1\n"
  "    for(idx=0;idx<pLen[myID];idx++){\n"
  "      youID = pLV[myID][idx];\n"
  "      if(youID==myID) {\n"
  "        found = 1;\n"
  "        pVD[myID]=iDist;\n"
  "        iCheckSz=0;\n"
  "        break;\n"
  "      }\n"
  "      if(pThru[youID] && !pVDTmp[youID]){\n"
  "        pVDTmp[youID]=iDist;\n"
  "        pCheck[iCheckSz++]=youID;\n"
  "      }\n"
  "    }\n"
  "\n"
  "    iTmpSz = 0;  jdx=0;\n"
  "\n"
  "    iDist++;\n"
  "  \n"
  "    //this does a breadth-first search but avoids recursion\n"
  "    while(iCheckSz>0){\n"
  "      iTmpSz = 0;\n"
  "      for(idx=0;idx<iCheckSz;idx++){\n"
  "        youID=pCheck[idx];\n"
  "        for(jdx=0;jdx<pLen[youID];jdx++){\n"
  "          youKidID=pLV[youID][jdx];\n"
  "          if(youKidID==myID){\n"
  "            pVD[myID]=iDist;\n"
  "            found = 1;\n"
  "            break;\n"
  "          }\n"
  "          if(pThru[youKidID] && !pVDTmp[youKidID]){ //found a new connection\n"
  "            pTmp[iTmpSz++] = youKidID; //save id of cell to search it's kids on next iteration\n"
  "            pVDTmp[youKidID]=iDist;\n"
  "          }\n"
  "        }\n"
  "      }\n"
  "      if(found) break;\n"
  "      iCheckSz = iTmpSz;\n"
  "      if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);\n"
  "      iDist++;\n"
  "    }\n"
  "    memset(pVDTmp,0,sizeof(int)*iCells);\n"
  "  }\n"
  "  \n"
  "  free(pTmp);\n"
  "  if(pUse) free(pUse); \n"
  "  free(pCheck);\n"
  "  FreeListVec(&pList);  \n"
  "  free(pVDTmp);\n"
  "\n"
  "  if( verbose > 0 ) printf(\"\\n\");\n"
  "\n"
  "  return 1.0;\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ":* usage GetPathEV(adjlist,outvec,myid,[startid,endid,maxdist])\n"
  ": adjlist == list of vectors specifying connectivity - adjacency list : from row -> to entry in column\n"
  ": outvec == vector of distances\n"
  ": myid == id of cell to start search from\n"
  ": startid == min id of cells search can terminate on or go through\n"
  ": endid   == max  '    '   '  '   '  '  '  '  ' '  '  '  '  '  ' \n"
  "FUNCTION GetPathEV () {\n"
  "  VERBATIM\n"
  "  ListVec* pList = AllocListVec(*hoc_objgetarg(1));\n"
  "  if(!pList){\n"
  "    printf(\"GetPathEV ERRA: problem initializing first arg!\\n\");\n"
  "    return 0.0;\n"
  "  }\n"
  " \n"
  "  int iCells = pList->isz; \n"
  "  if(iCells < 2){\n"
  "    printf(\"GetPathEV ERRB: size of List < 2 !\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  double** pLV = pList->pv;\n"
  "  unsigned int* pLen = pList->plen;\n"
  "\n"
  "  //init vector of distances to each cell , 0 == no path found\n"
  "  double* pVD; \n"
  "  int iVecSz = vector_arg_px(2,&pVD) , i = 0;\n"
  "  if(!pVD || iVecSz < iCells){\n"
  "    printf(\"GetPathEV ERRE: arg 2 must be a Vector with size %d\\n\",iCells);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }  \n"
  "  memset(pVD,0,sizeof(double)*iVecSz);//init to 0\n"
  "\n"
  "  //get id of cell to find paths from\n"
  "  int myID = (int) *getarg(3);\n"
  "  if(myID < 0 || myID >= iCells){\n"
  "    printf(\"GetPathEV ERRF: invalid id = %d\\n\",myID);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  //start/end id of cells to find path to\n"
  "  int iStartID = ifarg(4) ? (int)*getarg(4) : 0,\n"
  "      iEndID = ifarg(5) ? (int)*getarg(5) : iCells - 1,\n"
  "      iMaxDist = ifarg(6)? (int)*getarg(6): -1;\n"
  "\n"
  "  if(iStartID < 0 || iStartID >= iCells ||\n"
  "     iEndID < 0 || iEndID >= iCells ||\n"
  "     iStartID >= iEndID){\n"
  "       printf(\"GetPathEV ERRH: invalid ids start=%d end=%d numcells=%d\\n\",iStartID,iEndID,iCells);\n"
  "       FreeListVec(&pList);\n"
  "       return 0.0;\n"
  "     }\n"
  "\n"
  "  //check max distance\n"
  "  if(iMaxDist==0){\n"
  "    printf(\"GetPathEV ERRI: invalid maxdist=%d\\n\",iMaxDist);\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "\n"
  "  //init array of cells/neighbors to check\n"
  "  int* pCheck = (int*)malloc(sizeof(int)*iCells);\n"
  "  if(!pCheck){\n"
  "    printf(\"GetPathEV ERRG: out of memory!\\n\");\n"
  "    FreeListVec(&pList);\n"
  "    return 0.0;\n"
  "  }\n"
  "  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0;\n"
  "\n"
  "  pVD[myID]=1;\n"
  "\n"
  "  //mark neighbors of distance == 1\n"
  "  for(idx=0;idx<pLen[myID];idx++){\n"
  "    youID = pLV[myID][idx];\n"
  "    if(youID>=iStartID && youID<=iEndID && !pVD[youID]){\n"
  "      pVD[youID]=(double)iDist;\n"
  "      pCheck[iCheckSz++]=youID;\n"
  "    }\n"
  "  }\n"
  "\n"
  "  int* pTmp = (int*)malloc(sizeof(int)*iCells);\n"
  "  int iTmpSz = 0 , jdx=0;\n"
  "\n"
  "  iDist++;\n"
  "  \n"
  "  //this does a breadth-first search but avoids deep nesting of recursive version\n"
  "  while(iCheckSz>0 && (iMaxDist==-1 || iDist<=iMaxDist)){\n"
  "    iTmpSz = 0;\n"
  "    for(idx=0;idx<iCheckSz;idx++){\n"
  "      youID=pCheck[idx];\n"
  "      for(jdx=0;jdx<pLen[youID];jdx++){\n"
  "        youKidID=pLV[youID][jdx];\n"
  "        if(youKidID >= iStartID && youKidID <=iEndID && !pVD[youKidID]){ //found a new connection\n"
  "          pTmp[iTmpSz++] = youKidID; //save id of cell to search it's kids on next iteration\n"
  "          pVD[youKidID]=(double)iDist;\n"
  "        }\n"
  "      }\n"
  "    }\n"
  "    iCheckSz = iTmpSz;\n"
  "    if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);\n"
  "    iDist++;\n"
  "  }\n"
  "\n"
  "  pVD[myID]=0.0;\n"
  " \n"
  "  free(pCheck);\n"
  "  free(pTmp);\n"
  "  FreeListVec(&pList);\n"
  "\n"
  "  return 1.0;\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ":* FUNCTION Factorial()\n"
  "FUNCTION Factorial () {\n"
  "  VERBATIM\n"
  "  double N = (int)*getarg(1) , i = 0.0;\n"
  "  double val = 1.0;\n"
  "  if(N<=1) return 1.0;\n"
  "  if(N>=171){\n"
  "    double PI=3.1415926535897932384626433832795;\n"
  "    double E=2.71828183;\n"
  "    val=sqrt(2*PI*N)*(pow(N,N)/pow(E,N));\n"
  "  } else {\n"
  "    for(i=2.0;i<=N;i++) val*=i;\n"
  "  }\n"
  "  return (double) val;  \n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ":* FUNCTION perm()\n"
  ":count # of permutations from set of N elements with R selections\n"
  "FUNCTION perm () {\n"
  "  VERBATIM\n"
  "  if(ifarg(3)){\n"
  "    double N = (int)*getarg(1);\n"
  "    double R = (int)*getarg(2);\n"
  "    double b = *getarg(3);\n"
  "    double val = N/b;\n"
  "    int i = 0;\n"
  "    for(i=1;i<R;i++){\n"
  "      N--;\n"
  "      val*=(N/b);\n"
  "    }\n"
  "    return val;\n"
  "  } else {\n"
  "    int N = (int)*getarg(1);\n"
  "    int R = (int)*getarg(2);\n"
  "    int val = N;\n"
  "    int i = 0;\n"
  "    for(i=1;i<R;i++){\n"
  "      N--;\n"
  "      val*=N;\n"
  "    }\n"
  "    return (double)val;\n"
  "  }\n"
  "  ENDVERBATIM\n"
  "}\n"
  "\n"
  ":* install_intfsw\n"
  "PROCEDURE install () {\n"
  " if(INSTALLED==1){\n"
  "   printf(\"Already installed $Id: intfsw.mod,v 1.50 2009/02/26 18:24:34 samn Exp $ \\n\")\n"
  " } else {\n"
  " INSTALLED=1\n"
  " VERBATIM\n"
  " install_vector_method(\"gzmean\" ,gzmean);\n"
  " install_vector_method(\"nnmean\" ,nnmean);\n"
  " install_vector_method(\"copynz\" ,copynz);\n"
  " ENDVERBATIM\n"
  " printf(\"Installed $Id: intfsw.mod,v 1.50 2009/02/26 18:24:34 samn Exp $ \\n\")\n"
  " }\n"
  "}\n"
  ;
    hoc_reg_nmodl_filename(mech_type, nmodl_filename);
    hoc_reg_nmodl_text(mech_type, nmodl_file_text);
}
#endif
