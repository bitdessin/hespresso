#define TMB_LIB_INIT R_init_hespresso

#include <TMB.hpp>
#include "hobit.h"

template<class Type>
Type objective_function<Type>::operator()() {
    return hobit_mle(this);
}
