#pragma once
#include <cstdio>
#include "input.h"
#include "chain.h"
#include "MC.h"

// Writes observables to <out_prefix>_obs.dat and configurations to <out_prefix>_conf.xyz.
class Logger {
public:
    explicit Logger(const Input& in);
    ~Logger();
    void header();
    void sample(long sweep, const Chain& chain, const MC& mc);
    void dump_config(long sweep, const Chain& chain);
    void summary(const Chain& chain, const MC& mc);

private:
    const Input& in_;
    std::FILE* obs_  = nullptr;
    std::FILE* conf_ = nullptr;

    // running averages for the summary
    long   n_samples_   = 0;
    double sum_helix_   = 0.0;
    double sum_ebond_   = 0.0;
    double sum_ree_     = 0.0;
    double sum_rg2_     = 0.0;
};
