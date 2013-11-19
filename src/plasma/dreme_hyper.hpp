
#ifndef DREME_HYPER_HPP

/*  Return log of hypergeometric pvalue of #pos >= p
    p = positive successes
    P = positives
    n = negative successes
    N = negatives
    log_pthresh = short-circuit if p will be greater
    */
double getLogFETPvalue(double p, double P, double n, double N, double log_pthresh);

#endif

