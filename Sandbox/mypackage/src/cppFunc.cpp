#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame rescale_to_coverage(NumericVector cdf, double coverage, NumericVector Mid_Week) {

            double cov = coverage;
            int len_cdf = cdf.size();
            double max_cdf = cdf[len_cdf-1];
            double x = max_cdf/cov;
            NumericVector rescaled_cdf(len_cdf);

            NumericVector pdf_list (len_cdf);
            pdf_list[0] = 0;

            rescaled_cdf = cdf/x;
            pdf_list[Rcpp::Range(1, 259)] = rescaled_cdf[Rcpp::Range(1, 259)] - rescaled_cdf[Rcpp::Range(0, 258)];

            DataFrame df = DataFrame::create(Named("Mid_Week") = Mid_Week, _["RESCALED_CDF"] = rescaled_cdf, _["PDF"] = pdf_list);
            return df;

            }

// [[Rcpp::export]]
NumericMatrix incremental_cov_leakage(NumericMatrix x, NumericMatrix y) {
              int nrow = x.nrow()-1;
                        int ncol = x.ncol();

                        for (int w = 1; w <= nrow; w++) {
                        int t_from = w + 1;
                        for (int t = t_from; t < ncol; t++) {
                        int ww = w-1;
                        int tt = t-1;
                        NumericMatrix subtract_range1 = x(Range(0,ww), Range(tt,tt));
                        double sum_range1 = sum(subtract_range1);

                        NumericMatrix subtract_range2 = x(Range(0,ww), Range(t,t));
                        double sum_range2 = sum(subtract_range2);

                        double difference = sum_range1 - sum_range2;

                        x(w,t) = x(w,tt) - (y(t,t) - difference);
                        if (x(w,t) < 0) {
                        x(w,t) = 0;
                        }
                        }
                        }
                        return x;
                        }

// [[Rcpp::export]]
NumericMatrix impact_time_week(NumericMatrix x, NumericMatrix y, NumericMatrix z, NumericVector a, NumericVector b, NumericVector c) {
  int nrow = z.nrow();
            int ncol = z.ncol();
            NumericMatrix impact_tracked( nrow , ncol );

            for (int w = 0; w <= nrow; w++) {
            int k = w - 1;
            for (int t = w; t < ncol; t++) {
            impact_tracked(w,t) = z(w,t)*c[t-k] + y(w,t)*b[t-k] + x(w,t)*a[t-k];
            }
            }
            return impact_tracked;
            }

// [[Rcpp::export]]
NumericMatrix cumulative_covC(NumericMatrix x) {
            int nrow = x.nrow();
            int ncol = x.ncol();

            for (int i = 1; i < nrow; i++) {
            for (int j = 0; j < ncol; j++) {
            double x2 = x(i, j);
            double x1 = x(i-1,j);
            if ( x2 > x1 ){
            NumericVector v = {x2, x1};
            x(i,j) = min(v);
            }
            }
            }
            return x;
            }

// [[Rcpp::export]]
NumericMatrix incremental_covC(NumericMatrix x) {
            NumericMatrix m(3, 260);
            int nrow = x.nrow();
            int ncol = x.ncol()-1;

            for (int i = 1; i < nrow; i++) {
            for (int j = 1; j < ncol; j++) {

            if ( x(i-1, j) == x(i-1,j-1) ){
            m(i-1,j) = 0;
            } else if ( x(i, j) < x(i-1, j-1) ){
            m(i-1,j) = x(i-1,j) - x(i-1,j-1);
            } else {
            m(i-1,j) = x(i-1,j) - x(i,j);
            }
            }
            }
            return m;

            }
