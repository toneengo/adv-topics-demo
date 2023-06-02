#ifndef INTERVAL_CPP
#define INTERVAL_CPP

#include <vector>
#include <iostream>
#include <functional>
#include <math.h>

// Declaring outside of class to give access to include files
using ddfunc = std::function<double(double)>;

class Interval {
    public:
        enum Score {
            unassigned = -1,
            interesting = 0,
            ambiguous = 1,
            boring = 2
        };

    private:
        double start_;
        double end_;
        Score score_ = unassigned;

        std::pair<int, int> winner_;
        std::pair<int, int> transition_;

        // @return: The index of the two winning functions at a given point
        std::pair<int, int> sortFuncs(std::vector<ddfunc> const &funcs, double x, int &fcalls) {
            std::vector<std::pair<double, int>> order;

            for (int i = 0; i < funcs.size(); i++) {
                order.push_back({funcs[i](x), i}); // Keep track of both the call  the index
                ++fcalls;
            }

            std::sort(order.begin(), order.end(), std::greater<>());

            return {order[0].second, order[1].second};
        }

    public:
        // * Constructors
        Interval() {
            start_ = 0;
            end_ = 1;
        }
        Interval(double start, double end) {
            start_ = start;
            end_ = end;
        }

        // * Accessors
        double& start() { return start_; }
        double& end() { return end_; }
        Score& score() { return score_; }
        std::pair<int, int>& winner() { return winner_ ; }
        std::pair<int, int>& transition() { return transition_ ; }

        /*
            * Comparison function for intervals

            * @return true if intervals are within tolerance value, else false
        */
        static bool comp(Interval interval1, Interval interval2, double tol=1e-6) {
            if (
                (abs(interval1.start() - interval2.start()) < tol) && 
                (abs(interval1.end() - interval2.end()) < tol)
            ) {
                return true;
            }
            
            return false;
        };
        
        bool comp(double start, double end, double tol=1.48e-8) {
            if (abs(start_ - start) < tol && abs(end_ - end) < tol) {
                return true;
            }

            return false;
        }

        // String representation
        std::string toString() {
            return ("(" + std::to_string(start_) + ", " + std::to_string(end_) + ")");
        }

        /*
            * @return Interval score (interesting, ambiguous or boring) as an integer
        */ 
        Score linearInterpolationScorer(
            std::function<double(double)> f, 
            std::function<double(double)> fprime, 
            int &fcalls,
            double tol=1.48e-6
        ) {
            double diff = end_ - start_;
            std::pair<double, double> ans;
            std::vector<double> mStore;
            std::vector<double> bStore;
            int signCheck = 1;

            for (double point : {start_, end_}) {
                double x1 = f(point);
                double fder = fprime(point);
                
                signCheck *= x1;
                
                if (fder == 0) {
                    mStore.push_back(tol);
                }
                else {
                    mStore.push_back(fder);
                }

                bStore.push_back(x1 - point * fder);
            }

            fcalls += 4;

            // y = mx + b
            // mx + b = 0
            // x = -(b/m)
            ans.first = -1 * (bStore[0] / mStore[0]);
            ans.second = -1 * (bStore[1] / mStore[1]);

            // Sign change case
            if (signCheck < 0) {
                score_ = ambiguous;
                return ambiguous;
            }
            // Interesting case
            else if ((start_ < ans.first && ans.first < end_) && (start_ < ans.second && ans.second < end_)) {
                score_ = interesting;
                return interesting;
            }
            // Ambiguous case with extended interval
            else if ( 
                ((start_ < ans.first && ans.first < end_) || (start_ < ans.second && ans.second < end_)) ||
                ((start_ - diff < ans.first && ans.first < end_ + diff) || (start_ - diff < ans.second && ans.second < end_ + diff))
            ) {
                score_ = ambiguous;
                return ambiguous;
            }

            // Boring case
            score_ = boring;
            return boring;
        }

        Score multiDomainScorer(
            std::vector<ddfunc> const &funcs,
            int &fcalls,
            double tol=1e-6
        ) {
            if (score_ == unassigned) {
                std::pair<int, int> startWinners = sortFuncs(funcs, start_, fcalls);
                std::pair<int, int> endWinners = sortFuncs(funcs, end_, fcalls);

                if (startWinners.first == endWinners.first) {
                    score_ = boring;
                    return boring;
                }
                else if (startWinners.first == endWinners.second && startWinners.second == endWinners.first) {
                    score_ = interesting;
                    winner_ = {startWinners.first, startWinners.second};
                    transition_ = {startWinners.first, endWinners.first};
                    return interesting;
                }
                else {
                    score_ = ambiguous;
                    return ambiguous;
                }
            }

            return score_;
        }
};

#endif