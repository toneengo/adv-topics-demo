#ifndef INTERVAL_NODE_CPP
#define INTERVAL_NODE_CPP

#include <memory>

#include "Interval.cpp"

class IntervalNode {
    public:        
        using SharedPtr = std::shared_ptr<IntervalNode>;
        using SharedVec = std::vector<std::shared_ptr<IntervalNode>>;

    private:
        SharedPtr left_;
        SharedPtr right_;
        Interval interval_;

    public:
        // * Constructors
        IntervalNode() { interval_ = Interval(); }
        // IntervalNode(Interval interval) { interval_ = interval; }
        IntervalNode(double start, double end) { interval_ = Interval(start, end); }

        bool isLeaf() {
            if (left_ == nullptr && right_ == nullptr) {
                return true;
            }

            return false;
        }

        // * Accessors
        SharedPtr& left() {return left_; }
        SharedPtr& right() { return right_; } 
        Interval& interval() { return interval_; }
};

#endif