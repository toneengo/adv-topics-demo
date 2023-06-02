#ifndef INTERVAL_TREE_CPP
#define INTERVAL_TREE_CPP

#include <queue>

#include "IntervalNode.cpp"

// Declaring outside of class to give access to include files
using ddpair = std::pair<double, double>;

class IntervalTree {
    public:
        using SharedPtr = std::shared_ptr<IntervalTree>;

    private:
        IntervalNode::SharedPtr root_;

        /*
            * Splits a root into a subtree containing divided intervals
            
            * @param queue: Reference to a queue of intervals to divide
            * @param subtree: Reference to a vector of child nodes
            * @param subtree: Reference to the current subtree being divided
            * @param start: Start of the interval
            * @param end: End of the interval
            * @param level: Levels needed to traverse
            * @param lim: Tolerance limit
              
            * @return vector of references to subtree
        */
        void internalBFS(
            std::queue<ddpair>& intervalQueue,
            std::vector<IntervalNode::SharedPtr>& subtree, 
            IntervalNode::SharedPtr subtreeRoot,
            double start,
            double end,
            int levels, 
            double lim
        ) {
            // if interval is smaller than the limit or sufficient levels have been traversed
            if (0.5 * (end - start) > lim && levels > 0) {
                if (end == subtreeRoot->interval().end()) {
                    levels = levels - 1;
                }
                
                // get next node in the queue (starts as subtree root)
                IntervalNode::SharedPtr node = getIntervalNode(subtreeRoot, start, end);

                double middle = 0.5 * (start + end);

                if (isLeaf(node)) {
                    node->left() = std::make_shared<IntervalNode>(start, middle);
                    node->right() = std::make_shared<IntervalNode>(middle, end);
                }

                subtree.push_back(node->left());
                subtree.push_back(node->right());

                intervalQueue.push({start, middle});
                intervalQueue.push({middle, end});

                double newStart;
                double newEnd;

                std::tie(newStart, newEnd) = intervalQueue.front();

                intervalQueue.pop();

                internalBFS(
                    intervalQueue, 
                    subtree, 
                    subtreeRoot, 
                    newStart,
                    newEnd, 
                    levels, 
                    lim
                );
            }
        }        

    public:
        // * Constructors
        IntervalTree() {
            root_ = std::make_shared<IntervalNode>();
        };

        IntervalTree(double start, double end) {
            root_ = std::make_shared<IntervalNode>(start, end);
        }

        // * Accessors
        IntervalNode::SharedPtr& root() { return root_; }

        int size() { return size(root_); }
        
        static int size(IntervalNode::SharedPtr node) {
            if (node == nullptr) {
                return 0;
            }

            return (1 + size(node->left()) + size(node->right()));
        }

        /*
            * @return Shared pointer to node containing matching interval 
            * or nullptr if node does not exist
        */
        IntervalNode::SharedPtr getIntervalNode(IntervalNode::SharedPtr node, double start, double end) {
            if (node == nullptr) {
                return nullptr;
            }

            if (node->interval().comp(start, end)) {
                return node;
            }

            IntervalNode::SharedPtr left = getIntervalNode(node->left(), start, end);

            if (left != nullptr) {
                return left; 
            }

            IntervalNode::SharedPtr right = getIntervalNode(node->right(), start, end);

            if (right != nullptr) {
                return right; 
            }

            return nullptr;
        }

        bool isLeaf(IntervalNode::SharedPtr node) {
            if (node->left() == nullptr && node->right() == nullptr) {
                return true;
            }

            return false;
        }

        /*
            * @return A vector containing pointers to all subtree nodes
            * created during the bfs
        */
        IntervalNode::SharedVec intervalLevels(IntervalNode::SharedPtr& subtreeRoot, int levels, double lim=1e-6) {
            std::queue<ddpair> intervalQueue = {};            
            IntervalNode::SharedVec subtree = {};
            double start = subtreeRoot->interval().start();
            double end = subtreeRoot->interval().end();

            internalBFS(intervalQueue, subtree, subtreeRoot, start, end, levels, lim);

            return subtree;
        }
};

#endif