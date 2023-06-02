#ifndef MULTI_DOMAIN_ROOT_FINDER_CPP
#define MULTI_DOMAIN_ROOT_FINDER_CPP

#include <map>

#include "linearRootFinder.cpp"

int fact(int n) {
    if (n < 1) {
        return 1;
    }

    return n * fact(n - 1);
}

// n! / k!(n-k)!
int binom(int n, int k) {
    return fact(n) / (fact(k) * fact(n - k));
}

ddfunc bernstein(int v, int degree) {
    if (v < 0 || v > degree) {
        throw (std::domain_error("The index of the function must be above 0 and less than or equal to the degree"));
    }

    return [=](double x) { return binom(degree, v) * pow(x, v) * pow(1 - x, degree - v); };
}

ddfunc bernsteinPrime(int v, int degree) {
    return [=](double x) { return binom(degree, v) * ((v * pow(1 - x, degree - v) * pow(x, v - 1)) - ((degree - v) * pow(1 - x, degree - v - 1) * pow(x, v))); };
}

std::vector<ddfunc> bernsteinGenerator(int degree) {
    std::vector<ddfunc> fs;
    for (int i = 0; i < degree + 1; i++) {
        fs.push_back(bernstein(i, degree));
    }
    return fs;
}

std::vector<ddfunc> bernsteinPrimeGenerator(int degree) {
    std::vector<ddfunc> fprimes;
    for (int i = 0; i < degree + 1; i++) {
        fprimes.push_back(bernsteinPrime(i, degree));
    }
    return fprimes;
}

int highestValueEquation(std::vector<ddfunc> const &fs, double x) {
    std::vector<std::pair<double, int>> order;

    for (int i = 0; i < fs.size(); i++) {
        order.push_back({fs[i](x), i}); // Keep track of both the call and the index
    }

    std::sort(order.begin(), order.end(), std::greater<>());

    return order[0].second;
}

std::pair<std::vector<int>, std::vector<double>> dfsMultiDomain(
    double start,
    double end,
    std::vector<ddfunc> const fs,
    std::vector<ddfunc> const fprimes,
    int verificationLevels = 2,
    int subdivisionLevels = 1,
    double dfsTerminationDistance= 1 / pow(2, DFS_TERMINATION_LIMIT), // Changed to match reference
    int maxNewtonIterations = 250
) {
    std::vector<std::pair<double, int>> switchPointStore = {
        { end, highestValueEquation(fs, end) }
    };

    double x1;
    double x2;
    int fcalls = 0;

    std::stack<IntervalNode::SharedPtr> verificationStack;
    std::stack<IntervalNode::SharedPtr> ambiguousStack;

    IntervalNode::SharedVec visited;

    IntervalTree interpolationTree = IntervalTree(start, end);
    Interval::Score score = interpolationTree.root()->interval().multiDomainScorer(fs, fcalls);

    IntervalNode::SharedPtr currentNode = interpolationTree.root(); 

    if (score == Interval::interesting || score == Interval::boring) {
        verificationStack.push(interpolationTree.root());
    }
    else if (score == Interval::ambiguous) {
        ambiguousStack.push(interpolationTree.root());
    }

    while (!verificationStack.empty() || !ambiguousStack.empty()) {
        if (!verificationStack.empty()) {
            currentNode = verificationStack.top();
            verificationStack.pop();
            
            score = currentNode->interval().score();

            x1 = currentNode->interval().start();
            x2 = currentNode->interval().end();

            // Check whether current node has been visited
            if (std::find(visited.begin(), visited.end(),  currentNode) == visited.end() && x2 - x1 >= dfsTerminationDistance) {
                visited.push_back(currentNode);

                if (score == Interval::interesting) {
                    IntervalNode::SharedVec subtree = interpolationTree.intervalLevels(currentNode, verificationLevels);
                    IntervalNode::SharedVec leaves;
                    IntervalNode::SharedVec interestingLeaves;
                    IntervalNode::SharedVec boringLeaves;

                    for (IntervalNode::SharedPtr node : subtree) {
                        score = node->interval().multiDomainScorer(fs, fcalls);

                        if (node->isLeaf()) {
                            leaves.push_back(node);
                            if (score == Interval::interesting) {
                                interestingLeaves.push_back(node);
                            }
                            else if ((score = Interval::boring)) {
                                boringLeaves.push_back(node);
                            }
                        }
                    }

                    if (interestingLeaves.size() == 1 && boringLeaves.size() == leaves.size() - 1) {
                        IntervalNode::SharedPtr node = interestingLeaves.front();
                        x1 = node->interval().start();
                        x2 = node->interval().end();

                        ddfunc f = [&](double x) { 
                            return fs[node->interval().winner().first](x) - fs[node->interval().winner().second](x); 
                        };
                        ddfunc fprime = [&](double x) { 
                            return fprimes[node->interval().winner().first](x) - fprimes[node->interval().winner().second](x); 
                        };

                        double switchPoint = linearRootFinder(x1, x2, f, fprime, fcalls, verificationLevels, subdivisionLevels, dfsTerminationDistance, maxNewtonIterations);

                        // Check if valid value was return (invalid = start - 1)
                        if (switchPoint != node->interval().start() - 1) {
                            int domain = node->interval().winner().first;
                            switchPointStore.push_back({switchPoint, domain});
                        }
                    }
                    else {
                        ambiguousStack.push(currentNode);
                    }
                }
                
                else if (score == Interval::boring) {
                    IntervalNode::SharedVec subtree = interpolationTree.intervalLevels(currentNode, verificationLevels);
                    
                    for (IntervalNode::SharedPtr node : subtree) {
                        score = node->interval().multiDomainScorer(fs, fcalls);

                        if (score != Interval::boring) {
                            ambiguousStack.push(currentNode);
                            visited.pop_back();
                            break;
                        }
                    }
                }
            }
        }
        else if (!ambiguousStack.empty()) {
            currentNode = ambiguousStack.top();
            ambiguousStack.pop();

            x1 = currentNode->interval().start();
            x2 = currentNode->interval().end();

            // Check whether current node has been visited
            if (std::find(visited.begin(), visited.end(),  currentNode) == visited.end() && x2 - x1 >= dfsTerminationDistance) {
                visited.push_back(currentNode);

                IntervalNode::SharedVec subtree = interpolationTree.intervalLevels(currentNode, subdivisionLevels);

                for (IntervalNode::SharedPtr node : subtree) {
                    score = node->interval().multiDomainScorer(fs, fcalls);
                    
                    if (score == Interval::interesting || score == Interval::boring) {
                            verificationStack.push(node);
                    }
                    else if (score == Interval::ambiguous) {
                            ambiguousStack.push(node);
                    }
                }
            }
        }
    }

    std::sort(switchPointStore.begin(), switchPointStore.end());

    std::vector<int> domains;
    std::vector<double> switchPoints = { start };

    for (std::pair<double, int> p : switchPointStore)
    {
        switchPoints.push_back(p.first);
        domains.push_back(p.second);
    }

    return {domains, switchPoints};
}

#endif
