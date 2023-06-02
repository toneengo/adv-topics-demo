#ifndef LINEAR_ROOT_FINDER_CPP
#define LINEAR_ROOT_FINDER_CPP

#include <stack>
#include <algorithm>

#include "IntervalTree.cpp"
#include "Comp.cpp"

const int DFS_TERMINATION_LIMIT = 35;

// @return Tuple of { root, derivative }
std::tuple<double, double> newton(
    ddfunc const &f, 
    ddfunc const &fprime,
    double guess,
    int &fcalls,
    double tol = 1.48e-8,
    int maxiter = 100
) {
    
    if (tol <= 0) {
        throw std::domain_error("Tolerance (tol) is too low (" + std::to_string(tol) + " <= 0)");
    }
    if (maxiter < 1) {
        throw std::domain_error("Max iterations (maxiter) is too low (" + std::to_string(maxiter) + " < 1)");
    }

    double x1 = guess;
    double fval = -1;
    double fder = -1;

    for (int iter = 0; iter < maxiter; iter++) {
        fval = f(x1);
        fder = fprime(x1);

        fcalls += 2;

        // If zero is found
        if (fval == 0) {
            return {x1, fder};
        }

        // If derivative is zero
        if(fder == 0) {
            std::cout << "Derivative was zero\n";
            std::cout << "Failed to converge after " << iter << " iterations, value is " << fval << '\n';
            return {-1, 0};
        }

        double newtonStep = fval / fder;
        double x2 = x1 - newtonStep;

        fder = fprime(x1); // Unnecessary call?
        fcalls += 1;

        if (Comp::dcomp(x1, x2, tol)) {
            return {x2, fder};
        }

        x1 = x2;
    }

    std::cout << "Failed to converge after " + std::to_string(maxiter + 1) + " iterations, value is " + std::to_string(x1) << '\n';

    return {-1, 0};
}

double linearRootFinder(
    double start,
    double end,
    ddfunc const &f,
    ddfunc const &fprime,
    int &totalCalls,
    int verificationLevels = 3,
    int subdivisionLevels = 2,
    double dfsTerminationDistance = 1 / pow(2, DFS_TERMINATION_LIMIT), // Same a reference but without setting constant
    int maxNewtonIterations = 250
) {
    double x1;
    double x2;

    if (f(start) == 0) {
        return start;
    }

    if (f(end) == 0) {
        return end;
    }

    std::stack<IntervalNode::SharedPtr> verificationStack;
    std::stack<IntervalNode::SharedPtr> ambiguousStack;

    IntervalNode::SharedVec visited;

    IntervalTree interpolationTree = IntervalTree(start, end);
    Interval::Score score = interpolationTree.root()->interval().linearInterpolationScorer(f, fprime, totalCalls);

    if (score == Interval::interesting || score == Interval::boring) {
        verificationStack.push(interpolationTree.root());
    }
    else if (score == Interval::ambiguous) {
        ambiguousStack.push(interpolationTree.root());
    }

    while (!verificationStack.empty() || !ambiguousStack.empty()) {
        if (!verificationStack.empty()) {
            IntervalNode::SharedPtr currentNode = verificationStack.top();
            verificationStack.pop();
            
            score = currentNode->interval().score();

            x1 = currentNode->interval().start();
            x2 = currentNode->interval().end();

            if (
                // Check whether current node has been visited
                std::find(visited.begin(), visited.end(),  currentNode) == visited.end() &&
                x2 - x1 >= dfsTerminationDistance
            ) {
                visited.push_back(currentNode);

                if (score == Interval::interesting) {

                    // { root, derivative, total function calls }
                    double root;
                    double derivative;
                    
                    std::tie(root, derivative) = newton(
                        f,
                        fprime,
                        x1,
                        totalCalls,
                        1e-8,
                        maxNewtonIterations
                    );

                    // Tests if valid value was returned
                    if (!(derivative == 0)) {
                        if (root >= start && root <= end) {
                            return root ;
                        }
                    }
                }
                else if (score == Interval::boring) {
                    IntervalNode::SharedVec children = interpolationTree.intervalLevels(currentNode, subdivisionLevels);

                    for (IntervalNode::SharedPtr child : children) {
                        score = child->interval().linearInterpolationScorer(f, fprime, totalCalls);

                        if (score == Interval::interesting) {
                            verificationStack.push(child);
                        }
                        else if (score == Interval::ambiguous) {
                            ambiguousStack.push(child);
                        }
                        // Pass is score is boring
                    }
                }
            }
        }
        else if (!ambiguousStack.empty()) {
            IntervalNode::SharedPtr currentNode = ambiguousStack.top();
            ambiguousStack.pop();

            x1 = currentNode->interval().start();
            x2 = currentNode->interval().end();

            if (
                // Check whether current node has been visited
                std::find(visited.begin(), visited.end(),  currentNode) == visited.end() &&
                x2 - x1 >= dfsTerminationDistance
            ) {
                visited.push_back(currentNode);

                IntervalNode::SharedVec children = interpolationTree.intervalLevels(currentNode, subdivisionLevels);

                for (IntervalNode::SharedPtr child : children) {
                        score = child->interval().linearInterpolationScorer(f, fprime, totalCalls);
                        
                        if (score == Interval::interesting || score == Interval::boring) {
                                verificationStack.push(child);
                        }
                        else if (score == Interval::ambiguous) {
                                ambiguousStack.push(child);
                        }
                    }
            }
        }
    }

    return (start - 1);
}

#endif
