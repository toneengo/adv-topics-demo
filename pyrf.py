from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from scipy.special import binom
import random
import math
from uncertainties.umath import *
from matplotlib.colors import ListedColormap
from treelib import Tree

#matplotlib widget
DFS_LIMIT = 35
LEVEL_2D = 30
SEED = 3
DEGREE = 10

x_points = np.linspace(0, 1, 200)


def bernstein_generator(n, i):
    return binom(n, i) * (x_points**i) * (1 - x_points)**(n - i)


def bernstein_eqn(n, i, x):
    return binom(n, i) * (x**i) * (1 - x)**(n - i)


def bernstein_graph(n):
    color = iter(cm.rainbow(np.linspace(0, 1, n + 1)))
    for j in range(n + 1):
        c = next(color)
        y = bernstein_generator(n, j)
        axs[0].plot(x_points, y, c=c)


def bernstein_graph_sum(n, m, c, name):
    y = []
    for j in m:
        if len(y):
            y = np.add(y, bernstein_generator(n, j))
        else:
            y = bernstein_generator(n, j)
    axs[0].plot(x_points, y, c=c, label=name)


def bernstein_eqn_sum_value(n, m, x):
    y = 0
    for j in m:
        y += bernstein_eqn(n, j, x)
    return y
# ###################################################################


def Bernstein(n, i):
    def fxn(x):
        return binom(n, i) * (x**i) * (1 - x)**(n - i)
    return fxn


def Bernstein_prime(n, i):
    def fxn(x):
        try:
            return binom(n, i) * (i * (x**(i - 1)) * ((1 - x) **
                                                      (n - i)) - (x**i) * (n - i) * (1 - x)**(n - i - 1))
        except ZeroDivisionError:
            return 0
    return fxn


n = DEGREE
for i in range(n + 1):
    globals()['fxn' + '_{}'.format(i)] = Bernstein(n, i)
    globals()['fxn_prime' + '_{}'.format(i)] = Bernstein_prime(n, i)


def eqn_creator(to_sum):
    def fxn_sum(x):
        y = 0
        for i in to_sum:
            y += globals()['fxn' + '_{}'.format(i)](x)
        return y

    def fxn_prime_sum(x):
        y_prime = 0
        for i in to_sum:
            y_prime += globals()['fxn_prime' + '_{}'.format(i)](x)
        return y_prime
    return fxn_sum, fxn_prime_sum


a = np.arange(0, DEGREE + 1)
random.Random(SEED).shuffle(a)
newarr = np.array_split(a, 4)

eqn_1, eqn_1_prime = eqn_creator(newarr[0])
eqn_2, eqn_2_prime = eqn_creator(newarr[1])
eqn_3, eqn_3_prime = eqn_creator(newarr[2])
eqn_4, eqn_4_prime = eqn_creator(newarr[3])

def makelist():
    fs = [eqn_1, eqn_2, eqn_3, eqn_4]
    fprimes = [eqn_1_prime, eqn_2_prime, eqn_3_prime, eqn_4_prime]
    return fs, fprimes
# *****************************************************************************


def eqn_1_2(x):
    return eqn_1(x) - eqn_2(x)


def eqn_1_3(x):
    return eqn_1(x) - eqn_3(x)


def eqn_1_4(x):
    return eqn_1(x) - eqn_4(x)


def eqn_2_3(x):
    return eqn_2(x) - eqn_3(x)


def eqn_2_4(x):
    return eqn_2(x) - eqn_4(x)


def eqn_3_4(x):
    return eqn_3(x) - eqn_4(x)


def eqn_1_2_prime(x):
    return eqn_1_prime(x) - eqn_2_prime(x)


def eqn_1_3_prime(x):
    return eqn_1_prime(x) - eqn_3_prime(x)


def eqn_1_4_prime(x):
    return eqn_1_prime(x) - eqn_4_prime(x)


def eqn_2_3_prime(x):
    return eqn_2_prime(x) - eqn_3_prime(x)


def eqn_2_4_prime(x):
    return eqn_2_prime(x) - eqn_4_prime(x)


def eqn_3_4_prime(x):
    return eqn_3_prime(x) - eqn_4_prime(x)

# *****************************************************************************
# ##################################################################


guesses_collector = []


def newton(
    func,
    x0,
    fprime=None,
    args=(),
    tol=1e-8,
    maxiter=100,
    rtol=0.0,
    disp=True,
):

    if tol <= 0:
        raise ValueError("tol too small (%g <= 0)" % tol)
    if maxiter < 1:
        raise ValueError("maxiter must be greater than 0")

    # Convert to float (don't use float(x0); this works also for complex x0)
    p0 = 1.0 * x0
    funcalls = 0
    fderp = None
    # Newton-Raphson method
    for itr in range(maxiter):
        # first evaluate fval
        fval = func(x=p0, *args)

        funcalls += 1
        if fderp is None:
            fder = fprime(x=p0, *args)
            funcalls += 1
        else:
            fder = fderp
        # If fval is 0, a root has been found, then terminate
        if fval == 0:
            return p0, funcalls, fder
        if fder == 0:
            msg = "Derivative was zero."
            if disp:
                msg += " Failed to converge after %d iterations, value is %s." % (
                    itr + 1, p0, )
                pass
            return None, funcalls, None
        newton_step = fval / fder
        p = p0 - newton_step
        fder = fprime(x=p, *args)
        funcalls += 1
        #if np.isclose(p, p0, rtol=rtol, atol=tol):
        if abs(p-p0) < tol:
            return p, funcalls, fderp
        guesses_collector.append((p0, p))
        p0 = p
        
    if disp:
        msg = "Failed to converge after %d iterations, value is %s." % (
            itr + 1, p)
        return None, funcalls, None

    return p, funcalls, fderp


class Interval(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.score = None
        self.visited = None

        # Some extra variables used by multi_domain
        self.winner_pos = None
        self.transition = None

    def _linear_interpolation_scorer(self, f, fprime):
        start = self.start
        end = self.end
        ans = []
        m_store = []
        b_store = []
        sign_check = 1
        points = [start, end]
        for point in points:
            current_point = point
            z1 = f(x=current_point)
            der = fprime(x=current_point)
            sign_check *= z1
            if der == 0:
                m_store.append(1e-6)
            else:
                m_store.append(der)
            b_store.append(-current_point * der + z1)
        ans.extend([-1 * (b_store[0] / m_store[0]), -
                    1 * (b_store[1] / m_store[1])])
        diff = end - start

        if start <= ans[0] <= end and start <= ans[1] <= end:
            self.score = 2
            return 2
        if (start <= ans[0] <= end or start <= ans[1] <= end) or (
                start - diff < ans[1] < start or end < ans[0] < end + diff):
            self.score = 1
            return 1
        else:
            if sign_check < 0:
                # There is a sign change
                self.score = 1
                return 1
            else:
                self.score = -1
                return -1

    # This function is customised to 4 functions.
    def _multi_domain_scorer(self, eqn_1, eqn_2, eqn_3, eqn_4):
        # returns token, winner/None, global token

        start = self.start
        end = self.end

        # Top 2 of each end of the interval.
        start_winner = np.argpartition(
            np.array(
                [
                    eqn_1(
                        x=start), eqn_2(
                        x=start), eqn_3(
                        x=start), eqn_4(
                            x=start)]), -2)[
                                :-3:-1]
        end_winner = np.argpartition(np.array([eqn_1(x=end), eqn_2(
            x=end), eqn_3(x=end), eqn_4(x=end)]), -2)[:-3:-1]

        # Explanation of conditions below...
        reverse_end_winner = end_winner[::-1]
        reverse_ans = np.array_equal(start_winner, reverse_end_winner)
        # ans = np.array_equal(start_winner, end_winner)
        top_answer_same = False
        if start_winner[0] == end_winner[0]:
            top_answer_same = True

        if top_answer_same:
            self.score = -1
            return -1
        elif reverse_ans:
            Transition = [start_winner[0], end_winner[0]]
            # Sorting is done because...
            fxn_to_compare = np.sort(start_winner)
            winner_pos = []
            winner_pos.append(
                "eqn_" + str(fxn_to_compare[0] + 1) + "_" +
                str(fxn_to_compare[1] + 1))

            winner_pos.append(
                "eqn_" + str(fxn_to_compare[0] + 1) + "_" +
                str(fxn_to_compare[1] + 1) + "_" + "prime")

            self.score = 1
            self.winner_pos = winner_pos
            self.transition = Transition
            return 1
        else:
            self.score = 0
            return 0

    def _interval_levels(self, level, lim, tree):
        x1 = self.start
        x2 = self.end

        queue = []

        def internal_bfs(x1, x2, level, lim):
            sub_tree = tree.subtree(f"({self.start},{self.end})")
            if 0.5 * (x2 - x1) >= lim and (sub_tree.size()
                                           < ((2 ** (level + 2)) - 2) / 2):

                # Check if node has children or not before subdivision.
                if tree.get_node(f"({x1},{x2})").is_leaf():
                    # Creating children
                    tree.create_node(
                        parent=f"({x1},{x2})",
                        identifier=f"({x1},{0.5 * (x1 + x2)})",
                        data=Interval(
                            x1, 0.5 * (x1 + x2)))
                    tree.create_node(
                        parent=f"({x1},{x2})",
                        identifier=f"({0.5 * (x1 + x2)},{x2})",
                        data=Interval(
                            0.5 * (x1 + x2), x2))

                queue.append((x1, 0.5 * (x1 + x2)))
                queue.append((0.5 * (x1 + x2), x2))
                s, e = queue[0][0], queue[0][1]
                queue.pop(0)
                internal_bfs(s, e, level, lim)
        internal_bfs(x1, x2, level, lim)


def linear_root_finder(
    start,
    end,
    f,
    fprime,
    levels_1=2,
    levels_2=1,
    dfs_termination_distance=1/(2**DFS_LIMIT),
    max_newton_iter=250,

):
    root_store = []
    x1 = start
    x2 = end

    # Q: Currently, there are only 2 stacks, verification_stack handles 2 and -1
    # cases both, however is there any benefit to create a seperate stack which
    # handles -1 cases?
    verification_stack = []
    stack = []

    visited = {}
    total_calls = 0
    span_collector = []
    newton_calls_collector = []
    newton_found_root_outside_boundary = []

    interpolation_tree = Tree()
    root_interval = interpolation_tree.create_node(
        identifier=f"({x1},{x2})", data=Interval(
            x1, x2))

    token = root_interval.data._linear_interpolation_scorer(f, fprime)
    total_calls += 4
    if token == 2 or token == -1:
        verification_stack.append((f"({x1},{x2})", token))
    elif token == 1:
        stack.append((f"({x1},{x2})", token))

    while len(stack) or len(verification_stack):
        if verification_stack:
            interval_id, token = verification_stack[-1][0], verification_stack[-1][1]
            verification_stack.pop()
            current_interval = interpolation_tree.get_node(interval_id)
            x1 = current_interval.data.start
            x2 = current_interval.data.end

            if current_interval.data.visited is None and (
                    x2 - x1 >= dfs_termination_distance):
                current_interval.data.visited = True

                if current_interval.data.score == 2:

                    # Directly apply newton
                    cur_root, newton_calls, gradient = newton(
                        f,
                        current_interval.data.start,
                        fprime,
                        maxiter=max_newton_iter,
                        tol=1e-8,
                    )
                    total_calls += newton_calls
                    newton_calls_collector.append(newton_calls)

                    if cur_root is not None:
                        if cur_root >= start and cur_root <= end:
                            root_store.append(cur_root)
                        return root_store, total_calls
                else:
                    # Verification for boring
                    current_interval.data._interval_levels(
                        levels_2, dfs_termination_distance, interpolation_tree)
                    # We only go a level deep to verify.
                    children = interpolation_tree.children(interval_id)
                    for node_itr in children:
                        # Check if score is already alloted.
                        if node_itr.data.score is None:
                            internal_token = node_itr.data._linear_interpolation_scorer(
                                f, fprime)
                            total_calls += 4
                        if node_itr.data.score != -1:
                            if node_itr.data.score == 2:
                                verification_stack.append(
                                    (node_itr.identifier, node_itr.data.score))
                            else:
                                stack.append(
                                    (node_itr.identifier, node_itr.data.score))

        elif stack:
            interval_id, token = stack[-1][0], stack[-1][1]
            stack.pop()
            current_interval = interpolation_tree.get_node(interval_id)
            x1 = current_interval.data.start
            x2 = current_interval.data.end

            if current_interval.data.visited is None and (
                    x2 - x1 >= dfs_termination_distance):
                current_interval.data.visited = True

                current_interval.data._interval_levels(
                    levels_2, dfs_termination_distance, interpolation_tree)

                # We only go a level deep.
                internal_tree = interpolation_tree.subtree(interval_id)
                for node_itr in internal_tree.all_nodes_itr():
                    # Check if score is already alloted.
                    if node_itr.data.score is None:
                        internal_token = node_itr.data._linear_interpolation_scorer(
                            f, fprime)
                        total_calls += 4
                    if node_itr.data.score == 2 or node_itr.data.score == -1:
                        verification_stack.append(
                            (node_itr.identifier, node_itr.data.score))
                    elif node_itr.data.score == 1:
                        stack.append((node_itr.identifier, node_itr.data.score))

    return root_store, total_calls


def highest_value_eqn(X):
    return np.argpartition(np.array([eqn_1(x=X), eqn_2(
        x=X), eqn_3(x=X), eqn_4(x=X)]), -2)[-1]


def dfs_multi_domain(
        start, end, a=[], b=[], levels_1=2, levels_2=1,
        dfs_termination_distance=1 / (2 ** DFS_LIMIT)):
    switch_points_store = {}
    x1 = start
    x2 = end
    verification_stack = []
    stack = []
    visited = {}
    visited_2_debug = []
    total_calls = 0
    span_collector = []

    multi_domain_tree = Tree()
    root_interval = multi_domain_tree.create_node(
        identifier=f"({x1},{x2})", data=Interval(
            x1, x2))

    token = root_interval.data._multi_domain_scorer(eqn_1, eqn_2, eqn_3, eqn_4)
    total_calls += 6

    if token == 1 or token == -1:
        verification_stack.append(
            (f"({x1},{x2})", token))
    else:
        stack.append((f"({x1},{x2})", token))

    while len(stack) or len(verification_stack):
        if verification_stack:
            interval_id, token = verification_stack[-1][0], verification_stack[-1][1]
            visited_2_debug.append(verification_stack.pop())
            current_interval = multi_domain_tree.get_node(interval_id)
            x1 = current_interval.data.start
            x2 = current_interval.data.end

            if current_interval.data.visited is None and (
                    x2 - x1 >= dfs_termination_distance):

                current_interval.data.visited = True

                if token == 1:
                    current_interval.data._interval_levels(
                        levels_1, dfs_termination_distance, multi_domain_tree)
                    internal_tree = multi_domain_tree.subtree(interval_id)
                    leaf_scores = []
                    leaf_objs = []
                    for node_itr in internal_tree.all_nodes_itr():
                        # Check if score is already alloted.
                        if node_itr.data.score is None:
                            internal_token = node_itr.data._multi_domain_scorer(
                                eqn_1,
                                eqn_2,
                                eqn_3,
                                eqn_4)
                            total_calls += 6
                        if node_itr.is_leaf():
                            leaf_objs.append(node_itr)
                            leaf_scores.append(node_itr.data.score)

                    leaf_scores = np.array(leaf_scores)
                    leaves_with_score_1 = np.where(leaf_scores == 1)[0]
                    leaves_with_score_minus_1 = np.where(leaf_scores == -1)[0]
                    if len(leaves_with_score_1) == 1 and len(
                            leaves_with_score_minus_1) == len(leaf_scores) - 1:
                        leaf_with_score_1 = leaf_objs[leaves_with_score_1[0]]

                        winner_pos = leaf_with_score_1.data.winner_pos
                        leaf_x1, leaf_x2 = leaf_with_score_1.data.start, leaf_with_score_1.data.end

                        switch_point_array, root_finding_calls = linear_root_finder(
                            leaf_x1, leaf_x2, globals()[winner_pos[0]],
                            globals()[winner_pos[1]],
                            levels_1=levels_1, levels_2=levels_2,
                            dfs_termination_distance=dfs_termination_distance)

                        total_calls += root_finding_calls

                        if switch_point_array:
                            span_collector.append([x1, x2])
                            visited[(x1, x2)] = switch_point_array[0]
                            if switch_point_array[0] >= start and switch_point_array[0] <= end:
                                switch_points_store[switch_point_array[0]
                                                    ] = leaf_with_score_1.data.transition

                    else:
                        stack.append((interval_id, token))
                        current_interval.data.visited = None

                # Verifies if all children are -1 (boring) and ends further
                # subdivision.
                elif token == -1:
                    current_interval.data._interval_levels(
                        levels_1, dfs_termination_distance, multi_domain_tree)
                    internal_tree = multi_domain_tree.subtree(interval_id)
                    for node_itr in internal_tree.all_nodes_itr():
                        # Check if score is already alloted.
                        if node_itr.data.score is None:
                            internal_token = node_itr.data._multi_domain_scorer(
                                eqn_1,
                                eqn_2,
                                eqn_3,
                                eqn_4)
                            total_calls += 6
                        if node_itr.data.score != -1:
                            # Goes to the normal stack now
                            stack.append((interval_id, token))
                            current_interval.data.visited = None
                            break

        else:
            interval_id, token = stack[-1][0], stack[-1][1]
            visited_2_debug.append(stack.pop())
            current_interval = multi_domain_tree.get_node(interval_id)
            x1 = current_interval.data.start
            x2 = current_interval.data.end

            if current_interval.data.visited is None and (
                    x2 - x1 >= dfs_termination_distance):
                current_interval.data.visited = True

                current_interval.data._interval_levels(
                    levels_2, dfs_termination_distance, multi_domain_tree)

                # We only go a level deep, always use children() here as
                # subtree() will get everything.
                children = multi_domain_tree.children(interval_id)
                for node_itr in children:
                    # Check if score is already alloted.
                    if node_itr.data.score is None:
                        internal_token = node_itr.data._multi_domain_scorer(
                            eqn_1,
                            eqn_2,
                            eqn_3,
                            eqn_4)
                        total_calls += 6
                    if node_itr.data.score == 1 or node_itr.data.score == -1:
                        verification_stack.append(
                            (node_itr.identifier, node_itr.data.score))
                    elif node_itr.data.score == 0:
                        stack.append((node_itr.identifier, node_itr.data.score))

    switch_points_store = {
        k: switch_points_store[k] for k in sorted(switch_points_store)}
    points_of_transition = list(switch_points_store.keys())
    values_ans = list(switch_points_store.values())

    if len(values_ans):
        domain_in_each_region = [i[0] for i in values_ans]
        domain_in_each_region.append(values_ans[-1][1])

        points_of_transition.insert(0, start)
        points_of_transition.append(end)

        return domain_in_each_region, points_of_transition
    else:
        # If there are no transition points, then that means it's all the same domain.
        # So return the domain of the midpoint of that line and [start, end]
        # as the default transition points.
        domain_in_each_region = highest_value_eqn(
            (start + end) / 2)
        return [domain_in_each_region], [
            start, end]