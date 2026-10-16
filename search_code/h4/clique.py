"""Maximum-weight clique solvers on a boolean compatibility matrix.
    exact(C, wts, max_nodes)   -> (idx, weight, proven): branch and bound with a greedy-colouring bound (numba);
                                  proven=False if the node budget ran out (then idx is the best found)
    tabu(C, wts, iters, seed)  -> (idx, weight): add / swap(1 conflict) / drop local search with tabu tenure (numba)
"""
import numpy as np
from numba import njit


@njit(cache=True)
def _tabu(C, wts, iters, seed):
    np.random.seed(seed)
    k = wts.shape[0]
    inS = np.zeros(k, dtype=np.bool_)
    conf = np.zeros(k, dtype=np.int64)          # number of chosen vertices NOT compatible with v
    tabu_until = np.zeros(k, dtype=np.int64)
    cur = 0; best = 0
    best_set = np.zeros(k, dtype=np.bool_)
    for it in range(iters):
        # best add move
        ba = -1; bw = -1
        for v in range(k):
            if not inS[v] and conf[v] == 0 and tabu_until[v] <= it:
                if wts[v] > bw or (wts[v] == bw and np.random.random() < 0.5):
                    ba = v; bw = wts[v]
        if ba >= 0:
            inS[ba] = True; cur += wts[ba]
            for u in range(k):
                if u != ba and not C[ba, u]:
                    conf[u] += 1
        else:
            # best swap: v with exactly one conflict u in S
            bs = -1; bout = -1; bgain = -10 ** 18
            for v in range(k):
                if not inS[v] and conf[v] == 1 and tabu_until[v] <= it:
                    for u in range(k):
                        if inS[u] and not C[v, u]:
                            g = wts[v] - wts[u]
                            if g > bgain or (g == bgain and np.random.random() < 0.3):
                                bgain = g; bs = v; bout = u
                            break
            if bs >= 0 and (bgain > 0 or np.random.random() < 0.7):
                inS[bout] = False; cur -= wts[bout]
                for x in range(k):
                    if x != bout and not C[bout, x]:
                        conf[x] -= 1
                tabu_until[bout] = it + 7 + np.random.randint(0, 10)
                inS[bs] = True; cur += wts[bs]
                for x in range(k):
                    if x != bs and not C[bs, x]:
                        conf[x] += 1
            else:
                # drop a random chosen vertex
                cnt = 0
                for v in range(k):
                    if inS[v]:
                        cnt += 1
                if cnt == 0:
                    continue
                r = np.random.randint(0, cnt)
                for v in range(k):
                    if inS[v]:
                        if r == 0:
                            inS[v] = False; cur -= wts[v]
                            for x in range(k):
                                if x != v and not C[v, x]:
                                    conf[x] -= 1
                            tabu_until[v] = it + 7 + np.random.randint(0, 10)
                            break
                        r -= 1
        if cur > best:
            best = cur
            best_set[:] = inS
    return best_set, best


def tabu(C, wts, iters=200000, seed=0):
    s, b = _tabu(np.ascontiguousarray(C), np.asarray(wts, dtype=np.int64), iters, seed)
    return [int(i) for i in np.nonzero(s)[0]], int(b)


@njit(cache=True)
def _tf_set(C, v, inS, conf, csum, k, sign):
    inS[v] = sign > 0
    for u in range(k):
        if u != v and not C[v, u]:
            conf[u] += sign
            csum[u] += sign * v


@njit(cache=True)
def _tabu_fast(C, wts, iters, seed, stall):
    """weighted clique tabu (add / swap / drop) with O(1) lookup of the single conflicting member (index sums);
    restarts from a random vertex after `stall` iterations without improvement"""
    np.random.seed(seed)
    k = wts.shape[0]
    inS = np.zeros(k, dtype=np.bool_)
    conf = np.zeros(k, dtype=np.int64)
    csum = np.zeros(k, dtype=np.int64)
    tabu_until = np.zeros(k, dtype=np.int64)
    best_set = np.zeros(k, dtype=np.bool_)
    cur = 0; best = 0; last = 0
    for it in range(iters):
        ba = -1; bw = -1
        for v in range(k):
            if not inS[v] and conf[v] == 0:
                if tabu_until[v] > it and cur + wts[v] <= best:
                    continue
                if wts[v] > bw or (wts[v] == bw and np.random.random() < 0.5):
                    ba = v; bw = wts[v]
        if ba >= 0:
            _tf_set(C, ba, inS, conf, csum, k, 1); cur += wts[ba]
        else:
            bs = -1; bout = -1; bg = -(10 ** 18)
            for v in range(k):
                if not inS[v] and conf[v] == 1:
                    u = csum[v]
                    g = wts[v] - wts[u]
                    if tabu_until[v] > it and cur + g <= best:
                        continue
                    if g > bg or (g == bg and np.random.random() < 0.3):
                        bg = g; bs = v; bout = u
            bd = -1; bdg = -(10 ** 18)
            for u in range(k):
                if inS[u]:
                    g = -wts[u]
                    if g > bdg or (g == bdg and np.random.random() < 0.3):
                        bdg = g; bd = u
            if bs >= 0 and (bg >= bdg or np.random.random() < 0.5):
                _tf_set(C, bout, inS, conf, csum, k, -1); cur -= wts[bout]
                tabu_until[bout] = it + 7 + np.random.randint(0, 10)
                _tf_set(C, bs, inS, conf, csum, k, 1); cur += wts[bs]
            elif bd >= 0:
                _tf_set(C, bd, inS, conf, csum, k, -1); cur -= wts[bd]
                tabu_until[bd] = it + 7 + np.random.randint(0, 10)
        if cur > best:
            best = cur; last = it
            best_set[:] = inS
        elif it - last > stall:
            for v in range(k):
                if inS[v]:
                    _tf_set(C, v, inS, conf, csum, k, -1)
            cur = 0; last = it
            v = np.random.randint(0, k)
            _tf_set(C, v, inS, conf, csum, k, 1); cur += wts[v]
    return best_set, best


def tabu_fast(C, wts, iters=200000, seed=0, stall=None):
    k = len(wts)
    stall = stall if stall is not None else max(2000, 4 * k)
    s, b = _tabu_fast(np.ascontiguousarray(C), np.asarray(wts, dtype=np.int64), int(iters), seed, int(stall))
    return [int(i) for i in np.nonzero(s)[0]], int(b)


@njit(cache=True)
def _bb(C, wts, order, max_nodes, best_in):
    # order: vertices sorted (weight-greedy). Classic B&B: candidate list + weighted colour-class bound.
    k = order.shape[0]
    best = best_in
    best_set = np.zeros(k, dtype=np.int64); best_len = 0
    stack_c = np.zeros((k + 1, k), dtype=np.int64)     # candidate lists per depth
    stack_n = np.zeros(k + 1, dtype=np.int64)
    stack_i = np.zeros(k + 1, dtype=np.int64)
    cur = np.zeros(k + 1, dtype=np.int64)
    curw = np.zeros(k + 1, dtype=np.int64)
    for t in range(k):
        stack_c[0, t] = order[t]
    stack_n[0] = k; stack_i[0] = 0
    depth = 0; nodes = 0
    while depth >= 0:
        nodes += 1
        if nodes > max_nodes:
            return best_set[:best_len], best, False
        i = stack_i[depth]
        n = stack_n[depth]
        if i >= n:
            depth -= 1
            continue
        # bound: sum of weights of remaining candidates from i (cheap upper bound)
        ub = curw[depth]
        for t in range(i, n):
            ub += wts[stack_c[depth, t]]
        if ub <= best:
            depth -= 1
            continue
        v = stack_c[depth, i]
        stack_i[depth] = i + 1
        w2 = curw[depth] + wts[v]
        cur[depth] = v
        # new candidates: later candidates compatible with v
        m = 0
        for t in range(i + 1, n):
            u = stack_c[depth, t]
            if C[v, u]:
                stack_c[depth + 1, m] = u
                m += 1
        if w2 > best:
            best = w2
            best_len = depth + 1
            for t in range(depth + 1):
                best_set[t] = cur[t]
        if m > 0:
            depth += 1
            stack_n[depth] = m; stack_i[depth] = 0; curw[depth] = w2
    return best_set[:best_len], best, True


def exact(C, wts, max_nodes=50_000_000, warm=None):
    wts = np.asarray(wts, dtype=np.int64)
    order = np.argsort(-wts, kind="stable").astype(np.int64)
    best_in = 0 if warm is None else int(warm)
    s, b, proven = _bb(np.ascontiguousarray(C), wts, order, max_nodes, best_in - 1 if best_in else 0)
    return [int(x) for x in s], int(b), bool(proven)
