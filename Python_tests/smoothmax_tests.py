import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.stats import sem, t
import timeit
import numpy as np

# ================================
# Operator Definitions
# ================================


import numpy as np


def smoothmax(a, b, k=10):
    a = np.asarray(a)
    b = np.asarray(b)
    delta = b - a
    exp_term = np.exp(-k * delta)

    denom1 = 1.0 + exp_term
    denom2 = 1.0 - exp_term + 1e-12  # epsilon to avoid divide-by-zero

    term1 = 0.5 * (delta / denom1)
    term2 = 0.5 * (delta / denom2)

    smooth = a + term1 + term2

    # Use np.where to handle a == b case element-wise or for scalars
    return np.where(a == b, a + 1.0 / (2 * k), smooth)


#def smoothmax(a, b, k=10):
#     delta = b - a
#     if a == b:
#         return a + 1 / (2 * k)
#     exp_term = np.exp(-k * delta)
#     denom1 = 1.0 + exp_term
#     denom2 = 1.0 - exp_term
#     term1 = 0.5 * (delta / denom1)
#     term2 = 0.5 * (delta / denom2)
#     return a + term1 + term2

def softplus_max(a, b, k=10):
    delta = b - a
    softplus = np.log1p(np.exp(k * delta)) / k
    return a + softplus

def logsumexp_max(a, b, k=10):
    max_ab = max(a, b)
    return (np.log(np.exp(k * (a - max_ab)) + np.exp(k * (b - max_ab))) + k * max_ab) / k

def smooth_absmax(a, b, k=10):
    delta = abs(a - b)
    return (a + b) / 2 + np.log1p(np.exp(-k * delta)) / (2 * k)

def quadratic_softmax(a, b, k=10):
    epsilon = 1.0 / k  # sharpness control
    delta = b - a
    if True or abs(delta) < epsilon:
        r = (a + b) / 2 + (delta**2) / (8 * epsilon)
        return r
    else:
        return max(a, b)



def average_max(a, b, k=None):
    return (a + b) / 2.0

# ================================
# Chaining Utility
# ================================

def chain_operator(x_list, operator, k=10):
    result = x_list[-1]
    for x in reversed(x_list[:-1]):
        result = operator(x, result, k) if k is not None else operator(x, result)
    return result


def runtime_benchmark(operators, num_runs=10000, k=10):
    a, b = 1.23, 4.56
    print("\n=== Runtime Benchmark (Avg µs per call over {} runs) ===".format(num_runs))
    for name, (func, uses_k, _) in operators.items():
        if uses_k:
            stmt = f"func({a}, {b}, {k})"
        else:
            stmt = f"func({a}, {b})"
        setup = f"from __main__ import func"
        # Use a lambda to time dynamically bound function
        timer = timeit.Timer(lambda: func(a, b, k) if uses_k else func(a, b))
        total_time = timer.timeit(number=num_runs)
        avg_time_us = (total_time / num_runs) * 1e6
        print(f"{name:<20}: {avg_time_us:.2f} µs per call")

# ==========================================================
# Multivariate Aggregation Behavior (Error Distribution)
# ==========================================================

def aggregation_error_distribution(operators, k=10, num_trials=1000, vector_size=50):
    np.random.seed(42)
    print(f"\n=== Aggregation Error Distribution (n = {vector_size}, {num_trials} trials) ===")
    for name, (func, uses_k, _) in operators.items():
        errors = []
        for _ in range(num_trials):
            x_list = np.random.uniform(0, 100, size=vector_size)
            true_max = np.max(x_list)
            approx = chain_operator(list(x_list), func, k=k if uses_k else None)
            errors.append(abs(approx - true_max))
        mean = np.mean(errors)
        std = np.std(errors)
        p95 = np.percentile(errors, 95)
        print(f"{name:<20}: Mean = {mean:.3e}, Std = {std:.3e}, 95th% = {p95:.3e}")




if __name__ == "__main__":
    # ================================
    # Benchmark Parameters
    # ================================

    # operator_registry = {
    #     'Smoothmax': (smoothmax, True, 30),
    #     'Softplus': (softplus_max, True, 42),
    #     'LogSumExp': (logsumexp_max, True, 45),
    #     #'Average': (average_max, False, 1)
    # }
    operator_registry = {
        'Smoothmax (Ours)': (smoothmax, True, 30),
        'LogSumExp': (logsumexp_max, True, 45),
        'Smooth AbsMax': (smooth_absmax, True, 35),
        # 'Quadratic Softmax': (quadratic_softmax, True, 15),
        # 'Average': (average_max, False, 1)
    }


    # Benchmark configuration
    max_n = 100
    step = 10
    num_trials = 100
    k_value = 5
    confidence = 0.95
    n_values = range(2, max_n + 1, step)
    np.random.seed(42)
    base_data = np.random.uniform(0, 100, size=max_n * 2)


    # ==========================================================
    # Run Add-ons
    # ==========================================================

    # Run these after operator_registry is defined
    runtime_benchmark(operator_registry, num_runs=10000, k=k_value)
    aggregation_error_distribution(operator_registry, k=k_value, num_trials=num_trials, vector_size=50)

    # ================================
    # Results: Chaining Depth & Random Input Distribution
    # ================================

    chaining_results = {name: {'mean': [], 'ci_low': [], 'ci_high': []} for name in operator_registry}

    for n in n_values:
        for name, (op_func, uses_k, _) in operator_registry.items():
            errors = []
            for _ in range(num_trials):
                x_list = random.sample(list(base_data), n)
                true_max = max(x_list)
                approx = chain_operator(x_list, op_func, k=k_value if uses_k else None)
                errors.append(abs(approx - true_max))
            mean = np.mean(errors)
            std_err = sem(errors)
            h = std_err * t.ppf((1 + confidence) / 2., num_trials - 1)
            chaining_results[name]['mean'].append(mean)
            chaining_results[name]['ci_low'].append(mean - h)
            chaining_results[name]['ci_high'].append(mean + h)


    # ================================
    # Results: Precision and Sharpness Sweep
    # ================================

    a = 1.000
    b = 1.0001
    true_max = max(a, b)
    k_values_extended = np.logspace(0.1, 2, 500)

    precision_results = {name: [] for name in operator_registry}

    for k in k_values_extended:
        for name, (op_func, uses_k, _) in operator_registry.items():
            args = (a, b, k)
            approx = op_func(*args)
            precision_results[name].append(np.abs(approx - true_max))

    # ================================
    # Plotting: Chaining Depth
    # ================================

    plt.figure(figsize=(12, 7))
    for name in operator_registry:
        mean = np.array(chaining_results[name]['mean'])
        ci_low = np.array(chaining_results[name]['ci_low'])
        ci_high = np.array(chaining_results[name]['ci_high'])
        plt.plot(n_values, mean, label=f'{name} Mean Error')
        plt.fill_between(n_values, ci_low, ci_high, alpha=0.2)

    plt.xlabel('Number of Inputs Chained')
    plt.ylabel(f"Absolute Error (Mean ± {int(confidence * 100)}% CI)")
    plt.yscale('log')
    plt.title(f'Chaining Error vs. Input Count (k={k_value}, {num_trials} trials)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # ================================
    # Plotting: Precision vs. Sharpness
    # ================================

    plt.figure(figsize=(12, 7))
    for name, errors in precision_results.items():
        cost = operator_registry[name][2]
        plt.plot(k_values_extended, np.array(errors), label=f'{name}')

    # plt.axhline(0, color='black', linestyle=':')
    plt.xlabel('Sharpness (k)')
    plt.ylabel('Scaled Approximation Error')
    # plt.xscale('log')
    plt.yscale('symlog', linthresh=1e-8)
    plt.title('Precision vs. Sharpness')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

