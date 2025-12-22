import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.stats import sem, t
import timeit
from matplotlib.patches import Patch

# ================================
# Operator Definitions
# ================================

limited_type = lambda x:x

def smoothmax(a, b, k=10):
    # a = np.asarray(a)
    # b = np.asarray(b)
    delta = b - a
    exp_term = limited_type(np.exp(-k * delta))

    denom1 = 1.0 + exp_term
    denom2 = 1.0 - exp_term

    term1 = 0.5 * (delta / denom1)
    term2 = 0.5 * (delta / np.where(np.abs(denom2) < 0.0001, 1.0, denom2)) # Avoid division by 0, just for the warnings

    smooth = limited_type(a + term1 + term2)

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
    softplus = limited_type(np.log1p(np.exp(k * delta))) / k
    return limited_type(a + softplus)

def logsumexp_max(a, b, k=10):
    max_ab = max(a, b)
    return limited_type(limited_type(np.log(limited_type(limited_type(np.exp(k * (a - max_ab))) + limited_type(np.exp(k * (b - max_ab))))) + k * max_ab) / k)

def smooth_absmax(a, b, k=10):
    delta = abs(a - b)
    return limited_type((a + b) / 2 + limited_type(np.log1p(np.exp(-k * delta))) / (2 * k))

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
    for name, (func, uses_k, _, _) in operators.items():
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
    for name, (func, uses_k, _, _) in operators.items():
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
        'LogSumExp': (logsumexp_max, True, 45, "orange"),
        'Smooth AbsMax': (smooth_absmax, True, 35, "limegreen"),
        'Smoothmax (Ours)': (smoothmax, True, 30, "royalblue"),
        # 'Quadratic Softmax': (quadratic_softmax, True, 15, "red"),
        # 'Average': (average_max, False, 1)
    }


    # Benchmark configuration
    max_n = 500
    step = 10
    num_trials = 100
    k_value = 50
    confidence = 0.95
    n_values = range(2, max_n + 1, step)
    highest_value = 100
    np.random.seed(42)
    base_data = np.random.uniform(0, highest_value, size=max_n * 2)

    limited_types = [np.float64, np.float32, np.float16, lambda x:x]
    limited_type = limited_types[2]
    #
    # # ==========================================================
    # # Run Add-ons
    # # ==========================================================
    #
    # # Run these after operator_registry is defined
    # runtime_benchmark(operator_registry, num_runs=10000, k=k_value)
    # aggregation_error_distribution(operator_registry, k=k_value, num_trials=num_trials, vector_size=50)

    # ================================
    # Results: Chaining Depth & Random Input Distribution
    # ================================
    #
    # chaining_results = {name: {'mean': [], 'ci_low': [], 'ci_high': []} for name in operator_registry}
    #
    # for n in n_values:
    #     for name, (op_func, uses_k, _, _) in operator_registry.items():
    #         errors = []
    #         for _ in range(num_trials):
    #             x_list = random.sample(list(base_data), n)
    #             true_max = max(x_list)
    #             approx = chain_operator(x_list, op_func, k=k_value if uses_k else None)
    #             errors.append(abs(approx - true_max))
    #         mean = np.mean(errors)
    #         std_err = sem(errors)
    #         h = std_err * t.ppf((1 + confidence) / 2., num_trials - 1)
    #         chaining_results[name]['mean'].append(mean)
    #         chaining_results[name]['ci_low'].append(mean - h)
    #         chaining_results[name]['ci_high'].append(mean + h)
    #
    #
    # # ================================
    # # Plotting: Chaining Depth
    # # ================================
    #
    # plt.figure(figsize=(12, 7))
    # for name in operator_registry:
    #     mean = np.array(chaining_results[name]['mean'])
    #     ci_low = np.array(chaining_results[name]['ci_low'])
    #     ci_high = np.array(chaining_results[name]['ci_high'])
    #     plt.plot(n_values, mean, label=f'{name} Mean Error', color=operator_registry[name][3])
    #     plt.fill_between(n_values, ci_low, ci_high, alpha=0.2)
    #
    # plt.xlabel('Number of Inputs Chained')
    # plt.ylabel(f"Absolute Error (Mean ± {int(confidence * 100)}% CI)")
    # plt.yscale('log')
    # plt.title(f'Chaining Error vs. Input Count (k={k_value}, {num_trials} trials)')
    # plt.legend()
    # plt.grid(True)
    # plt.tight_layout()
    # plt.show()

    # =======================================
    # Plotting: Precision vs. Sharpness (Dual Scale)
    # =======================================
    #
    # diff_values = [0.0001, 0.01, 1.0, 100.0]
    #
    # for i, diff in enumerate(diff_values):
    #     a = 0.01
    #     b = a + diff
    #     true_max = max(a, b)
    #
    #     # Define k-ranges
    #     k_values_low = np.logspace(0.5, 2, 500)
    #     k_values_high = np.arange(1000, 100000, 500)
    #
    #     # Container for both low and high results
    #     k_ranges = [("Low k", k_values_low), ("High k", k_values_high)]
    #
    #     fig, axs = plt.subplots(1, 2, figsize=(16, 6))
    #     fig.suptitle(f'Precision vs. Sharpness (a = {a}, b = {b}, diff = {diff})')
    #
    #     for ax, (k_label, k_values) in zip(axs, k_ranges):
    #         precision_results = {name: [] for name in operator_registry}
    #
    #         for k in k_values:
    #             for name, (op_func, uses_k, _, _) in operator_registry.items():
    #                 approx = op_func(a, b, k)
    #                 precision_results[name].append(abs(approx - true_max))
    #
    #         for name, errors in precision_results.items():
    #             color = operator_registry[name][3]
    #             ax.plot(k_values, errors, label=name, color=color)
    #
    #         ax.set_title(f"{k_label} values")
    #         ax.set_xlabel('Sharpness (k)')
    #         # ax.set_xscale('log' if k_label == "Low k" else 'linear')
    #         ax.set_yscale('symlog', linthresh=1e-8)
    #         ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    #
    #     axs[0].set_ylabel('Approximation Error')
    #     axs[1].legend(title="Operator", loc='upper right')
    #     plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    #     plt.show()

    # ============================================
    # Plot: Error Distribution by Precision
    # ============================================

    # Prepare error data per (type, operator)
    precision_error_data = {}  # {(type_label, operator_name): [errors]}

    type_labels = ['float64', 'float32', 'float16']
    vector_size = 50
    num_trials = 200

    # Plotting
    # plt.figure(figsize=(14, 7))
    fig, axes = plt.subplots(1, 3, sharey=True)

    for iPlot, k_value in enumerate([5, 10, 50]):

        for t_label, ltype in zip(type_labels, limited_types):
            limited_type = ltype  # override global limited_type for this block
            for name, (op_func, uses_k, _, _) in operator_registry.items():
                errors = []
                for _ in range(num_trials):
                    x_list = np.random.uniform(0, highest_value, size=vector_size)
                    true_max = np.max(x_list)
                    k = k_value # if name == "Smoothmax (Ours)" else k_value * 2
                    approx = chain_operator(list(x_list), op_func, k=k if uses_k else None)
                    error = abs(float(approx) - float(true_max))  # force to Python float for plotting
                    errors.append(error)
                precision_error_data[(t_label, name)] = errors

        # Reformat for plotting: group by float type
        boxplot_data = []
        boxplot_labels = []
        group_labels = []

        for t_label in type_labels:
            for name in operator_registry:
                boxplot_data.append(precision_error_data[(t_label, name)])
                boxplot_labels.append(name)
                group_labels.append(t_label)

        positions = []
        width = 0.2
        ticks = []
        tick_labels = []
        colors = []

        op_names = list(operator_registry.keys())
        num_ops = len(op_names)

        # Offset each boxplot group slightly so they don’t overlap
        for i, t_label in enumerate(type_labels):
            for j, op_name in enumerate(op_names):
                index = i * num_ops + j
                pos = i + (j - (num_ops - 1) / 2) * width
                positions.append(pos)
                ticks.append(i)
                tick_labels.append(t_label)
                colors.append(operator_registry[op_name][3])

        bplot = axes[iPlot].boxplot(boxplot_data, positions=positions, widths=width, patch_artist=True, showfliers=False)

        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)

        # Add labels manually
        axes[iPlot].set_xticks(ticks=list(range(len(type_labels))), labels=type_labels)
        axes[iPlot].set_yscale('log')

        legend_patches = [Patch(label=op_name, facecolor=operator_registry[op_name][3]) for op_name in op_names]
        plt.legend(legend_patches, op_names, title='SmoothMax Operator', loc='upper left')
        axes[iPlot].grid(True, which='both', linestyle='--', linewidth=0.5)

        if iPlot == 0:
            axes[iPlot].set_ylabel('Absolute Error')
        if iPlot == 1:
            axes[iPlot].set_title(f'Error distribution by floating-point precision (chaining {vector_size} times the smooth max operators)')
            axes[iPlot].set_xlabel('Floating Point Precision Type')
        axes[iPlot].set_xlabel(f"k = {k_value}")

    plt.tight_layout()
    plt.show()

