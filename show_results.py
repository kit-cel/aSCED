import numpy as np 

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
matplotlib.rcParams['text.usetex'] = True



def plot_error_rates(Sim, *args):
    '''
    Plots error rates of simulation and reference data.

    Parameters:
    ---------------------------------------------------------------
    Sim                   - Simulation instance.
    args    (str)         - Key to access data in Sim.error_rates.
            (dict)        - Reference data.
            (tuple)       - Tuple containing reference data & label.
    '''

    types = []          
    ref_curves = []     
    ref_labels = []    
    for arg in args:
        if isinstance(arg, str):
            types.append(arg)
        elif isinstance(arg, dict):
            ref_curves.append(arg)
            ref_labels.append('')
        elif isinstance(arg, tuple) and len(arg) == 2:
            ref_curves.append(arg[0])
            ref_labels.append(arg[1])

    # Get plot mode (SNR or Iterations as parameters)
    if all('-SNR' in t for t in types):
        mode = "-SNR"
    elif all('-Itr' in t for t in types):
        mode = '-Itr'

    ax = plt.subplot()

    ref_type = types[0]
    params = sorted(Sim.error_rates[ref_type].keys())

    # Plot reference curves
    markers = ['x', '+', '*', 's', '^', 'D', 'v', 'p', 'd', '|', '_']
    for ind, ref_curve in enumerate(ref_curves):
        if isinstance(ref_curve, dict):
            ref_params = sorted(ref_curve.keys())
            ref_values = [ref_curve[p] for p in ref_params]

            label = ref_labels[ind]
            marker = markers[ind % len(markers)]

            ax.semilogy(ref_params, ref_values, marker=marker, label=label)
        else:
            print(f"Invalid type of reference curve - must be dict")

    # Plot error rates of simulation
    for type in types:
        if type in Sim.error_rates:
            data = Sim.error_rates[type]
            values = [data[p] for p in params]

            marker = '' if 'preFEC' in type else 'o'
            linestyle = '-' if 'FER' in type else ':'
            label = type.replace(mode, '')

            ax.semilogy(params, values, color='black', marker=marker, linestyle=linestyle, label=label)
        else:
            print(f'No data found for {type}!')

    if mode == '-SNR':
        ax.set_xlabel(f'${Sim.SNR_type}$ (dB)')
    elif mode == '-Itr':
        ax.set_xlabel('Iterations')
        ax.set_xticks(np.linspace(min(params), max(params), 5, dtype=int))

    ax.set_xlim(left=min(params))
    ax.set_ylim(top=1)
    ax.set_ylabel('Error rate')
    ax.legend()
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.show()



def plot_error_rates(*args):
    '''
    Plots error rate curves from dictionaries.

    Parameters:
    ---------------------------------------------------------------
    args    (dict)        - Reference data
            (tuple)       - Tuple containing reference data & label.
    '''
    ax = plt.subplot()

    markers = ['x', '+', '*', 's', '^', 'D', 'v', 'p', 'd', '|', '_']

    for ind, arg in enumerate(args):
        if isinstance(arg, dict):
            ref_curve = arg
            label = ''
        elif isinstance(arg, tuple) and len(arg) == 2:
            ref_curve, label = arg
        else:
            print(f"Invalid argument {arg} - must be dict or (dict, label)")
            continue

        ref_params = sorted(ref_curve.keys())
        ref_values = [ref_curve[p] for p in ref_params]

        marker = markers[ind % len(markers)]
        ax.semilogy(ref_params, ref_values, marker=marker, label=label)

    ax.set_xlim(left=min(ref_params))
    ax.set_ylim(top=1)
    ax.set_ylabel('Error rate')
    ax.legend()
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.show()
