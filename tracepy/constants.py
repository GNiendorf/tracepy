# Constants used by optimizer
SURV_CONST = 100 # Weight of failed propagation.
MAX_RMS = 999 # Maximum RMS penalty for trace error.

# Constants used for ray objects
MAX_INTERSECTION_ITERATIONS = 1e4 # Max iter before failed intersection search.
MAX_REFRACTION_ITERATIONS = 1e5 # Max iter before failed refraction.
INTERSECTION_CONVERGENCE_TOLERANCE = 1e-6 # Tolerance for intersection search.
REFRACTION_CONVERGENCE_TOLERANCE = 1e-15 # Tolerance for refraction.

# Constants used for plotting
PLOT_ROUNDING_ACC = 14 # Rounding accuracy for spot diagrams and plots.