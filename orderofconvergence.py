import numpy as np
import matplotlib.pyplot as plt



a=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0.058559963267599796, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0.05157971775812126, 0.06009227720622905, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0.05157971775812126, 0.05277287954152942, 0.025194307120084563, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0.05157971775812126, 0.05277287954152942, 0.0563042489366149, 0.09491919773555096, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0.05157971775812126, 0.05277287954152942, 0.0563042489366149, 0.04426096219067894, 0.023324230233886124, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0.05157971775812126, 0.05277287954152942, 0.0563042489366149, 0.04426096219067894, 0.02482125852054806, 0.026838855054359713, 0, 0, 0, 0, 0, 0, 0, 0], [0.05157971775812126, 0.05277287954152942, 0.0563042489366149, 0.04426096219067894, 0.02482125852054806, 0.0247007239828591, 0.014654605597178868, 0, 0, 0, 0, 0, 0, 0], [0.05157971775812126, 0.05277287954152942, 0.0563042489366149, 0.04426096219067894, 0.02482125852054806, 0.0247007239828591, 0.014244909806994742, 0.015151865704734587, 0, 0, 0, 0, 0, 0], [0.05157971775812126, 0.05277287954152942, 0.0563042489366149, 0.04426096219067894, 0.02482125852054806, 0.0247007239828591, 0.014244909806994742, 0.016749479017404426, 0.08490772088113624, 0, 0, 0, 0, 0], [0.05157971775812126, 0.05277287954152942, 0.0563042489366149, 0.04426096219067894, 0.02482125852054806, 0.0247007239828591, 0.014244909806994742, 0.016749479017404426, 0.010449711093808017, 0.003580329721013816, 0, 0, 0, 0], [0.05157971775812126, 0.05277287954152942, 0.0563042489366149, 0.04426096219067894, 0.02482125852054806, 0.0247007239828591, 0.014244909806994742, 0.016749479017404426, 0.010449711093808017, 0.12651226735676774, 0.35156761954661786, 0, 0, 0], [0.05157971775812126, 0.05277287954152942, 0.0563042489366149, 0.04426096219067894, 0.02482125852054806, 0.0247007239828591, 0.014244909806994742, 0.016749479017404426, 0.010449711093808017, 0.12651226735676774, 0.02431743700065669, 0.31989192531005023, 0, 0], [0.05157971775812126, 0.05277287954152942, 0.0563042489366149, 0.04426096219067894, 0.02482125852054806, 0.0247007239828591, 0.014244909806994742, 0.016749479017404426, 0.010449711093808017, 0.12651226735676774, 0.02431743700065669, 0.12154616870319461, 0.01776478738292895, 0], [0.05157971775812126, 0.05277287954152942, 0.0563042489366149, 0.04426096219067894, 0.02482125852054806, 0.0247007239828591, 0.014244909806994742, 0.016749479017404426, 0.010449711093808017, 0.12651226735676774, 0.02431743700065669, 0.12154616870319461, 0.2737512248682149, 0.2193152972031051]]

b=[0.05157971775812126, 0.05277287954152942, 0.0563042489366149, 0.04426096219067894, 0.02482125852054806, 0.0247007239828591, 0.014244909806994742, 0.016749479017404426, 0.010449711093808017, 0.12651226735676774, 0.02431743700065669, 0.12154616870319461, 0.2737512248682149, 0.13457949602380026, 0.02340951519880686]

b =  np.array(b)
a_values =  a
def brusselator(t, y):
    A, B = y
    k1 = k1_value
    k2 = k2_value
    dA_dt = k1 - (k2 + 1) * A + A**2 * B
    dB_dt = k2 * A - A**2 * B
    return np.array([dA_dt, dB_dt])

def rk5_3x5(f, t_span, y0, h, a_values, b):

    t0, tf = t_span
    t = t0
    y = y0
    t_values = [t0]
    y_values = [y0]
    while t < tf:
        k_values = [h * f(t, y)]
        for i in range(1, len(a_values)):
            print(f"Iteration: {i}, a values: {a[i][:i]}")  # Print all a values 
            ti = t + h * sum(a_values[i-1][:i])
            yi = y + sum(a * k for a, k in zip(a[i][:i], k_values))

            k_values.append(h * f(ti, yi))
        y = y + b.dot(k_values)
        t += h
        t_values.append(t)
        y_values.append(y)
    return np.array(t_values), np.array(y_values)



k1_value = 1.0
k2_value = 2.0

# Initial conditions
y0 = np.array([1.5, 3.0])

# Time span
t_span = (1, 50)

# Step size
#h = 0.001

def l2_norm_error(y_numerical, y_reference):
    error = y_numerical - y_reference
    l2_norm = np.linalg.norm(error, axis=1)
    return np.sqrt(np.mean(l2_norm**2))

def calculate_order_of_convergence(errors, hs):
    orders = []
    for i in range(1, len(errors)):
        order = np.log(errors[i] / errors[i-1]) / np.log(hs[i] / hs[i-1])
        orders.append(order)
    return orders

from scipy.interpolate import interp1d

def interpolate_solution(t_values, y_values, t_values_ref):
    # Interpolate numerical solution to match reference solution time points
    y_interpolated = []
    for i in range(y_values.shape[1]):
        f = interp1d(t_values, y_values[:, i], kind='cubic')
        t_max = min(t_values.max(), t_values_ref.max())  # Adjusted end time
        t_new = np.linspace(t_values.min(), t_max, len(t_values_ref))
        y_interpolated.append(f(t_new))
    return np.array(y_interpolated).T

# Call the Runge-Kutta method
#t_values, y_values = rk5_3x5(brusselator, t_span, y0, h, a_values, b)
errors = []

hs = np.linspace(0.01,0.001,20)  # Step sizes to test # Step sizes to test
#hs = np.logspace(np.log10(0.1), np.log10(0.01), 10)
t_values_ref, y_values_ref = rk5_3x5(brusselator, t_span, y0, 0.0001, a_values, b)
for h in hs:
    # Run numerical method
    t_values, y_values = rk5_3x5(brusselator, t_span, y0, h, a_values, b)
    
    # Generate reference solution
    #t_values_ref, y_values_ref = rk5_3x5(brusselator, t_span, y0, 0.001, a_values, b)
    y_values_interpolated = interpolate_solution(t_values, y_values, t_values_ref)
    
    # Calculate L2 norm error
    error = l2_norm_error(y_values_interpolated, y_values_ref)
    errors.append(error)
    
    
orders = calculate_order_of_convergence(errors, hs)

print("Step sizes:", hs)
print("Errors:", errors)
print("Orders of convergence:", orders)

orders = np.diff(np.log(errors)) / np.diff(np.log(hs))



