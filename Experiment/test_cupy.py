import numpy as np
import cupy as cp
import time
from os import environ
environ['OMP_NUM_THREADS'] = '12'


iterations = 1  # Number of iterations

# Create NumPy matrices
a_cpu = np.random.randint(400, size=(1000,50)) 
b_cpu = np.random.randint(-10,10, size=(50,274000)) 

# Create CuPy matrices on the GPU
a_gpu = cp.asarray(a_cpu)
b_gpu = cp.asarray(b_cpu)

# Time CPU calculation
time_cpu = []
for _ in range(iterations):
    start_cpu = time.time()
    c_cpu = np.matmul(a_cpu, b_cpu)
    time_cpu.append(time.time() - start_cpu)

# Time GPU calculation
time_gpu = []
for _ in range(iterations):
    start_gpu = time.time()
    c_gpu = cp.matmul(a_gpu, b_gpu)
    time_gpu.append(time.time() - start_gpu)

# Ensure all GPU calculations have finished
cp.cuda.Stream.null.synchronize()



print(f"CPU time (for {iterations} iterations): {np.mean(time_cpu):.4f} seconds")
print(f"GPU time (for {iterations} iterations): {np.mean(time_gpu):.4f} seconds")
print(f"Speedup: {np.mean(time_cpu) / np.mean(time_gpu):.2f}x")