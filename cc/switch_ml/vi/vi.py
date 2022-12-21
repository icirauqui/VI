import numpy as np
import scipy as sp

import seaborn as sns
from matplotlib import pyplot as plt
#%matplotlib inline

# Generate data
N = 100
mu_arr = np.array([1, 10])
sigma_arr = np.array([1, 1])
x = np.append(np.random.normal(mu_arr[0], sigma_arr[0], N), 
              np.random.normal(mu_arr[1], sigma_arr[1], N))
print(x[:10])


fig, ax = plt.subplots(figsize=(12,3))
ax.scatter(x[:N], np.zeros(N), c='green', marker=2, s=150)
ax.scatter(x[N:], np.zeros(N), c='orange', marker='x', s=150)
ax.set_ylim([-1e-4, 1e-4])
_ = ax.set_yticks([])
sns.histplot(x[:N], color='green')
sns.histplot(x[N:], color='orange')