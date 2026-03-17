import time
from functools import wraps

# decarator to measure number of calls and wall time in seconds of function
"""
@func_profile
def abc(xxx):
    np.log(np.sin(xxx)+1j+np.cos(xxx)+xxx**2)

x = np.linspace(0,10,1000000)

for i in range(10):
    abc(x)

# number of calls to abc
print(abc.call_count)

# total time spent in abc
print(abc.total_time)

"""
def func_profile(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        
        # Increment call count (initialize if first call)
        if not hasattr(wrapper, "call_count"):
            wrapper.call_count = 0
        wrapper.call_count += 1

        result = func(*args, **kwargs)  # Call original function

        # Accumulate time (initialize if first call)
        if not hasattr(wrapper, "total_time"):
            wrapper.total_time = 0.0
        wrapper.total_time += time.perf_counter() - start_time

        return result
    return wrapper
