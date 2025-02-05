class SecondOrderRungeKutta:
    def __init__(self, dydt, initial_t, initial_y, h):
        """
        Initializes the solver.

        :param dydt: The differential equation as a Python function f(t, y).
        :param initial_t: Initial value of t.
        :param initial_y: Initial value of y.
        :param h: Step size.
        """
        self.dydt = dydt
        self.t = initial_t
        self.y = initial_y
        self.h = h

    def step(self):
        """
        Advances the solution by one step using the second order Runge-Kutta method.
        """
        k1 = self.h * self.dydt(self.t, self.y)
        k2 = self.h * self.dydt(self.t + 0.5 * self.h, self.y + 0.5 * k1)
        self.y += k2
        self.t += self.h

    def solve(self, final_t):
        """
        Solves the ODE from the initial to the final time.

        :param final_t: The final value of t for the simulation.
        :return: A tuple of lists (ts, ys) with time steps and y values.
        """
        ts = [self.t]
        ys = [self.y]

        while self.t < final_t:
            self.step()
            ts.append(self.t)
            ys.append(self.y)

        return ts, ys
