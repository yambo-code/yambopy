

def add_fig_kwargs(func):
    """
    Decorator that adds keyword arguments for functions returning matplotlib
    figures.

    Taken from pymatgen:
    http://pymatgen.org/
    https://github.com/materialsproject/pymatgen
    """
    def wrapper(*args, **kwargs):
        # pop the kwds used by the decorator.
        title = kwargs.pop("title", None)
        show = kwargs.pop("show", True)
        savefig = kwargs.pop("savefig", None)

        # Call func and return immediately if None is returned.
        fig = func(*args, **kwargs)
        if fig is None:
            return fig

        # Operate on matplotlib figure.
        if title is not None:
            fig.suptitle(title)
        if savefig:
            fig.savefig(savefig)
        if show:
            import matplotlib.pyplot as plt
            plt.show()
        return fig
    return wrapper
