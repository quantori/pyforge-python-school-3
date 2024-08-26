def mark_as_tested(func):
    """
    Decorator to add "TESTED" to the function's docstring.

    Helps me keep track of which functions have been tested and are reliable.

    Because at this point, it's getting hard to test every function in the project.
    """
    if func.__doc__:
        func.__doc__ += "\n\nTESTED"
    else:
        func.__doc__ = "TESTED"
    return func